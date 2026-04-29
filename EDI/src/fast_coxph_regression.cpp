#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>

using namespace Rcpp;

namespace {

struct CoxData {
    const Eigen::VectorXd y;
    const Eigen::VectorXd dead;
    // Store X in row-major order so row[i] is contiguous
    std::vector<double> X_rowmaj;  // n x p, row-major
    std::vector<int> idx_asc;
    std::vector<double> unique_event_times;
    std::vector<int> event_counts;
    int n;
    int p;

    CoxData(const Eigen::VectorXd& y_in, const Eigen::VectorXd& dead_in, const Eigen::MatrixXd& X_in) :
        y(y_in), dead(dead_in), n(y_in.size()), p(X_in.cols()),
        X_rowmaj(y_in.size() * X_in.cols()) {

        // Copy X into row-major layout for cache-friendly row access in inner loop
        for (int i = 0; i < n; ++i)
            for (int q = 0; q < p; ++q)
                X_rowmaj[i * p + q] = X_in(i, q);

        idx_asc.resize(n);
        std::iota(idx_asc.begin(), idx_asc.end(), 0);
        std::sort(idx_asc.begin(), idx_asc.end(), [&](int i, int j) {
            if (y[i] == y[j]) return dead[i] > dead[j];
            return y[i] < y[j];
        });

        for (int i : idx_asc) {
            if (dead[i] > 0.5) {
                if (unique_event_times.empty() || y[i] > unique_event_times.back()) {
                    unique_event_times.push_back(y[i]);
                    event_counts.push_back(0);
                }
            }
        }

        int k = -1;
        for (int i = 0; i < n; ++i) {
            int id = idx_asc[i];
            if (k + 1 < (int)unique_event_times.size() && y[id] == unique_event_times[k+1])
                k++;
            if (k >= 0 && y[id] == unique_event_times[k] && dead[id] > 0.5)
                event_counts[k]++;
        }
    }

    inline const double* row(int i) const { return X_rowmaj.data() + i * p; }
};

// Compute Cox partial log-likelihood, gradient, and Hessian.
// Uses pre-allocated raw C arrays to avoid any Eigen heap allocation in the inner loop.
double compute_cox_ll_grad_hess_fast(
        const CoxData& data,
        const std::vector<double>& beta,
        std::vector<double>& grad,   // length p, reset inside
        std::vector<double>& hess,   // length p*p, reset inside
        bool estimate_only,
        // Pre-allocated work buffers (passed in to avoid re-allocating each Newton step)
        std::vector<double>& exp_eta,   // length n
        std::vector<double>& r_x_exp,   // length p
        std::vector<double>& r_xx_exp,  // length p*p
        std::vector<double>& sum_x_dk,  // length p
        std::vector<double>& e_z        // length p
) {
    const int n = data.n;
    const int p = data.p;

    // --- Compute eta = X * beta and exp_eta ---
    double max_eta = -1e300;
    for (int i = 0; i < n; ++i) {
        const double* xi = data.row(i);
        double eta_i = 0.0;
        for (int q = 0; q < p; ++q) eta_i += xi[q] * beta[q];
        exp_eta[i] = eta_i;  // store eta first; will exp in next pass
        if (eta_i > max_eta) max_eta = eta_i;
    }
    for (int i = 0; i < n; ++i) exp_eta[i] = std::exp(exp_eta[i] - max_eta);

    // --- Initialise running risk-set accumulators (full risk set) ---
    double r_exp = 0.0;
    for (int q = 0; q < p; ++q) r_x_exp[q] = 0.0;
    if (!estimate_only)
        for (int qq = 0; qq < p * p; ++qq) r_xx_exp[qq] = 0.0;

    for (int i = 0; i < n; ++i) {
        const double* xi = data.row(i);
        double wi = exp_eta[i];
        r_exp += wi;
        for (int q = 0; q < p; ++q) r_x_exp[q] += xi[q] * wi;
        if (!estimate_only)
            for (int q = 0; q < p; ++q)
                for (int r = 0; r < p; ++r)
                    r_xx_exp[q * p + r] += xi[q] * xi[r] * wi;
    }

    // --- Reset output grad / hess ---
    for (int q = 0; q < p; ++q) grad[q] = 0.0;
    if (!estimate_only)
        for (int qq = 0; qq < p * p; ++qq) hess[qq] = 0.0;

    double neg_ll = 0.0;
    int j = 0;  // pointer into sorted observations

    for (size_t k = 0; k < data.unique_event_times.size(); ++k) {
        double tk = data.unique_event_times[k];

        // Remove observations whose time < tk from the risk set
        while (j < n && data.y[data.idx_asc[j]] < tk) {
            int id = data.idx_asc[j];
            const double* xi = data.row(id);
            double wi = exp_eta[id];
            r_exp -= wi;
            for (int q = 0; q < p; ++q) r_x_exp[q] -= xi[q] * wi;
            if (!estimate_only)
                for (int q = 0; q < p; ++q)
                    for (int r = 0; r < p; ++r)
                        r_xx_exp[q * p + r] -= xi[q] * xi[r] * wi;
            ++j;
        }

        // Accumulate sum over deaths at this event time
        double sum_eta_dk = 0.0;
        for (int q = 0; q < p; ++q) sum_x_dk[q] = 0.0;
        int count_dk = data.event_counts[k];

        int m = j;
        while (m < n && data.y[data.idx_asc[m]] == tk) {
            int id = data.idx_asc[m];
            if (data.dead[id] > 0.5) {
                const double* xi = data.row(id);
                // eta[id] = log(exp_eta[id]) + max_eta  (recover original eta)
                sum_eta_dk += std::log(std::max(exp_eta[id], 1e-300)) + max_eta;
                for (int q = 0; q < p; ++q) sum_x_dk[q] += xi[q];
            }
            ++m;
        }

        if (count_dk > 0) {
            double safe_r = std::max(r_exp, 1e-100);
            double inv_r = 1.0 / safe_r;
            neg_ll -= (sum_eta_dk - count_dk * (std::log(safe_r) + max_eta));

            // e_z = r_x_exp / r_exp
            for (int q = 0; q < p; ++q) e_z[q] = r_x_exp[q] * inv_r;

            // grad -= (sum_x_dk - count_dk * e_z)
            for (int q = 0; q < p; ++q)
                grad[q] -= sum_x_dk[q] - count_dk * e_z[q];

            if (!estimate_only) {
                // hess += count_dk * (r_xx_exp * inv_r - e_z * e_z^T)
                for (int q = 0; q < p; ++q)
                    for (int r = 0; r < p; ++r)
                        hess[q * p + r] += count_dk * (r_xx_exp[q * p + r] * inv_r - e_z[q] * e_z[r]);
            }
        }
    }

    return neg_ll;
}

// Newton-Raphson for (possibly stratified) Cox model.
// strata_data: vector of CoxData, one per stratum.  For unstratified Cox pass a single element.
// Returns beta (converged), H^{-1} (inverse information), neg_ll, converged, iterations.
struct CoxFitResult {
    std::vector<double> beta;
    Eigen::MatrixXd vcov;   // p x p inverse information (empty if estimate_only)
    Eigen::MatrixXd hess_mat; // p x p information matrix at convergence
    double neg_ll;
    bool converged;
    int iterations;
};

CoxFitResult cox_newton_raphson(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> start_beta,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol)
{
    const int p = strata_data[0].p;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    if (start_beta.isNotNull()) {
        NumericVector sb(start_beta);
        for (int q = 0; q < p; ++q) beta[q] = sb[q];
    }
    beta = apply_fixed_values(beta, fixed_spec);

    std::vector<double> grad_vec(p), hess_vec(p * p);

    double old_ll = 1e300;
    int iter = 0;

    for (iter = 0; iter < maxit; ++iter) {
        // Accumulate over strata
        std::fill(grad_vec.begin(), grad_vec.end(), 0.0);
        std::fill(hess_vec.begin(), hess_vec.end(), 0.0);
        double ll = 0.0;

        std::vector<double> beta_vec(p);
        for (int q = 0; q < p; ++q) beta_vec[q] = beta[q];

        for (const CoxData& sd : strata_data) {
            std::vector<double> g_s(p, 0.0), h_s(p * p, 0.0);
            std::vector<double> exp_eta(sd.n), r_x_exp(p), r_xx_exp(p * p), sum_x_dk(p), e_z(p);
            ll += compute_cox_ll_grad_hess_fast(sd, beta_vec, g_s, h_s, false,
                                                 exp_eta, r_x_exp, r_xx_exp, sum_x_dk, e_z);
            for (int q = 0; q < p; ++q) grad_vec[q] += g_s[q];
            for (int qq = 0; qq < p * p; ++qq) hess_vec[qq] += h_s[qq];
        }

        if (std::abs(old_ll - ll) < tol) break;
        old_ll = ll;

        Eigen::MatrixXd H = Eigen::Map<const Eigen::MatrixXd>(hess_vec.data(), p, p);
        Eigen::VectorXd g = Eigen::Map<const Eigen::VectorXd>(grad_vec.data(), p);

        // Apply fixed parameter constraints to NR step
        for (int i = 0; i < fixed_spec.fixed_idx.size(); ++i) {
            int idx = fixed_spec.fixed_idx[i];
            g[idx] = 0.0;
            H.row(idx).setZero();
            H.col(idx).setZero();
            H(idx, idx) = 1.0;
        }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
        if (ldlt.info() != Eigen::Success) break;
        Eigen::VectorXd delta = ldlt.solve(g);
        beta -= delta;
    }

    CoxFitResult res;
    res.beta.assign(beta.data(), beta.data() + p);
    res.neg_ll = old_ll;
    res.converged = (iter < maxit);
    res.iterations = iter;

    if (!estimate_only) {
        std::fill(hess_vec.begin(), hess_vec.end(), 0.0);
        std::vector<double> beta_vec(p);
        for (int q = 0; q < p; ++q) beta_vec[q] = beta[q];
        for (const CoxData& sd : strata_data) {
            std::vector<double> g_s(p, 0.0), h_s(p * p, 0.0);
            std::vector<double> exp_eta(sd.n), r_x_exp(p), r_xx_exp(p * p), sum_x_dk(p), e_z(p);
            compute_cox_ll_grad_hess_fast(sd, beta_vec, g_s, h_s, false,
                                           exp_eta, r_x_exp, r_xx_exp, sum_x_dk, e_z);
            for (int qq = 0; qq < p * p; ++qq) hess_vec[qq] += h_s[qq];
        }
        Eigen::MatrixXd H_final = Eigen::Map<const Eigen::MatrixXd>(hess_vec.data(), p, p);
        res.hess_mat = H_final;
        Eigen::MatrixXd H_free = subset_matrix(H_final, fixed_spec.free_idx, fixed_spec.free_idx);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(H_free);
        if (lu.isInvertible()) {
            Eigen::MatrixXd vcov_free = lu.inverse();
            res.vcov = expand_free_covariance(p, fixed_spec, vcov_free, true);
        } else {
            res.vcov = Eigen::MatrixXd::Constant(p, p, NA_REAL);
        }
    }

    return res;
}

// ── Functor wrapping stratified Cox partial likelihood for L-BFGS ─────────
class StratifiedCoxObjective {
    const std::vector<CoxData>& m_strata;
    const int m_p;
public:
    StratifiedCoxObjective(const std::vector<CoxData>& strata, int p)
        : m_strata(strata), m_p(p) {}

    // Returns negative partial log-likelihood; fills analytic gradient.
    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        std::vector<double> beta(par.data(), par.data() + m_p);
        std::vector<double> g(m_p, 0.0), h_dummy(m_p * m_p, 0.0);
        double nll = 0.0;
        for (const CoxData& sd : m_strata) {
            std::vector<double> g_s(m_p, 0.0), h_s(m_p * m_p, 0.0);
            std::vector<double> exp_eta(sd.n), r_x_exp(m_p), r_xx_exp(m_p * m_p), sum_x_dk(m_p), e_z(m_p);
            nll += compute_cox_ll_grad_hess_fast(sd, beta, g_s, h_s, /*estimate_only=*/true,
                                                  exp_eta, r_x_exp, r_xx_exp, sum_x_dk, e_z);
            for (int q = 0; q < m_p; ++q) g[q] += g_s[q];
        }
        grad = Eigen::Map<Eigen::VectorXd>(g.data(), m_p);
        return nll;
    }

    // Analytic Hessian (information matrix) — used by Newton-Raphson path.
    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        std::vector<double> beta(par.data(), par.data() + m_p);
        std::vector<double> h(m_p * m_p, 0.0);
        for (const CoxData& sd : m_strata) {
            std::vector<double> g_s(m_p, 0.0), h_s(m_p * m_p, 0.0);
            std::vector<double> exp_eta(sd.n), r_x_exp(m_p), r_xx_exp(m_p * m_p), sum_x_dk(m_p), e_z(m_p);
            compute_cox_ll_grad_hess_fast(sd, beta, g_s, h_s, /*estimate_only=*/false,
                                           exp_eta, r_x_exp, r_xx_exp, sum_x_dk, e_z);
            for (int qq = 0; qq < m_p * m_p; ++qq) h[qq] += h_s[qq];
        }
        return Eigen::Map<Eigen::MatrixXd>(h.data(), m_p, m_p);
    }
};

// Fit via L-BFGS; returns same CoxFitResult as cox_newton_raphson.
CoxFitResult cox_lbfgs(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> start_beta,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol)
{
    const int p = strata_data[0].p;
    Eigen::VectorXd par = Eigen::VectorXd::Zero(p);
    if (start_beta.isNotNull()) {
        NumericVector sb(start_beta);
        for (int q = 0; q < p; ++q) par[q] = sb[q];
    }
    par = apply_fixed_values(par, fixed_spec);

    StratifiedCoxObjective obj(strata_data, p);
    LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, tol, "lbfgs", "lbfgs");

    CoxFitResult res;
    res.beta.assign(fit.params.data(), fit.params.data() + p);
    res.neg_ll    = fit.value;
    res.converged = fit.converged;
    res.iterations = fit.niter;

    if (!estimate_only && fit.converged) {
        Eigen::MatrixXd H = obj.hessian(fit.params);
        res.hess_mat = H;
        Eigen::MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
        Eigen::FullPivLU<Eigen::MatrixXd> lu(H_free);
        if (lu.isInvertible()) {
            Eigen::MatrixXd vcov_free = lu.inverse();
            res.vcov = expand_free_covariance(p, fixed_spec, vcov_free, true);
        } else {
            res.vcov = Eigen::MatrixXd::Constant(p, p, NA_REAL);
        }
    }
    return res;
}

// Dispatch: pick NR or L-BFGS based on optimization_alg string.
CoxFitResult cox_fit(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> start_beta,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol,
    const std::string& optimization_alg)
{
    if (optimization_alg == "lbfgs")
        return cox_lbfgs(strata_data, start_beta, fixed_spec, estimate_only, maxit, tol);
    return cox_newton_raphson(strata_data, start_beta, fixed_spec, estimate_only, maxit, tol);
}

// Compute Lin-Wei-Amato sandwich vcov given converged beta and cluster IDs.
//
// Score residual: U_i[q] = (x_i[q] - e(y_i)[q])*dead_i
//                          - exp(eta_i) * sum_{k: t_k<=y_i} (x_i[q] - e_k[q]) * d_k/R_k
// where e_k = weighted mean of x in risk set at t_k.  sum_i U_i = 0 at MLE.
//
// Sandwich: H^{-1} * (sum_c V_c * V_c^T) * H^{-1},  V_c = sum_{i in c} U_i
Eigen::MatrixXd compute_robust_vcov(
    const std::vector<CoxData>& strata_data,
    const std::vector<double>& beta,
    const Eigen::MatrixXd& H_inv,
    const std::vector<int>& cluster)
{
    int n_total = 0;
    for (const CoxData& sd : strata_data) n_total += sd.n;
    const int p = beta.size();

    Eigen::MatrixXd U(n_total, p);
    U.setZero();

    int row_offset = 0;
    for (const CoxData& sd : strata_data) {
        const int ns = sd.n;

        double max_eta = -1e300;
        std::vector<double> exp_eta(ns);
        for (int i = 0; i < ns; ++i) {
            const double* xi = sd.row(i);
            double eta_i = 0.0;
            for (int q = 0; q < p; ++q) eta_i += xi[q] * beta[q];
            exp_eta[i] = eta_i;
            if (eta_i > max_eta) max_eta = eta_i;
        }
        for (int i = 0; i < ns; ++i) exp_eta[i] = std::exp(exp_eta[i] - max_eta);

        double r_exp = 0.0;
        std::vector<double> r_x_exp(p, 0.0);
        for (int i = 0; i < ns; ++i) {
            const double* xi = sd.row(i);
            r_exp += exp_eta[i];
            for (int q = 0; q < p; ++q) r_x_exp[q] += xi[q] * exp_eta[i];
        }

        const int n_events = (int)sd.unique_event_times.size();
        // dk/Rk, e_k (weighted mean x at t_k), and dk*e_k/Rk per event time
        std::vector<double> dk_over_Rk(n_events, 0.0);
        std::vector<std::vector<double>> ek(n_events, std::vector<double>(p, 0.0));
        std::vector<std::vector<double>> dk_ek_over_Rk(n_events, std::vector<double>(p, 0.0));

        int j = 0;
        for (int k = 0; k < n_events; ++k) {
            double tk = sd.unique_event_times[k];
            while (j < ns && sd.y[sd.idx_asc[j]] < tk) {
                int id = sd.idx_asc[j];
                const double* xi = sd.row(id);
                r_exp -= exp_eta[id];
                for (int q = 0; q < p; ++q) r_x_exp[q] -= xi[q] * exp_eta[id];
                ++j;
            }
            int dk = sd.event_counts[k];
            if (dk > 0) {
                double safe_r = std::max(r_exp, 1e-100);
                dk_over_Rk[k] = (double)dk / safe_r;
                for (int q = 0; q < p; ++q) {
                    ek[k][q] = r_x_exp[q] / safe_r;
                    dk_ek_over_Rk[k][q] = (double)dk * ek[k][q] / safe_r;
                }
            }
        }

        // Cumulative A(t) = sum_{k: t_k<=t} d_k/R_k  and  B_q(t) = sum d_k*e_k[q]/R_k
        std::vector<double> cum_A(n_events, 0.0);
        std::vector<std::vector<double>> cum_B(n_events, std::vector<double>(p, 0.0));
        if (n_events > 0) {
            cum_A[0] = dk_over_Rk[0];
            for (int q = 0; q < p; ++q) cum_B[0][q] = dk_ek_over_Rk[0][q];
            for (int k = 1; k < n_events; ++k) {
                cum_A[k] = cum_A[k-1] + dk_over_Rk[k];
                for (int q = 0; q < p; ++q)
                    cum_B[k][q] = cum_B[k-1][q] + dk_ek_over_Rk[k][q];
            }
        }

        // Score residuals: U_i[q] = (x_i[q] - e(y_i)[q])*dead_i
        //                           - exp_i * (x_i[q]*A(y_i) - B_q(y_i))
        for (int i = 0; i < ns; ++i) {
            double yi = sd.y[i];
            const double* xi = sd.row(i);
            double ei = exp_eta[i];
            double di = sd.dead[i];

            // Binary search: last event time <= y_i
            int k_last = -1;
            if (n_events > 0) {
                int lo = 0, hi = n_events - 1;
                while (lo <= hi) {
                    int mid = (lo + hi) / 2;
                    if (sd.unique_event_times[mid] <= yi) { k_last = mid; lo = mid + 1; }
                    else hi = mid - 1;
                }
            }
            double A_i = (k_last >= 0) ? cum_A[k_last] : 0.0;

            // e(y_i) for dead subjects: the event time y_i equals t_{k_last} when dead_i=1
            // (since y_i is one of the unique event times in this case)
            for (int q = 0; q < p; ++q) {
                double B_iq    = (k_last >= 0) ? cum_B[k_last][q] : 0.0;
                double e_yi_q  = (di > 0.5 && k_last >= 0) ? ek[k_last][q] : 0.0;
                U(row_offset + i, q) = (xi[q] - e_yi_q) * di
                                       - ei * (xi[q] * A_i - B_iq);
            }
        }
        row_offset += ns;
    }

    // Aggregate score residuals by cluster: V_c = sum_{i in c} U_i
    // cluster is length n_total, one entry per row of U (in strata-stacked order)
    std::map<int, Eigen::VectorXd> cluster_scores;
    for (int i = 0; i < n_total; ++i) {
        int c = cluster[i];
        auto it = cluster_scores.find(c);
        if (it == cluster_scores.end()) {
            cluster_scores[c] = U.row(i).transpose();
        } else {
            it->second += U.row(i).transpose();
        }
    }

    // B = sum_c V_c * V_c^T
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(p, p);
    for (auto& kv : cluster_scores) {
        B += kv.second * kv.second.transpose();
    }

    // Robust vcov = H^{-1} * B * H^{-1}
    return H_inv * B * H_inv;
}

} // namespace

//' @title Fast Cox Proportional Hazards Regression (C++)
//' @description High-performance Cox regression fitting using Breslow's method for ties.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators (1=event, 0=censored).
//' @param X A numeric matrix of predictors.
//' @param start_beta Optional starting values for coefficients.
//' @param estimate_only If TRUE, only return coefficients and likelihood.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param cluster Optional vector of cluster IDs for robust variance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm ("newton_raphson" or "lbfgs").
//' @return A list containing coefficients, vcov (optional), and convergence status.
//' @export
// [[Rcpp::export]]
List fast_coxph_regression_cpp(const Eigen::VectorXd& y,
                               const Eigen::VectorXd& dead,
                               const Eigen::MatrixXd& X,
                               Nullable<NumericVector> start_beta = R_NilValue,
                               bool estimate_only = false,
                               int maxit = 20,
                               double tol = 1e-9,
                               Nullable<IntegerVector> cluster = R_NilValue,
                               Nullable<IntegerVector> fixed_idx = R_NilValue,
                               Nullable<NumericVector> fixed_values = R_NilValue,
                               std::string optimization_alg = "newton_raphson") {
    int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    std::vector<CoxData> strata_data;
    strata_data.emplace_back(y, dead, X);

    CoxFitResult fit = cox_fit(strata_data, start_beta, fixed_spec, estimate_only, maxit, tol, optimization_alg);

    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = fit.beta[q];

    if (estimate_only) {
        return List::create(
            Named("coefficients") = coef_r,
            Named("converged") = fit.converged,
            Named("neg_ll") = fit.neg_ll,
            Named("iterations") = fit.iterations
        );
    }

    // Use robust sandwich vcov if cluster provided
    Eigen::MatrixXd vcov_mat;
    if (cluster.isNotNull()) {
        IntegerVector cl(cluster);
        std::vector<int> cl_vec(cl.begin(), cl.end());
        vcov_mat = compute_robust_vcov(strata_data, fit.beta, fit.vcov, cl_vec);
    } else {
        vcov_mat = fit.vcov;
    }

    return List::create(
        Named("coefficients") = coef_r,
        Named("vcov") = vcov_mat,
        Named("converged") = fit.converged,
        Named("neg_ll") = fit.neg_ll,
        Named("iterations") = fit.iterations
    );
}

//' @title Fast Stratified Cox Proportional Hazards Regression (C++)
//' @description High-performance stratified Cox regression fitting.
//' @param y A numeric vector of survival times.
//' @param dead A numeric vector of event indicators.
//' @param X A numeric matrix of predictors.
//' @param strata_r An integer vector of strata IDs.
//' @param start_beta Optional starting values for coefficients.
//' @param estimate_only If TRUE, only return coefficients and likelihood.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, vcov, and convergence status.
//' @export
// [[Rcpp::export]]
List fast_stratified_coxph_regression_cpp(
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& dead,
    const Eigen::MatrixXd& X,
    const Rcpp::IntegerVector& strata_r,
    Nullable<NumericVector> start_beta = R_NilValue,
    bool estimate_only = false,
    int maxit = 20,
    double tol = 1e-9,
    Nullable<IntegerVector> fixed_idx = R_NilValue,
    Nullable<NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "newton_raphson")
{
    const int n = y.size();
    const int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    std::vector<int> strata(strata_r.begin(), strata_r.end());
    std::vector<int> unique_strata = strata;
    std::sort(unique_strata.begin(), unique_strata.end());
    unique_strata.erase(std::unique(unique_strata.begin(), unique_strata.end()), unique_strata.end());

    // Build one CoxData per stratum
    std::vector<CoxData> strata_data;
    strata_data.reserve(unique_strata.size());
    for (int s : unique_strata) {
        std::vector<int> idx;
        for (int i = 0; i < n; ++i) if (strata[i] == s) idx.push_back(i);

        const int ns = (int)idx.size();
        Eigen::VectorXd y_s(ns), dead_s(ns);
        Eigen::MatrixXd X_s(ns, p);
        for (int ii = 0; ii < ns; ++ii) {
            y_s[ii]    = y[idx[ii]];
            dead_s[ii] = dead[idx[ii]];
            X_s.row(ii) = X.row(idx[ii]);
        }
        strata_data.emplace_back(y_s, dead_s, X_s);
    }

    CoxFitResult fit = cox_fit(strata_data, start_beta, fixed_spec, estimate_only, maxit, tol, optimization_alg);

    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = fit.beta[q];

    if (estimate_only) {
        return List::create(
            Named("coefficients") = coef_r,
            Named("converged") = fit.converged,
            Named("neg_ll") = fit.neg_ll,
            Named("iterations") = fit.iterations
        );
    }

    return List::create(
        Named("coefficients") = coef_r,
        Named("vcov") = fit.vcov,
        Named("converged") = fit.converged,
        Named("neg_ll") = fit.neg_ll,
        Named("iterations") = fit.iterations
    );
}

//' @title Compute Cox PH Score (C++)
//' @description Calculates the score vector (gradient of log partial likelihood) for unstratified Cox regression.
//' @param y Survival times.
//' @param dead Event indicators.
//' @param X Predictor matrix.
//' @param beta Coefficient vector.
//' @return Score vector.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd get_coxph_score_cpp(const Eigen::VectorXd& y,
                                     const Eigen::VectorXd& dead,
                                     const Eigen::MatrixXd& X,
                                     const Eigen::VectorXd& beta) {
    std::vector<CoxData> strata_data;
    strata_data.emplace_back(y, dead, X);
    StratifiedCoxObjective obj(strata_data, X.cols());
    Eigen::VectorXd grad(beta.size());
    obj(beta, grad);
    return -grad;
}

//' @title Compute Cox PH Hessian (C++)
//' @description Calculates the Hessian of log partial likelihood for unstratified Cox regression.
//' @param y Survival times.
//' @param dead Event indicators.
//' @param X Predictor matrix.
//' @param beta Coefficient vector.
//' @return Hessian matrix.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_coxph_hessian_cpp(const Eigen::VectorXd& y,
                                       const Eigen::VectorXd& dead,
                                       const Eigen::MatrixXd& X,
                                       const Eigen::VectorXd& beta) {
    std::vector<CoxData> strata_data;
    strata_data.emplace_back(y, dead, X);
    StratifiedCoxObjective obj(strata_data, X.cols());
    return -obj.hessian(beta);
}

//' @title Compute Stratified Cox PH Score (C++)
//' @description Calculates the score vector (gradient of log partial likelihood) for stratified Cox regression.
//' @param y Survival times.
//' @param dead Event indicators.
//' @param X Predictor matrix.
//' @param strata_r Integer vector of strata IDs.
//' @param beta Coefficient vector.
//' @return Score vector.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd get_stratified_coxph_score_cpp(const Eigen::VectorXd& y,
                                                const Eigen::VectorXd& dead,
                                                const Eigen::MatrixXd& X,
                                                const Rcpp::IntegerVector& strata_r,
                                                const Eigen::VectorXd& beta) {
    const int n = y.size();
    const int p = X.cols();
    std::vector<int> strata(strata_r.begin(), strata_r.end());
    std::vector<int> unique_strata = strata;
    std::sort(unique_strata.begin(), unique_strata.end());
    unique_strata.erase(std::unique(unique_strata.begin(), unique_strata.end()), unique_strata.end());

    std::vector<CoxData> strata_data;
    strata_data.reserve(unique_strata.size());
    for (int s : unique_strata) {
        std::vector<int> idx;
        for (int i = 0; i < n; ++i) if (strata[i] == s) idx.push_back(i);
        const int ns = (int)idx.size();
        Eigen::VectorXd y_s(ns), dead_s(ns);
        Eigen::MatrixXd X_s(ns, p);
        for (int ii = 0; ii < ns; ++ii) {
            y_s[ii] = y[idx[ii]]; dead_s[ii] = dead[idx[ii]];
            X_s.row(ii) = X.row(idx[ii]);
        }
        strata_data.emplace_back(y_s, dead_s, X_s);
    }
    StratifiedCoxObjective obj(strata_data, p);
    Eigen::VectorXd grad(beta.size());
    obj(beta, grad);
    return -grad;
}

//' @title Compute Stratified Cox PH Hessian (C++)
//' @description Calculates the Hessian of log partial likelihood for stratified Cox regression.
//' @param y Survival times.
//' @param dead Event indicators.
//' @param X Predictor matrix.
//' @param strata_r Integer vector of strata IDs.
//' @param beta Coefficient vector.
//' @return Hessian matrix.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_stratified_coxph_hessian_cpp(const Eigen::VectorXd& y,
                                                   const Eigen::VectorXd& dead,
                                                   const Eigen::MatrixXd& X,
                                                   const Rcpp::IntegerVector& strata_r,
                                                   const Eigen::VectorXd& beta) {
    const int n = y.size();
    const int p = X.cols();
    std::vector<int> strata(strata_r.begin(), strata_r.end());
    std::vector<int> unique_strata = strata;
    std::sort(unique_strata.begin(), unique_strata.end());
    unique_strata.erase(std::unique(unique_strata.begin(), unique_strata.end()), unique_strata.end());

    std::vector<CoxData> strata_data;
    strata_data.reserve(unique_strata.size());
    for (int s : unique_strata) {
        std::vector<int> idx;
        for (int i = 0; i < n; ++i) if (strata[i] == s) idx.push_back(i);
        const int ns = (int)idx.size();
        Eigen::VectorXd y_s(ns), dead_s(ns);
        Eigen::MatrixXd X_s(ns, p);
        for (int ii = 0; ii < ns; ++ii) {
            y_s[ii] = y[idx[ii]]; dead_s[ii] = dead[idx[ii]];
            X_s.row(ii) = X.row(idx[ii]);
        }
        strata_data.emplace_back(y_s, dead_s, X_s);
    }
    StratifiedCoxObjective obj(strata_data, p);
    return -obj.hessian(beta);
}
