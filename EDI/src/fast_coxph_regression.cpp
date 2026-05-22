#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>

using namespace Rcpp;

namespace {

using RowMajorMatrixMap =
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

struct CoxData {
    Eigen::VectorXd y;
    Eigen::VectorXd dead;
    std::vector<double> X_rowmaj;
    std::vector<int> idx_asc;
    std::vector<double> unique_event_times;
    std::vector<int> event_counts;
    int n;
    int p;

    CoxData(const Eigen::VectorXd& y_in, const Eigen::VectorXd& dead_in, const Eigen::MatrixXd& X_in) :
        y(y_in), dead(dead_in), n(y_in.size()), p(X_in.cols()),
        X_rowmaj(y_in.size() * X_in.cols()) {

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
    inline RowMajorMatrixMap matrix_map() const { return RowMajorMatrixMap(X_rowmaj.data(), n, p); }
};

struct CoxWorkspace {
    Eigen::VectorXd eta;
    std::vector<double> exp_eta;
    std::vector<double> r_x_exp;
    std::vector<double> r_xx_exp;
    std::vector<double> sum_x_dk;
    std::vector<double> e_z;
    std::vector<double> grad;
    std::vector<double> hess;

    CoxWorkspace() = default;

    CoxWorkspace(int n, int p) :
        eta(n),
        exp_eta(n),
        r_x_exp(p),
        r_xx_exp(p * p),
        sum_x_dk(p),
        e_z(p),
        grad(p),
        hess(p * p) {}
};

double compute_cox_neg_ll_only(const CoxData& data, const std::vector<double>& beta, CoxWorkspace& workspace) {
    const int n = data.n;
    const int p = data.p;
    Eigen::Map<const Eigen::VectorXd> beta_map(beta.data(), p);
    workspace.eta.noalias() = data.matrix_map() * beta_map;
    const double max_eta = workspace.eta.maxCoeff();
    Eigen::Map<Eigen::VectorXd>(workspace.exp_eta.data(), n) =
        (workspace.eta.array() - max_eta).exp().matrix();
    double r_exp = Eigen::Map<const Eigen::VectorXd>(workspace.exp_eta.data(), n).sum();
    double neg_ll = 0.0;
    int j = 0;
    for (size_t k = 0; k < data.unique_event_times.size(); ++k) {
        const double tk = data.unique_event_times[k];
        while (j < n && data.y[data.idx_asc[j]] < tk) {
            r_exp -= workspace.exp_eta[data.idx_asc[j]];
            ++j;
        }
        double sum_eta_dk = 0.0;
        const int count_dk = data.event_counts[k];
        for (int m = j; m < n && data.y[data.idx_asc[m]] == tk; ++m) {
            int id = data.idx_asc[m];
            if (data.dead[id] > 0.5) sum_eta_dk += workspace.eta[id];
        }
        if (count_dk > 0) {
            const double safe_r = std::max(r_exp, 1e-100);
            neg_ll -= (sum_eta_dk - count_dk * (std::log(safe_r) + max_eta));
        }
    }
    return neg_ll;
}

std::vector<CoxWorkspace> make_cox_workspaces(const std::vector<CoxData>& strata_data) {
    std::vector<CoxWorkspace> workspaces;
    workspaces.reserve(strata_data.size());
    for (const CoxData& sd : strata_data) {
        workspaces.emplace_back(sd.n, sd.p);
    }
    return workspaces;
};

double compute_cox_ll_grad_hess_fast(
        const CoxData& data,
        const std::vector<double>& beta,
        std::vector<double>& grad,
        std::vector<double>& hess,
        bool estimate_only,
        CoxWorkspace& workspace
) {
    const int n = data.n;
    const int p = data.p;
    std::vector<double>& exp_eta = workspace.exp_eta;
    std::vector<double>& r_x_exp = workspace.r_x_exp;
    std::vector<double>& r_xx_exp = workspace.r_xx_exp;
    std::vector<double>& sum_x_dk = workspace.sum_x_dk;
    std::vector<double>& e_z = workspace.e_z;
    Eigen::Map<Eigen::VectorXd> r_x_exp_map(r_x_exp.data(), p);
    Eigen::Map<Eigen::MatrixXd> r_xx_exp_map(r_xx_exp.data(), p, p);
    Eigen::Map<Eigen::VectorXd> sum_x_dk_map(sum_x_dk.data(), p);
    Eigen::Map<Eigen::VectorXd> e_z_map(e_z.data(), p);
    Eigen::Map<Eigen::VectorXd> grad_map(grad.data(), p);
    Eigen::Map<Eigen::MatrixXd> hess_map(hess.data(), p, p);

    Eigen::Map<const Eigen::VectorXd> beta_map(beta.data(), p);
    workspace.eta.noalias() = data.matrix_map() * beta_map;
    const double max_eta = workspace.eta.maxCoeff();
    Eigen::Map<Eigen::VectorXd>(exp_eta.data(), n) =
        (workspace.eta.array() - max_eta).exp().matrix();

    Eigen::Map<const Eigen::VectorXd> exp_eta_map(exp_eta.data(), n);
    double r_exp = exp_eta_map.sum();
    r_x_exp_map.noalias() = data.matrix_map().transpose() * exp_eta_map;
    if (!estimate_only) {
        r_xx_exp_map.noalias() = weighted_crossprod(data.matrix_map(), exp_eta_map);
    }

    grad_map.setZero();
    if (!estimate_only) hess_map.setZero();

    double neg_ll = 0.0;
    int j = 0;

    for (size_t k = 0; k < data.unique_event_times.size(); ++k) {
        double tk = data.unique_event_times[k];

        while (j < n && data.y[data.idx_asc[j]] < tk) {
            int id = data.idx_asc[j];
            Eigen::Map<const Eigen::VectorXd> xi(data.row(id), p);
            const double wi = exp_eta[id];
            r_exp -= wi;
            r_x_exp_map.noalias() -= wi * xi;
            if (!estimate_only) r_xx_exp_map.noalias() -= wi * (xi * xi.transpose());
            ++j;
        }

        double sum_eta_dk = 0.0;
        sum_x_dk_map.setZero();
        int count_dk = data.event_counts[k];

        int m = j;
        while (m < n && data.y[data.idx_asc[m]] == tk) {
            int id = data.idx_asc[m];
            if (data.dead[id] > 0.5) {
                Eigen::Map<const Eigen::VectorXd> xi(data.row(id), p);
                sum_eta_dk += workspace.eta[id];
                sum_x_dk_map.noalias() += xi;
            }
            ++m;
        }

        if (count_dk > 0) {
            double safe_r = std::max(r_exp, 1e-100);
            double inv_r = 1.0 / safe_r;
            neg_ll -= (sum_eta_dk - count_dk * (std::log(safe_r) + max_eta));
            e_z_map.noalias() = r_x_exp_map * inv_r;
            grad_map.noalias() -= sum_x_dk_map - count_dk * e_z_map;
            if (!estimate_only) {
                hess_map.noalias() += count_dk * (r_xx_exp_map * inv_r - e_z_map * e_z_map.transpose());
            }
        }
    }

    return neg_ll;
}

struct CoxFitResult {
    std::vector<double> beta;
    Eigen::MatrixXd vcov;
    Eigen::MatrixXd hess_mat;
    double neg_ll;
    bool converged;
    int iterations;
};

CoxFitResult cox_newton_raphson(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> warm_start_beta,
    bool smart_cold_start,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    const int p = strata_data[0].p;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        NumericVector sb(warm_start_beta);
        for (int q = 0; q < p; ++q) beta[q] = sb[q];
    } else if (smart_cold_start) {
        int total_n = 0;
        for (const auto& sd : strata_data) total_n += sd.n;
        Eigen::MatrixXd X_full(total_n, p);
        Eigen::VectorXd log_y(total_n);
        int offset = 0;
        for (const auto& sd : strata_data) {
            X_full.block(offset, 0, sd.n, p) = sd.matrix_map();
            for (int i = 0; i < sd.n; ++i) log_y[offset + i] = std::log(std::max(sd.y[i], 1e-8));
            offset += sd.n;
        }
        beta = -safe_ols_solve(X_full, log_y);
    }
    beta = apply_fixed_values(beta, fixed_spec);

    std::vector<double> grad_vec(p), hess_vec(p * p);
    std::vector<CoxWorkspace> workspaces = make_cox_workspaces(strata_data);
    Eigen::VectorXd beta_candidate(p);
    std::vector<double> beta_cand_vec(p);

    double old_ll = 1e300;
    int iter = 0;

    for (iter = 0; iter < maxit; ++iter) {
        std::fill(grad_vec.begin(), grad_vec.end(), 0.0);
        std::fill(hess_vec.begin(), hess_vec.end(), 0.0);
        double ll = 0.0;

        std::vector<double> beta_vec(p);
        for (int q = 0; q < p; ++q) beta_vec[q] = beta[q];

        for (std::size_t s = 0; s < strata_data.size(); ++s) {
            const CoxData& sd = strata_data[s];
            CoxWorkspace& ws = workspaces[s];
            // estimate_only skips vcov materialization, but Newton-Raphson still
            // needs the Hessian to update the coefficients correctly.
            ll += compute_cox_ll_grad_hess_fast(sd, beta_vec, ws.grad, ws.hess, false, ws);
            for (int q = 0; q < p; ++q) grad_vec[q] += ws.grad[q];
            if (iter > 0 || warm_start_fisher_info.isNull()) {
                for (int qq = 0; qq < p * p; ++qq) hess_vec[qq] += ws.hess[qq];
            }
        }

        if (std::abs(old_ll - ll) < tol) break;
        old_ll = ll;

        Eigen::MatrixXd H;
        if (iter == 0 && warm_start_fisher_info.isNotNull()) {
            H = as<Eigen::MatrixXd>(warm_start_fisher_info);
        } else {
            H = Eigen::Map<const Eigen::MatrixXd>(hess_vec.data(), p, p);
        }
        Eigen::VectorXd g = Eigen::Map<const Eigen::VectorXd>(grad_vec.data(), p);

        for (int i = 0; i < (int)fixed_spec.fixed_idx.size(); ++i) {
            int idx = fixed_spec.fixed_idx[i];
            g[idx] = 0.0;
            H.row(idx).setZero();
            H.col(idx).setZero();
            H(idx, idx) = 1.0;
        }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
        if (ldlt.info() != Eigen::Success) break;
        const Eigen::VectorXd delta = ldlt.solve(g);

        // Step-halving line search: ensure the neg log-likelihood decreases.
        // Each halving costs ~1/3 of a full iteration (eta + exp + event scan,
        // no gradient or Hessian), so up to 10 halvings is cheap.
        static const int kMaxHalvings = 10;
        auto eval_ll_candidate = [&](double step) -> double {
            for (int q = 0; q < p; ++q) beta_candidate[q] = beta[q] - step * delta[q];
            for (int fi = 0; fi < (int)fixed_spec.fixed_idx.size(); ++fi)
                beta_candidate[fixed_spec.fixed_idx[fi]] = fixed_spec.fixed_values[fi];
            for (int q = 0; q < p; ++q) beta_cand_vec[q] = beta_candidate[q];
            double ll_c = 0.0;
            for (std::size_t s = 0; s < strata_data.size(); ++s)
                ll_c += compute_cox_neg_ll_only(strata_data[s], beta_cand_vec, workspaces[s]);
            return ll_c;
        };

        double step = 1.0;
        double ll_candidate = eval_ll_candidate(step);
        for (int h = 0; h < kMaxHalvings && ll_candidate >= ll; ++h) {
            step *= 0.5;
            ll_candidate = eval_ll_candidate(step);
        }
        beta = beta_candidate;
    }

    CoxFitResult res;
    res.beta.assign(beta.data(), beta.data() + p);
    res.neg_ll = old_ll;
    res.converged = (iter < maxit);
    res.iterations = iter;
    res.hess_mat = Eigen::Map<const Eigen::MatrixXd>(hess_vec.data(), p, p);

    if (!estimate_only) {
        Eigen::MatrixXd H_free = subset_matrix(res.hess_mat, fixed_spec.free_idx, fixed_spec.free_idx);
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

class StratifiedCoxObjective {
    const std::vector<CoxData>& m_strata;
    const int m_p;
    mutable std::vector<CoxWorkspace> m_workspaces;
public:
    StratifiedCoxObjective(const std::vector<CoxData>& strata, int p)
        : m_strata(strata), m_p(p), m_workspaces(make_cox_workspaces(strata)) {}

    double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
        std::vector<double> beta(par.data(), par.data() + m_p);
        std::vector<double> g(m_p, 0.0);
        double nll = 0.0;
        for (std::size_t s = 0; s < m_strata.size(); ++s) {
            const CoxData& sd = m_strata[s];
            CoxWorkspace& ws = m_workspaces[s];
            nll += compute_cox_ll_grad_hess_fast(sd, beta, ws.grad, ws.hess, true, ws);
            for (int q = 0; q < m_p; ++q) g[q] += ws.grad[q];
        }
        grad = Eigen::Map<Eigen::VectorXd>(g.data(), m_p);
        return nll;
    }

    Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
        std::vector<double> beta(par.data(), par.data() + m_p);
        std::vector<double> h(m_p * m_p, 0.0);
        for (std::size_t s = 0; s < m_strata.size(); ++s) {
            const CoxData& sd = m_strata[s];
            CoxWorkspace& ws = m_workspaces[s];
            compute_cox_ll_grad_hess_fast(sd, beta, ws.grad, ws.hess, false, ws);
            for (int qq = 0; qq < m_p * m_p; ++qq) h[qq] += ws.hess[qq];
        }
        return Eigen::Map<Eigen::MatrixXd>(h.data(), m_p, m_p);
    }
};

CoxFitResult cox_lbfgs(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> warm_start_beta,
    bool smart_cold_start,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    const int p = strata_data[0].p;
    Eigen::VectorXd par = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        NumericVector sb(warm_start_beta);
        for (int q = 0; q < p; ++q) par[q] = sb[q];
    } else if (smart_cold_start) {
        int total_n = 0;
        for (const auto& sd : strata_data) total_n += sd.n;
        Eigen::MatrixXd X_full(total_n, p);
        Eigen::VectorXd log_y(total_n);
        int offset = 0;
        for (const auto& sd : strata_data) {
            X_full.block(offset, 0, sd.n, p) = sd.matrix_map();
            for (int i = 0; i < sd.n; ++i) log_y[offset + i] = std::log(std::max(sd.y[i], 1e-8));
            offset += sd.n;
        }
        par = -safe_ols_solve(X_full, log_y);
    }
    par = apply_fixed_values(par, fixed_spec);

    StratifiedCoxObjective obj(strata_data, p);
    Eigen::MatrixXd info_start;
    Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, tol, "lbfgs", "lbfgs", 0, info_start_ptr);

    CoxFitResult res;
    res.beta.assign(fit.params.data(), fit.params.data() + p);
    res.neg_ll    = fit.value;
    res.converged = fit.converged;
    res.iterations = fit.niter;

    if (!estimate_only && fit.converged) {
        res.hess_mat = obj.hessian(fit.params);
        Eigen::MatrixXd H_free = subset_matrix(res.hess_mat, fixed_spec.free_idx, fixed_spec.free_idx);
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

CoxFitResult cox_fit(
    const std::vector<CoxData>& strata_data,
    Nullable<NumericVector> warm_start_beta,
    bool smart_cold_start,
    const FixedParamSpec& fixed_spec,
    bool estimate_only,
    int maxit,
    double tol,
    const std::string& optimization_alg,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    if (optimization_alg == "lbfgs")
        return cox_lbfgs(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol);
    return cox_newton_raphson(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, warm_start_fisher_info);
}

Eigen::MatrixXd compute_robust_vcov(
    const std::vector<CoxData>& strata_data,
    const std::vector<double>& beta,
    const Eigen::MatrixXd& H_inv,
    const std::vector<int>& cluster)
{
    int n_total = 0;
    for (const CoxData& sd : strata_data) n_total += sd.n;
    const int p = beta.size();
    Eigen::Map<const Eigen::VectorXd> beta_map(beta.data(), p);

    Eigen::MatrixXd U(n_total, p);
    U.setZero();

    int row_offset = 0;
    for (const CoxData& sd : strata_data) {
        const int ns = sd.n;
        Eigen::VectorXd eta = sd.matrix_map() * beta_map;
        std::vector<double> exp_eta(ns);
        Eigen::Map<Eigen::VectorXd>(exp_eta.data(), ns) =
            (eta.array() - eta.maxCoeff()).exp().matrix();

        double r_exp = 0.0;
        std::vector<double> r_x_exp(p, 0.0);
        for (int i = 0; i < ns; ++i) {
            const double* xi = sd.row(i);
            r_exp += exp_eta[i];
            for (int q = 0; q < p; ++q) r_x_exp[q] += xi[q] * exp_eta[i];
        }

        const int n_events = (int)sd.unique_event_times.size();
        std::vector<double> dk_over_Rk(n_events, 0.0);
        std::vector<std::vector<double>> ek(n_events, std::vector<double>(p, 0.0));
        std::vector<std::vector<double>> dk_ek_over_Rk(n_events, std::vector<double>(p, 0.0));

        int j = 0;
        for (int k = 0; k < n_events; ++k) {
            double tk = sd.unique_event_times[k];
            while (j < ns && sd.y[sd.idx_asc[j]] < tk) {
                int id = sd.idx_asc[j];
                r_exp -= exp_eta[id];
                for (int q = 0; q < p; ++q) r_x_exp[q] -= sd.row(id)[q] * exp_eta[id];
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

        for (int i = 0; i < ns; ++i) {
            double yi = sd.y[i];
            const double* xi = sd.row(i);
            double ei = exp_eta[i];
            double di = sd.dead[i];

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
            for (int q = 0; q < p; ++q) {
                double B_iq    = (k_last >= 0) ? cum_B[k_last][q] : 0.0;
                double e_yi_q  = (di > 0.5 && k_last >= 0) ? ek[k_last][q] : 0.0;
                U(row_offset + i, q) = (xi[q] - e_yi_q) * di - ei * (xi[q] * A_i - B_iq);
            }
        }
        row_offset += ns;
    }

    std::map<int, Eigen::VectorXd> cluster_scores;
    for (int i = 0; i < n_total; ++i) {
        int c = cluster[i];
        if (cluster_scores.find(c) == cluster_scores.end()) cluster_scores[c] = U.row(i).transpose();
        else cluster_scores[c] += U.row(i).transpose();
    }

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(p, p);
    for (auto& kv : cluster_scores) B += kv.second * kv.second.transpose();
    return H_inv * B * H_inv;
}

} // namespace

// [[Rcpp::export]]
SEXP build_cox_data_cache_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead) {
    auto* data = new std::vector<CoxData>();
    data->emplace_back(y, dead, X);
    return Rcpp::XPtr<std::vector<CoxData>>(data, true);
}

// [[Rcpp::export]]
List fast_coxph_regression_prebuilt_cpp(
    SEXP cox_data_xptr,
    Nullable<NumericVector> warm_start_beta = R_NilValue,
    bool smart_cold_start = true,
    bool estimate_only = false,
    int maxit = 20,
    double tol = 1e-9,
    Nullable<IntegerVector> fixed_idx = R_NilValue,
    Nullable<NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "newton_raphson",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    Rcpp::XPtr<std::vector<CoxData>> data_ptr(cox_data_xptr);
    const std::vector<CoxData>& strata_data = *data_ptr;
    const int p = strata_data[0].p;
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    CoxFitResult fit = cox_fit(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, optimization_alg, warm_start_fisher_info);
    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = fit.beta[q];
    if (estimate_only) {
        return List::create(_["coefficients"] = coef_r, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
    }
    return List::create(_["coefficients"] = coef_r, _["vcov"] = fit.vcov, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
}

// [[Rcpp::export]]
List fast_coxph_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead,
                               Nullable<NumericVector> warm_start_beta = R_NilValue,
                               bool smart_cold_start = true,
                               bool estimate_only = false,
                               int maxit = 20,
                               double tol = 1e-9,
                               Nullable<IntegerVector> cluster = R_NilValue,
                               Nullable<IntegerVector> fixed_idx = R_NilValue,
                               Nullable<NumericVector> fixed_values = R_NilValue,
                               std::string optimization_alg = "newton_raphson",
                               Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    std::vector<CoxData> strata_data;
    strata_data.emplace_back(y, dead, X);
    CoxFitResult fit = cox_fit(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, optimization_alg, warm_start_fisher_info);
    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = fit.beta[q];
    if (estimate_only) {
        return List::create(_["coefficients"] = coef_r, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
    }
    Eigen::MatrixXd vcov_mat = (cluster.isNotNull()) ? compute_robust_vcov(strata_data, fit.beta, fit.vcov, std::vector<int>(IntegerVector(cluster).begin(), IntegerVector(cluster).end())) : fit.vcov;
    return List::create(_["coefficients"] = coef_r, _["vcov"] = vcov_mat, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
}

// [[Rcpp::export]]
List fast_stratified_coxph_regression_cpp(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& dead,
    const Rcpp::IntegerVector& strata,
    Nullable<NumericVector> warm_start_beta = R_NilValue,
    bool smart_cold_start = true,
    bool estimate_only = false,
    int maxit = 20,
    double tol = 1e-9,
    Nullable<IntegerVector> fixed_idx = R_NilValue,
    Nullable<NumericVector> fixed_values = R_NilValue,
    std::string optimization_alg = "newton_raphson",
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue)
{
    const int n = y.size();
    const int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    // Efficiently group by strata
    std::map<int, std::vector<int>> strata_map;
    for (int i = 0; i < n; ++i) strata_map[strata[i]].push_back(i);

    std::vector<CoxData> strata_data;
    strata_data.reserve(strata_map.size());
    for (auto const& [sid, idx] : strata_map) {
        const int ns = (int)idx.size();
        Eigen::VectorXd y_s(ns), dead_s(ns);
        Eigen::MatrixXd X_s(ns, p);
        for (int ii = 0; ii < ns; ++ii) {
            int id = idx[ii];
            y_s[ii] = y[id]; dead_s[ii] = dead[id]; X_s.row(ii) = X.row(id);
        }
        strata_data.emplace_back(y_s, dead_s, X_s);
    }

    CoxFitResult fit = cox_fit(strata_data, warm_start_beta, smart_cold_start, fixed_spec, estimate_only, maxit, tol, optimization_alg, warm_start_fisher_info);
    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = fit.beta[q];
    if (estimate_only) {
        return List::create(_["coefficients"] = coef_r, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
    }
    return List::create(_["coefficients"] = coef_r, _["vcov"] = fit.vcov, _["converged"] = fit.converged, _["neg_ll"] = fit.neg_ll, _["iterations"] = fit.iterations, _["fisher_information"] = fit.hess_mat);
}

// [[Rcpp::export]]
Eigen::VectorXd get_coxph_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Eigen::VectorXd& beta) {
    std::vector<CoxData> strata_data; strata_data.emplace_back(y, dead, X);
    StratifiedCoxObjective obj(strata_data, X.cols());
    Eigen::VectorXd grad(beta.size()); obj(beta, grad); return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_coxph_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Eigen::VectorXd& beta) {
    std::vector<CoxData> strata_data; strata_data.emplace_back(y, dead, X);
    StratifiedCoxObjective obj(strata_data, X.cols()); return -obj.hessian(beta);
}

// [[Rcpp::export]]
Eigen::VectorXd get_stratified_coxph_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Rcpp::IntegerVector& strata, const Eigen::VectorXd& beta) {
    const int n = y.size(); const int p = X.cols();
    std::map<int, std::vector<int>> strata_map; for (int i = 0; i < n; ++i) strata_map[strata[i]].push_back(i);
    std::vector<CoxData> strata_data; strata_data.reserve(strata_map.size());
    for (auto const& [sid, idx] : strata_map) {
        int ns = (int)idx.size(); Eigen::VectorXd ys(ns), ds(ns); Eigen::MatrixXd Xs(ns, p);
        for (int ii = 0; ii < ns; ++ii) { int id = idx[ii]; ys[ii] = y[id]; ds[ii] = dead[id]; Xs.row(ii) = X.row(id); }
        strata_data.emplace_back(ys, ds, Xs);
    }
    StratifiedCoxObjective obj(strata_data, p); Eigen::VectorXd grad(beta.size()); obj(beta, grad); return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_stratified_coxph_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& dead, const Rcpp::IntegerVector& strata, const Eigen::VectorXd& beta) {
    const int n = y.size(); const int p = X.cols();
    std::map<int, std::vector<int>> strata_map; for (int i = 0; i < n; ++i) strata_map[strata[i]].push_back(i);
    std::vector<CoxData> strata_data; strata_data.reserve(strata_map.size());
    for (auto const& [sid, idx] : strata_map) {
        int ns = (int)idx.size(); Eigen::VectorXd ys(ns), ds(ns); Eigen::MatrixXd Xs(ns, p);
        for (int ii = 0; ii < ns; ++ii) { int id = idx[ii]; ys[ii] = y[id]; ds[ii] = dead[id]; Xs.row(id) = X.row(id); } // Fixed typo here (row(ii))
        strata_data.emplace_back(ys, ds, Xs);
    }
    StratifiedCoxObjective obj(strata_data, p); return -obj.hessian(beta);
}
