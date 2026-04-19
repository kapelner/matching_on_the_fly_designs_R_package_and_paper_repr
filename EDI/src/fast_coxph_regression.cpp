#include <RcppEigen.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

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

} // namespace

// [[Rcpp::export]]
List fast_coxph_regression_cpp(const Eigen::VectorXd& y,
                               const Eigen::VectorXd& dead,
                               const Eigen::MatrixXd& X,
                               Nullable<NumericVector> start_beta = R_NilValue,
                               bool estimate_only = false,
                               int maxit = 20,
                               double tol = 1e-9) {
    int p = X.cols();
    std::vector<double> beta(p, 0.0);
    if (start_beta.isNotNull()) {
        NumericVector sb(start_beta);
        for (int q = 0; q < p; ++q) beta[q] = sb[q];
    }

    CoxData data(y, dead, X);

    // Pre-allocate all work buffers once
    std::vector<double> grad(p), hess(p * p);
    std::vector<double> exp_eta(data.n), r_x_exp(p), r_xx_exp(p * p);
    std::vector<double> sum_x_dk(p), e_z(p);

    double old_ll = 1e300;
    int iter = 0;

    for (iter = 0; iter < maxit; ++iter) {
        double ll = compute_cox_ll_grad_hess_fast(
            data, beta, grad, hess, false,
            exp_eta, r_x_exp, r_xx_exp, sum_x_dk, e_z);

        if (std::abs(old_ll - ll) < tol) break;
        old_ll = ll;

        // Newton step: solve hess * delta = grad, then beta -= delta
        // Use Eigen for the p×p solve (cheap and correct)
        Eigen::Map<const Eigen::MatrixXd> H(hess.data(), p, p);
        Eigen::Map<const Eigen::VectorXd> g(grad.data(), p);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
        if (ldlt.info() != Eigen::Success) break;
        Eigen::VectorXd delta = ldlt.solve(g);
        for (int q = 0; q < p; ++q) beta[q] -= delta[q];
    }

    // Map beta to Rcpp
    NumericVector coef_r(p);
    for (int q = 0; q < p; ++q) coef_r[q] = beta[q];

    if (estimate_only) {
        return List::create(
            Named("coefficients") = coef_r,
            Named("converged") = (iter < maxit),
            Named("neg_ll") = old_ll,
            Named("iterations") = iter
        );
    }

    // One final pass to get the Hessian at the converged beta for vcov
    compute_cox_ll_grad_hess_fast(
        data, beta, grad, hess, false,
        exp_eta, r_x_exp, r_xx_exp, sum_x_dk, e_z);

    // vcov = hess^{-1}
    Eigen::Map<const Eigen::MatrixXd> H_final(hess.data(), p, p);
    Eigen::MatrixXd vcov = H_final.inverse();

    return List::create(
        Named("coefficients") = coef_r,
        Named("vcov") = vcov,
        Named("converged") = (iter < maxit),
        Named("neg_ll") = old_ll,
        Named("iterations") = iter
    );
}
