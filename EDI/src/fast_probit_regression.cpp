#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

namespace {

// 1/sqrt(2) and 1/sqrt(2*pi) as compile-time constants
static const double kSqrt1_2    = 0.7071067811865476;
static const double k1_Sqrt2Pi  = 0.3989422804014327;

// Use std::erfc directly — avoids R dispatch overhead vs R::pnorm5 in the hot loop.
// Clamped to [-8, 8] so mu stays in (6e-16, 1-6e-16) and phi^2/(mu*(1-mu)) is finite.
inline double pnorm_fast(double x) {
    if (x >= 8.0) return 1.0 - 6e-16;
    if (x <= -8.0) return 6e-16;
    return 0.5 * std::erfc(-x * kSqrt1_2);
}

// phi(x) = exp(-x^2/2) / sqrt(2*pi)
inline double dnorm_fast(double x) {
    return k1_Sqrt2Pi * std::exp(-0.5 * x * x);
}

// Log-scale pnorm: falls back to R for |x| > 6 where direct log(erfc) loses precision.
inline double log_pnorm_lower(double x) {
    if (x > -6.0 && x < 6.0) return std::log(0.5 * std::erfc(-x * kSqrt1_2));
    return R::pnorm5(x, 0.0, 1.0, 1, 1);
}
inline double log_pnorm_upper(double x) {
    if (x > -6.0 && x < 6.0) return std::log(0.5 * std::erfc(x * kSqrt1_2));
    return R::pnorm5(x, 0.0, 1.0, 0, 1);
}

// Probit smart cold start: OLS on probit-transformed mustart = (y + 0.5) / 2
inline Eigen::VectorXd probit_smart_cold_start(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    const int n = X.rows();
    const int p = X.cols();
    Eigen::VectorXd y_adj(n);
    for (int i = 0; i < n; ++i) {
        const double mustart = (y[i] + 0.5) / 2.0;
        y_adj[i] = R::qnorm5(mustart, 0.0, 1.0, 1, 0);
    }
    Eigen::VectorXd beta;
    if (!edi_opt::robust_ols_solve(X, y_adj, beta)) return Eigen::VectorXd::Zero(p);
    return beta;
}

// Probit Fisher information at beta: X^T * diag(phi^2 / (mu*(1-mu))) * X
inline Eigen::MatrixXd probit_smart_hessian(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) {
    const Eigen::VectorXd eta = X * beta;
    const int n = X.rows();
    Eigen::VectorXd w(n);
    for (int i = 0; i < n; ++i) {
        const double e = std::max(-8.0, std::min(8.0, eta[i]));
        const double mu = pnorm_fast(e);
        const double phi = dnorm_fast(e);
        const double vm = std::max(1e-15, mu * (1.0 - mu));
        w[i] = phi * phi / vm;
    }
    return X.transpose() * w.asDiagonal() * X;
}

} // namespace

// Internal probit IRLS core.  Mirrors fast_logistic_regression_internal with probit link.
ModelResult fast_probit_regression_internal(
        const Eigen::MatrixXd& X_eigen,
        const Eigen::VectorXd& y_eigen,
        const Eigen::VectorXd& weights_eigen = Eigen::VectorXd(),
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        int maxit = 100,
        double tol = 1e-8,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
        bool estimate_only = false) {

    const int n = X_eigen.rows();
    const int p = X_eigen.cols();
    const bool use_weights = (weights_eigen.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = static_cast<int>(fixed_spec.free_idx.size());

    Eigen::VectorXd beta_start = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        beta_start = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
        if (static_cast<int>(beta_start.size()) != p)
            Rcpp::stop("warm_start_beta must have length equal to ncol(X)");
    } else if (smart_cold_start) {
        beta_start = probit_smart_cold_start(X_eigen, y_eigen);
    }
    beta_start = apply_fixed_values(beta_start, fixed_spec);

    std::vector<double> eta_fixed_v(n, 0.0);
    for (int k = 0; k < static_cast<int>(fixed_spec.fixed_idx.size()); ++k) {
        const int idx = fixed_spec.fixed_idx[k];
        const double val = fixed_spec.fixed_values[k];
        for (int i = 0; i < n; ++i) eta_fixed_v[i] += X_eigen(i, idx) * val;
    }

    std::vector<double> X_free_raw(n * p_free);
    for (int j = 0; j < p_free; ++j) {
        const int col_idx = fixed_spec.free_idx[j];
        for (int i = 0; i < n; ++i) X_free_raw[i * p_free + j] = X_eigen(i, col_idx);
    }

    std::vector<double> beta_free(p_free);
    for (int j = 0; j < p_free; ++j) beta_free[j] = beta_start[fixed_spec.free_idx[j]];

    std::vector<double> mu_v(n);
    std::vector<double> w_v(n);

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        X_free_map(X_free_raw.data(), n, p_free);
    Eigen::Map<Eigen::VectorXd> mu_map(mu_v.data(), n);
    Eigen::Map<Eigen::VectorXd> w_diag_map(w_v.data(), n);
    Eigen::Map<Eigen::VectorXd> beta_free_map(beta_free.data(), p_free);
    Eigen::Map<const Eigen::VectorXd> eta_fixed_map(eta_fixed_v.data(), n);

    Eigen::MatrixXd final_XtWX(p_free, p_free);
    Eigen::VectorXd score_free(p_free);
    Eigen::VectorXd adj_residual(n);

    bool converged = false;
    int iterations = 0;

    for (int iter = 0; iter < maxit; ++iter) {
        const Eigen::VectorXd eta = X_free_map * beta_free_map + eta_fixed_map;

        // mu = Phi(eta), w = phi^2/(mu*(1-mu)), adj_residual = phi/(mu*(1-mu)) * (y - mu)
        for (int i = 0; i < n; ++i) {
            const double e = eta[i];
            const double mu_i = pnorm_fast(e);
            const double phi_i = dnorm_fast(e);
            const double vm = std::max(1e-15, mu_i * (1.0 - mu_i));
            const double phi_over_vm = phi_i / vm;
            mu_v[i] = mu_i;
            w_v[i] = phi_i * phi_over_vm;
            adj_residual[i] = phi_over_vm * (y_eigen[i] - mu_i);
        }
        if (use_weights) {
            w_diag_map.array() *= weights_eigen.array();
            adj_residual.array() *= weights_eigen.array();
        }

        // Score: X^T * adj_residual
        score_free.noalias() = X_free_map.transpose() * adj_residual;
        if (score_free.norm() < tol) { converged = true; break; }
        ++iterations;

        if (iter == 0 && warm_start_fisher_info.isNotNull()) {
            Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(Rcpp::NumericMatrix(warm_start_fisher_info));
            if (info_full.rows() != p || info_full.cols() != p)
                Rcpp::stop("warm_start_fisher_info must be a p x p matrix");
            final_XtWX = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
        } else if (iter == 0 && smart_cold_start && warm_start_beta.isNull()) {
            final_XtWX = subset_matrix(
                probit_smart_hessian(X_eigen, beta_start),
                fixed_spec.free_idx, fixed_spec.free_idx
            );
        } else {
            final_XtWX.noalias() = X_free_map.transpose() * w_diag_map.asDiagonal() * X_free_map;
        }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(final_XtWX);
        if (ldlt.info() != Eigen::Success) break;
        const Eigen::VectorXd delta = ldlt.solve(score_free);
        if (!delta.allFinite()) break;

        beta_free_map += delta;
        if (delta.norm() < tol) { converged = true; break; }
    }

    ModelResult res;
    res.b = beta_start;
    for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = beta_free_map[j];
    if (!estimate_only) {
        res.mu = mu_map;
        res.XtWX = expand_free_covariance(p, fixed_spec, final_XtWX, false);
    }
    res.iterations = iterations;
    res.converged = converged;
    return res;
}

//' @title Compute Probit Regression Score
//' @description Score vector (gradient of log-likelihood) for a probit regression model.
//' @param X Numeric matrix of predictors (including intercept column).
//' @param y Binary numeric vector of responses.
//' @param beta Numeric vector of coefficients.
//' @return Numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_probit_regression_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& beta) {
    const Eigen::VectorXd eta = X * beta;
    const int n = X.rows();
    Eigen::VectorXd adj_res(n);
    for (int i = 0; i < n; ++i) {
        const double mu_i = pnorm_fast(eta[i]);
        const double phi_i = dnorm_fast(eta[i]);
        const double vm = std::max(1e-15, mu_i * (1.0 - mu_i));
        adj_res[i] = phi_i / vm * (y[i] - mu_i);
    }
    return X.transpose() * adj_res;
}

//' @title Compute Probit Regression Hessian
//' @description Hessian matrix (second derivatives of log-likelihood) for a probit regression model.
//' @param X Numeric matrix of predictors (including intercept column).
//' @param beta Numeric vector of coefficients.
//' @return Numeric matrix representing the Hessian (negative Fisher information).
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_probit_regression_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) {
    const Eigen::VectorXd eta = X * beta;
    const int n = X.rows();
    Eigen::VectorXd w(n);
    for (int i = 0; i < n; ++i) {
        const double mu_i = pnorm_fast(eta[i]);
        const double phi_i = dnorm_fast(eta[i]);
        const double vm = std::max(1e-15, mu_i * (1.0 - mu_i));
        w[i] = phi_i * phi_i / vm;
    }
    return -weighted_crossprod(X, w);
}

//' @title Fast Probit Regression (C++)
//' @description High-performance probit regression fitting using IRLS.
//' @param X Numeric matrix of predictors (including intercept column).
//' @param y Binary numeric vector of responses.
//' @param warm_start_beta Optional starting values for coefficients.
//' @param smart_cold_start Logical. If TRUE, use a probit-transformed OLS guess on cold start.
//' @param maxit Maximum number of IRLS iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional 1-based indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Accepted for API consistency; IRLS is always used.
//' @param warm_start_weights Accepted for API consistency; ignored for probit.
//' @param warm_start_fisher_info Optional initial Fisher information matrix for the first iteration.
//' @return A list with coefficients \code{b}, IRLS weights \code{w}, and \code{fisher_information}.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_probit_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y,
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        int maxit = 100, double tol = 1e-8,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
        bool estimate_only = false) {
    ModelResult res = fast_probit_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta,
        smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg,
        warm_start_weights, warm_start_fisher_info, estimate_only);
    if (estimate_only) {
        return List::create(
            Named("b") = res.b,
            Named("converged") = res.converged,
            Named("iterations") = res.iterations
        );
    }
    // Compute IRLS weights phi^2/(mu*(1-mu)) at the solution
    const int n = X.rows();
    const Eigen::VectorXd eta = X * res.b;
    Eigen::VectorXd weights_vec(n);
    for (int i = 0; i < n; ++i) {
        const double phi_i = dnorm_fast(eta[i]);
        const double vm = std::max(1e-15, res.mu[i] * (1.0 - res.mu[i]));
        weights_vec[i] = phi_i * phi_i / vm;
    }
    return List::create(
        Named("b") = res.b,
        Named("w") = weights_vec,
        Named("iterations") = res.iterations,
        Named("fisher_information") = res.XtWX
    );
}

//' @title Fast Weighted Probit Regression (C++)
//' @description High-performance weighted probit regression fitting using IRLS.
//' @param X Numeric matrix of predictors (including intercept column).
//' @param y Binary numeric vector of responses.
//' @param weights Numeric vector of observation weights.
//' @param warm_start_beta Optional starting values for coefficients.
//' @param smart_cold_start Logical. If TRUE, use a probit-transformed OLS guess on cold start.
//' @param maxit Maximum number of IRLS iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional 1-based indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Accepted for API consistency; IRLS is always used.
//' @param warm_start_weights Accepted for API consistency; ignored for probit.
//' @param warm_start_fisher_info Optional initial Fisher information matrix.
//' @return A list with coefficients, fitted values \code{mu}, and \code{fisher_information}.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_probit_regression_weighted_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y,
        const Eigen::VectorXd& weights,
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        int maxit = 100, double tol = 1e-8,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    ModelResult res = fast_probit_regression_internal(X, y, weights, warm_start_beta,
        smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg,
        warm_start_weights, warm_start_fisher_info);
    return List::create(
        Named("b") = res.b,
        Named("mu") = res.mu,
        Named("XtWX") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}

//' @title Fast Probit Regression with Variance (C++)
//' @description Probit regression with full variance-covariance matrix and likelihood statistics.
//' @param X Numeric matrix of predictors (including intercept column).
//' @param y Binary numeric vector of responses.
//' @param j 1-based index of the parameter for which to return the specific variance.
//' @param warm_start_beta Optional starting values for coefficients.
//' @param smart_cold_start Logical. If TRUE, use a probit-transformed OLS guess on cold start.
//' @param fixed_idx Optional 1-based indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Accepted for API consistency; IRLS is always used.
//' @param warm_start_weights Accepted for API consistency; ignored for probit.
//' @param warm_start_fisher_info Optional initial Fisher information matrix.
//' @return A list with coefficients, variance components, score, and likelihood statistics.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_probit_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int j = 2,
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
        bool smart_cold_start = true,
        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
        std::string optimization_alg = "irls",
        Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    ModelResult res = fast_probit_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta,
        smart_cold_start, 100, 1e-8, fixed_idx, fixed_values, optimization_alg,
        warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);

    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int overall_j) -> int {
        for (int jj = 0; jj < static_cast<int>(fixed_spec.free_idx.size()); ++jj)
            if (fixed_spec.free_idx[jj] == overall_j) return jj + 1; // 1-based
        return -1;
    };

    int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
    res.ssq_b_j = (free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;

    int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1; // column index 1 (0-based) = treatment
    res.ssq_b_2 = (free_2 > 0) ? compute_diagonal_inverse_entry(info_free, free_2) : NA_REAL;

    Eigen::VectorXd score = get_probit_regression_score_cpp(X, y, res.b);
    Eigen::MatrixXd information = res.XtWX;

    // Probit neg log-likelihood: -sum(y*log(Phi(eta)) + (1-y)*log(Phi(-eta)))
    const Eigen::VectorXd eta = X * res.b;
    double neg_loglik = 0.0;
    for (int i = 0; i < eta.size(); ++i) {
        const double log_p  = log_pnorm_lower(eta[i]);  // log(Phi(eta))
        const double log_q  = log_pnorm_upper(eta[i]);  // log(Phi(-eta))
        neg_loglik -= y[i] * log_p + (1.0 - y[i]) * log_q;
    }

    return List::create(
        Named("b") = res.b,
        Named("params") = res.b,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("score") = score,
        Named("observed_information") = information,
        Named("fisher_information") = information,
        Named("information") = information,
        Named("information_type") = "fisher",
        Named("hessian") = -information,
        Named("neg_loglik") = neg_loglik,
        Named("neg_ll") = neg_loglik,
        Named("loglik") = R_finite(neg_loglik) ? -neg_loglik : NA_REAL,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}
