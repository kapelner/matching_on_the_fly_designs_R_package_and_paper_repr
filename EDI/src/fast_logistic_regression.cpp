#include "_helper_functions.h"
#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

inline double plogis_manual(double x) {
    if (x > 20.0) return 1.0;
    if (x < -20.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}

inline Eigen::ArrayXd plogis_array_manual(const Eigen::ArrayXd& x) {
    const Eigen::ArrayXd x_clamped = x.max(-20.0).min(20.0);
    return 1.0 / (1.0 + (-x_clamped).exp());
}

// Manual LLT (Cholesky) decomposition and solver to avoid Eigen templates
// A is p*p symmetric, positive definite. b is p*1.
bool solve_llt_raw(double* x, const double* A, const double* b, int p) {
    std::vector<double> L(p * p, 0.0);
    for (int i = 0; i < p; i++) {
        for (int j = 0; j <= i; j++) {
            double s = 0;
            for (int k = 0; k < j; k++) s += L[i * p + k] * L[j * p + k];
            if (i == j) {
                double val = A[i * p + i] - s;
                if (val <= 1e-15) return false;
                L[i * p + j] = std::sqrt(val);
            } else {
                L[i * p + j] = (A[i * p + j] - s) / L[j * p + j];
            }
        }
    }
    std::vector<double> y_tmp(p);
    for (int i = 0; i < p; i++) {
        double s = 0;
        for (int k = 0; k < i; k++) s += L[i * p + k] * y_tmp[k];
        y_tmp[i] = (b[i] - s) / L[i * p + i];
    }
    for (int i = p - 1; i >= 0; i--) {
        double s = 0;
        for (int k = i + 1; k < p; k++) s += L[k * p + i] * x[k];
        x[i] = (y_tmp[i] - s) / L[i * p + i];
    }
    return true;
}

} // namespace

// Internal pure C++ logic
ModelResult fast_logistic_regression_internal(const Eigen::MatrixXd& X_eigen, 
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
                                              Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    int n = X_eigen.rows();
    int p = X_eigen.cols();
    bool use_weights = (weights_eigen.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    int p_free = fixed_spec.free_idx.size();
    Eigen::VectorXd beta_start = Eigen::VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        beta_start = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
        if (beta_start.size() != p) Rcpp::stop("warm_start_beta must have length equal to ncol(X)");
    } else if (smart_cold_start) {
        beta_start = edi_opt::logistic_smart_cold_start(X_eigen, y_eigen);
    }
    beta_start = apply_fixed_values(beta_start, fixed_spec);
    std::vector<double> beta_full(p, 0.0);
    for(int j = 0; j < p; ++j) beta_full[j] = beta_start[j];
    
    std::vector<double> beta_free(p_free, 0.0);
    for(int j=0; j<p_free; j++) beta_free[j] = beta_full[fixed_spec.free_idx[j]];

    std::vector<double> eta_fixed(n, 0.0);
    for(size_t k=0; k<fixed_spec.fixed_idx.size(); k++) {
        int idx = fixed_spec.fixed_idx[k];
        double val = fixed_spec.fixed_values[k];
        for(int i=0; i<n; i++) eta_fixed[i] += X_eigen(i, idx) * val;
    }

    std::vector<double> X_free_raw(n * p_free);
    for (int j = 0; j < p_free; ++j) {
        int col_idx = fixed_spec.free_idx[j];
        for(int i=0; i<n; i++) X_free_raw[i * p_free + j] = X_eigen(i, col_idx);
    }

    std::vector<double> mu(n);
    std::vector<double> w_diag(n);
    bool converged = false;

    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> X_free_map(X_free_raw.data(), n, p_free);
    Eigen::Map<Eigen::VectorXd> mu_map(mu.data(), n);
    Eigen::Map<Eigen::VectorXd> w_diag_map(w_diag.data(), n);
    Eigen::Map<Eigen::VectorXd> beta_free_map(beta_free.data(), p_free);
    Eigen::Map<const Eigen::VectorXd> eta_fixed_map(eta_fixed.data(), n);

    Eigen::MatrixXd final_XtWX_eigen(p_free, p_free);
    Eigen::VectorXd score_free_eigen(p_free);

    int iterations = 0;
    for (int iter = 0; iter < maxit; iter++) {
        iterations = iter + 1;
        
        // mu = plogis(X*beta + eta_fixed)
        mu_map.noalias() = X_free_map * beta_free_map + eta_fixed_map;
        mu_map.array() = plogis_array_safe(mu_map.array());

        // w_diag = weights * mu * (1-mu)
        if (iter == 0 && warm_start_weights.isNotNull()) {
            Eigen::VectorXd ww = as<Eigen::VectorXd>(warm_start_weights);
            if (ww.size() != n) Rcpp::stop("warm_start_weights must have length equal to nrow(X)");
            w_diag_map = ww;
        } else {
            w_diag_map.array() = mu_map.array() * (1.0 - mu_map.array());
            if (use_weights) w_diag_map.array() *= weights_eigen.array();
            w_diag_map.array() = w_diag_map.array().max(1e-10);
        }

        // score = X^T * (weights * (y - mu))
        Eigen::VectorXd diff = y_eigen - mu_map;
        if (use_weights) diff.array() *= weights_eigen.array();
        score_free_eigen.noalias() = X_free_map.transpose() * diff;

        // XtWX = X^T * diag(w_diag) * X
        if (iter == 0 && warm_start_fisher_info.isNotNull()) {
            Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(warm_start_fisher_info);
            if (info_full.rows() != p || info_full.cols() != p) Rcpp::stop("warm_start_fisher_info must be a p x p matrix");
            final_XtWX_eigen = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
        } else if (iter == 0 && smart_cold_start && warm_start_beta.isNull()) {
            Eigen::MatrixXd H_full = edi_opt::logistic_smart_hessian(X_eigen, beta_start);
            final_XtWX_eigen = subset_matrix(H_full, fixed_spec.free_idx, fixed_spec.free_idx);
        } else {
            final_XtWX_eigen.noalias() = X_free_map.transpose() * w_diag_map.asDiagonal() * X_free_map;
        }

        Eigen::LDLT<Eigen::MatrixXd> ldlt(final_XtWX_eigen);
        if (ldlt.info() != Eigen::Success) break;
        Eigen::VectorXd delta = ldlt.solve(score_free_eigen);
        if (!delta.allFinite()) break;

        beta_free_map += delta;
        if (delta.norm() < tol) { converged = true; break; }
    }

    ModelResult res;
    res.b = beta_start; // Copy fixed values
    for(int j=0; j<p_free; j++) res.b[fixed_spec.free_idx[j]] = beta_free_map[j];
    res.mu = mu_map;
    res.XtWX = expand_free_covariance(p, fixed_spec, final_XtWX_eigen, false);
    res.iterations = iterations;
    res.converged = converged;
    return res;
}

//' @title Compute Logistic Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a logistic regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& beta) {
    Eigen::VectorXd mu = plogis_array_manual((X * beta).array()).matrix();
    return X.transpose() * (y - mu);
}

//' @title Compute Logistic Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a logistic regression model.
//' @param X A numeric matrix of predictors.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) {
    Eigen::VectorXd mu = plogis_array_manual((X * beta).array()).matrix();
    Eigen::VectorXd w = mu.array() * (1.0 - mu.array());
    return -weighted_crossprod(X, w);
}

//' @title Compute Weighted Logistic Regression Score
//' @description Calculates the score vector for a weighted logistic regression model.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A numeric vector of weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the weighted score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_logistic_regression_weighted_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& weights, const Eigen::VectorXd& beta) {
    Eigen::VectorXd mu = plogis_array_manual((X * beta).array()).matrix();
    return X.transpose() * weights.cwiseProduct(y - mu);
}

//' @title Compute Weighted Logistic Regression Hessian
//' @description Calculates the Hessian matrix for a weighted logistic regression model.
//' @param X A numeric matrix of predictors.
//' @param weights A numeric vector of weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the weighted Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_regression_weighted_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& weights, const Eigen::VectorXd& beta) {
    Eigen::VectorXd mu = plogis_array_manual((X * beta).array()).matrix();
    Eigen::VectorXd w = weights.array() * mu.array() * (1.0 - mu.array());
    return -weighted_crossprod(X, w);
}

//' @title Fast Logistic Regression (C++)
//' @description High-performance logistic regression fitting using IRLS.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @return A list containing coefficients and weights.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rbinom(10, 1, 0.5)
//' fast_logistic_regression_cpp(X, y)
// [[Rcpp::export]]
List fast_logistic_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y,
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                  bool smart_cold_start = true,
                                  int maxit = 100, double tol = 1e-8,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "irls",
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    Eigen::VectorXd weights_vec(X.rows());
    for(int i=0; i<X.rows(); i++) weights_vec[i] = res.mu[i] * (1.0 - res.mu[i]);
    return List::create(
        Named("b") = res.b,
        Named("w") = weights_vec,
        Named("iterations") = res.iterations,
        Named("fisher_information") = res.XtWX
    );
}

//' @title Fast Weighted Logistic Regression (C++)
//' @description High-performance weighted logistic regression fitting using IRLS.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param weights A numeric vector of weights.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm to use.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, fitted values, and information matrix.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_logistic_regression_weighted_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& weights,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                           bool smart_cold_start = true,
                                           int maxit = 100, double tol = 1e-8,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls",
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    ModelResult res = fast_logistic_regression_internal(X, y, weights, warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    return List::create(
        Named("b") = res.b,
        Named("mu") = res.mu,
        Named("XtWX") = res.XtWX,
        Named("fisher_information") = res.XtWX,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations
    );
}

//' @title Fast Logistic Regression with Variance (C++)
//' @description Logistic regression with full variance-covariance matrix and score calculation.
//' @param X A numeric matrix of predictors.
//' @param y A binary numeric vector of responses.
//' @param j The 1-based index of the parameter for which to return specific variance.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @return A list containing coefficients, vcov, score, and likelihood statistics.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rbinom(10, 1, 0.5)
//' fast_logistic_regression_with_var_cpp(X, y)
// [[Rcpp::export]]
List fast_logistic_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int j = 2,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                           bool smart_cold_start = true,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls",
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, 100, 1e-8, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    
    int p_free = fixed_spec.free_idx.size();
    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    std::vector<double> info_free_v(p_free * p_free);
    for(int r=0; r<p_free; r++) for(int c=0; c<p_free; c++) info_free_v[r * p_free + c] = info_free(r, c);

    Eigen::MatrixXd cov_free = Eigen::MatrixXd::Zero(p_free, p_free);
    for(int col=0; col<p_free; col++) {
        std::vector<double> b(p_free, 0.0);
        b[col] = 1.0;
        std::vector<double> x_sol(p_free);
        solve_llt_raw(x_sol.data(), info_free_v.data(), b.data(), p_free);
        for(int r=0; r<p_free; r++) cov_free(r, col) = x_sol[r];
    }
    
    Eigen::MatrixXd vcov = expand_free_covariance(X.cols(), fixed_spec, cov_free, true);
    res.ssq_b_j = (j > 0 && j <= X.cols()) ? vcov(j - 1, j - 1) : NA_REAL;
    res.ssq_b_2 = (X.cols() >= 2) ? vcov(1, 1) : NA_REAL;
    Eigen::VectorXd score = get_logistic_regression_score_cpp(X, y, res.b);
    Eigen::MatrixXd information = res.XtWX;
    Eigen::VectorXd eta = X * res.b;
    double neg_loglik = 0.0;
    for (int i = 0; i < eta.size(); ++i) {
        double log_one_plus_exp_eta = (eta[i] > 0.0) ?
            eta[i] + std::log1p(std::exp(-eta[i])) :
            std::log1p(std::exp(eta[i]));
        neg_loglik += log_one_plus_exp_eta - y[i] * eta[i];
    }

    return List::create(
        Named("b") = res.b,
        Named("params") = res.b,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("vcov") = vcov,
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
