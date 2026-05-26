#include "_helper_functions.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
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
    return plogis_array_safe(x);
}

inline double log1pexp_stable(double x) {
    return (x > 0.0) ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

class LogisticLbfgsObjective : public Numer::MFuncGrad {
private:
    const RowMajorMatrixXd& X;
    const Eigen::VectorXd& y;
    const Eigen::VectorXd& weights;
    const Eigen::VectorXd& eta_fixed;
    const bool use_weights;

public:
    LogisticLbfgsObjective(const RowMajorMatrixXd& X_,
                           const Eigen::VectorXd& y_,
                           const Eigen::VectorXd& weights_,
                           const Eigen::VectorXd& eta_fixed_,
                           bool use_weights_) :
        X(X_), y(y_), weights(weights_), eta_fixed(eta_fixed_), use_weights(use_weights_) {}

    double f_grad(Numer::Constvec& beta, Numer::Refvec grad) {
        Eigen::VectorXd eta = X * beta + eta_fixed;
        Eigen::VectorXd resid(eta.size());
        double f = 0.0;

        for (int i = 0; i < eta.size(); ++i) {
            const double wi = use_weights ? weights[i] : 1.0;
            const double mui = plogis_manual(eta[i]);
            f += wi * (log1pexp_stable(eta[i]) - y[i] * eta[i]);
            resid[i] = wi * (mui - y[i]);
        }

        grad.noalias() = X.transpose() * resid;
        return f;
    }
};

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
                                              bool smart_cold_start = false,
                                              int maxit = 100, 
                                              double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls",
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                              bool estimate_only = false) {
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

    if (optimization_alg == "lbfgs") {
        RowMajorMatrixXd X_free_eigen(n, p_free);
        for (int j = 0; j < p_free; ++j) {
            X_free_eigen.col(j) = X_eigen.col(fixed_spec.free_idx[j]);
        }

        Eigen::VectorXd beta_free_eigen(p_free);
        for (int j = 0; j < p_free; ++j) {
            beta_free_eigen[j] = beta_free[j];
        }
        Eigen::VectorXd eta_fixed_eigen = Eigen::Map<const Eigen::VectorXd>(eta_fixed.data(), n);

        bool converged = true;
        int status = 0;
        double fopt = NA_REAL;
        if (p_free > 0) {
            LogisticLbfgsObjective nll(X_free_eigen, y_eigen, weights_eigen, eta_fixed_eigen, use_weights);
            status = Numer::optim_lbfgs(nll, beta_free_eigen, fopt, maxit, tol, tol);
            converged = (status >= 0) && beta_free_eigen.allFinite();
        }

        ModelResult res;
        res.b = beta_start;
        for (int j = 0; j < p_free; ++j) {
            res.b[fixed_spec.free_idx[j]] = beta_free_eigen[j];
        }

        if (!estimate_only) {
            Eigen::VectorXd eta = X_eigen * res.b;
            res.mu = plogis_array_safe(eta.array()).matrix();
            Eigen::VectorXd w_diag = res.mu.array() * (1.0 - res.mu.array());
            if (use_weights) w_diag.array() *= weights_eigen.array();
            w_diag.array() = w_diag.array().max(1e-10);

            Eigen::MatrixXd info_free(p_free, p_free);
            if (p_free > 0) {
                info_free.noalias() = X_free_eigen.transpose() * w_diag.asDiagonal() * X_free_eigen;
            } else {
                info_free.resize(0, 0);
            }
            res.XtWX = expand_free_covariance(p, fixed_spec, info_free, false);
        }
        res.iterations = NA_INTEGER;
        res.converged = converged;
        return res;
    }

    RowMajorMatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) {
        X_free.col(j) = X_eigen.col(fixed_spec.free_idx[j]);
    }

    Eigen::VectorXd mu(n);
    Eigen::VectorXd w_diag(n);
    Eigen::VectorXd beta_free_vec = Eigen::Map<Eigen::VectorXd>(beta_free.data(), p_free);
    Eigen::VectorXd eta_fixed_vec = Eigen::Map<Eigen::VectorXd>(eta_fixed.data(), n);
    bool converged = false;

    Eigen::MatrixXd XtWX = Eigen::MatrixXd::Zero(p_free, p_free);
    int iterations = 0;
    for (int iter = 0; iter < maxit; iter++) {
        iterations++;
        
        // mu = plogis(X*beta + eta_fixed)
        Eigen::VectorXd eta = X_free * beta_free_vec + eta_fixed_vec;
        mu.array() = plogis_array_safe(eta.array());

        // w_diag = weights * mu * (1-mu)
        if (iter == 0 && warm_start_weights.isNotNull()) {
            Eigen::VectorXd ww = as<Eigen::VectorXd>(warm_start_weights);
            if (ww.size() == n) w_diag = ww;
            else w_diag.array() = mu.array() * (1.0 - mu.array());
        } else {
            w_diag.array() = mu.array() * (1.0 - mu.array());
        }
        if (use_weights) w_diag.array() *= weights_eigen.array();
        w_diag.array() = w_diag.array().max(1e-10);

        // score = X^T * (weights * (y - mu))
        Eigen::VectorXd diff = y_eigen - mu;
        if (use_weights) diff.array() *= weights_eigen.array();
        
        Eigen::VectorXd score = X_free.transpose() * diff;
        if (score.norm() < tol) {
            converged = true;
        }

        // XtWX = X^T * diag(w_diag) * X
        if (iter == 0 && warm_start_fisher_info.isNotNull()) {
            Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(warm_start_fisher_info);
            XtWX = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
        } else if (iter == 0 && smart_cold_start && warm_start_beta.isNull()) {
            XtWX = subset_matrix(edi_opt::logistic_smart_hessian(X_eigen, beta_start), fixed_spec.free_idx, fixed_spec.free_idx);
        } else {
            XtWX = weighted_crossprod(X_free, w_diag);
        }

        if (converged) break;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) break;
        Eigen::VectorXd delta = ldlt.solve(score);
        if (!delta.allFinite()) break;

        beta_free_vec += delta;
        if (delta.norm() < tol) { converged = true; break; }
    }

    ModelResult res;
    res.b = beta_start;
    for(int j=0; j<p_free; j++) res.b[fixed_spec.free_idx[j]] = beta_free_vec[j];
    if (!estimate_only) {
        res.mu = mu;
        res.XtWX = expand_free_covariance(p, fixed_spec, XtWX, false);
    }
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
//' @description High-performance logistic regression fitting using L-BFGS or IRLS.
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
List fast_logistic_regression_cpp(SEXP X_sexp, SEXP y_sexp,
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                  bool smart_cold_start = false,
                                  int maxit = 100, double tol = 1e-8,
                                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                  std::string optimization_alg = "irls",
                                  Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                  bool estimate_only = false) {
    Eigen::MatrixXd X = as<Eigen::MatrixXd>(X_sexp);
    Eigen::VectorXd y = as<Eigen::VectorXd>(y_sexp);

    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
    
    if (estimate_only) {
        return List::create(
            Named("b") = res.b,
            Named("converged") = res.converged,
            Named("iterations") = res.iterations
        );
    }
    Eigen::VectorXd weights_vec = res.mu.array() * (1.0 - res.mu.array());
    return List::create(
        Named("b") = res.b,
        Named("w") = weights_vec,
        Named("iterations") = res.iterations,
        Named("fisher_information") = res.XtWX
    );
}

//' @title Fast Weighted Logistic Regression (C++)
//' @description High-performance weighted logistic regression fitting using L-BFGS or IRLS.
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
List fast_logistic_regression_weighted_cpp(SEXP X_sexp, SEXP y_sexp, SEXP weights_sexp,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                           bool smart_cold_start = false,
                                           int maxit = 100, double tol = 1e-8,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls",
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    Eigen::MatrixXd X = as<Eigen::MatrixXd>(X_sexp);
    Eigen::VectorXd y = as<Eigen::VectorXd>(y_sexp);
    Eigen::VectorXd weights = as<Eigen::VectorXd>(weights_sexp);
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
List fast_logistic_regression_with_var_cpp(SEXP X_sexp, SEXP y_sexp, int j = 2,
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                           bool smart_cold_start = false,
                                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                           std::string optimization_alg = "irls",
                                           Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                           Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    Eigen::MatrixXd X = as<Eigen::MatrixXd>(X_sexp);
    Eigen::VectorXd y = as<Eigen::VectorXd>(y_sexp);
    ModelResult res = fast_logistic_regression_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, 100, 1e-8, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    
    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int overall_j) -> int {
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == overall_j) return jj + 1; // 1-based for compute_diagonal_inverse_entry
        return -1;
    };

    int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
    res.ssq_b_j = (free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;

    int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;
    res.ssq_b_2 = (free_2 > 0) ? compute_diagonal_inverse_entry(info_free, free_2) : NA_REAL;
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
