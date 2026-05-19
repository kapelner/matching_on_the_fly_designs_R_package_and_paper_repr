#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

// Internal pure C++ logic
ModelResult fast_ols_internal(const Eigen::MatrixXd& X,
                              const Eigen::VectorXd& y,
                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
    const int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = fixed_spec.free_idx.size();
    Eigen::MatrixXd X_free(X.rows(), p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
    Eigen::VectorXd y_adj = y;
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        y_adj.noalias() -= X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }

    ModelResult res;
    Eigen::MatrixXd XtX_free = X_free.transpose() * X_free;
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X_free);
    Eigen::VectorXd beta_free = cod.solve(y_adj);
    res.b = Eigen::VectorXd::Zero(p);
    for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = beta_free[j];
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];
    if (!res.b.allFinite()) {
        res.b = Eigen::VectorXd::Constant(p, NA_REAL);
    }
    res.XtWX = expand_free_covariance(p, fixed_spec, XtX_free, false);
    return res;
}

//' @title Fast Ordinary Least Squares (C++)
//' @description High-performance OLS fitting using Eigen's Complete Orthogonal Decomposition.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients and the XtX matrix.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rnorm(10)
//' fast_ols_cpp(X, y)
// [[Rcpp::export]]
List fast_ols_cpp(const Eigen::MatrixXd& X,
                  const Eigen::VectorXd& y,
                  Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
    ModelResult res = fast_ols_internal(X, y, fixed_idx, fixed_values);
    return List::create(
        Named("b") = res.b,
        Named("XtX") = res.XtWX
    );
}

//' @title Fast OLS with Variance (C++)
//' @description OLS fitting with variance-covariance matrix and error variance estimation.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @return A list containing coefficients, vcov, ssq_b_j, and sigma2_hat.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ols_with_var_cpp(const Eigen::MatrixXd& X,
                           const Eigen::VectorXd& y,
                           int j = 2,
                           Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                           Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue) {
    int n = X.rows();
    int p = X.cols();
    ModelResult res = fast_ols_internal(X, y, fixed_idx, fixed_values);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);

    Eigen::VectorXd e = y - X * res.b;
    double sse = e.squaredNorm();
    res.sigma2_hat = sse / (n - fixed_spec.free_idx.size());

    Eigen::MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int k) -> int {
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == k) return jj + 1;
        return -1;
    };
    int free_j = (j > 0 && j <= p) ? free_idx_of(j - 1) : -1;
    double raw_j = (free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;
    res.ssq_b_j = R_finite(raw_j) ? res.sigma2_hat * raw_j : NA_REAL;
    int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;
    double raw_2 = (free_2 > 0) ? compute_diagonal_inverse_entry(info_free, free_2) : NA_REAL;
    res.ssq_b_2 = R_finite(raw_2) ? res.sigma2_hat * raw_2 : NA_REAL;

    return List::create(
        Named("b") = res.b,
        Named("XtX") = res.XtWX,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("sigma2_hat") = res.sigma2_hat
    );
}
