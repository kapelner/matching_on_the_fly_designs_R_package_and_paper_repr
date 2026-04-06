#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

// Internal pure C++ logic
ModelResult fast_ols_internal(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    ModelResult res;
    res.XtWX = X.transpose() * X; // XtWX is just XtX here
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
    res.b = cod.solve(y);
    if (!res.b.allFinite()) {
        res.b = Eigen::VectorXd::Constant(X.cols(), NA_REAL);
    }
    return res;
}

// [[Rcpp::export]]
List fast_ols_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    ModelResult res = fast_ols_internal(X, y);
    return List::create(
        Named("b") = res.b,
        Named("XtX") = res.XtWX
    );
}

// [[Rcpp::export]]
List fast_ols_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int j = 2) {
    int n = X.rows();
    int p = X.cols();
    ModelResult res = fast_ols_internal(X, y);

    Eigen::VectorXd e = y - X * res.b;
    double sse = e.squaredNorm();
    res.sigma2_hat = sse / (n - p);

    res.ssq_b_j = res.sigma2_hat * compute_diagonal_inverse_entry(res.XtWX, j);
    res.ssq_b_2 = (X.cols() >= 2) ? (res.sigma2_hat * compute_diagonal_inverse_entry(res.XtWX, 2)) : NA_REAL;

    return List::create(
        Named("b") = res.b,
        Named("XtX") = res.XtWX,
        Named("ssq_b_j") = res.ssq_b_j,
        Named("ssq_b_2") = res.ssq_b_2,
        Named("sigma2_hat") = res.sigma2_hat
    );
}
