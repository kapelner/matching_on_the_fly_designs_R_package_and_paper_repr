#include "_helper_functions.h"
using namespace Rcpp;

// [[Rcpp::export]]
List fast_ols_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
  // Use normal equations with LDLT for speed (X'X)beta = X'y
  Eigen::MatrixXd XtX = X.transpose() * X;
  Eigen::VectorXd Xty = X.transpose() * y;

  // LDLT is much faster than QR for small to medium problems
  Eigen::VectorXd beta = XtX.ldlt().solve(Xty);

  return List::create(
    Named("b") = beta,
    Named("XtX") = XtX  // Return XtX for reuse in variance computation
  );
}

// [[Rcpp::export]]
List fast_ols_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int j = 2) {
  int n = X.rows();
  int p = X.cols();
  List mod = fast_ols_cpp(X, y);
  Eigen::VectorXd beta = mod["b"];
  Eigen::MatrixXd XtX = mod["XtX"];

  // Residuals and RSS
  Eigen::VectorXd e = y - X * beta;
  double sse = e.squaredNorm();
  double sigma2_hat = sse / (n - p);

  // Compute variance using the XtX we already have
  mod["ssq_b_j"] = sigma2_hat * eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtX, j);
  return mod;
}
