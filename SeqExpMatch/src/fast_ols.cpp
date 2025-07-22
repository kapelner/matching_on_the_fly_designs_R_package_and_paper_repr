#include "_helper_functions.h"
using namespace Rcpp;

// [[Rcpp::export]]
List fast_ols_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
  Eigen::MatrixXd XtX = X.transpose() * X;
  Eigen::VectorXd Xty = X.transpose() * y;

  // Solve XtX * beta = Xty
  Eigen::VectorXd beta = XtX.ldlt().solve(Xty);  // LDLT is stable and fast

  return List::create(
    Named("b") = beta,
    Named("XtX") = XtX
  );
}

// [[Rcpp::export]]
List fast_ols_with_sd_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
  int n = X.rows();
  int p = X.cols();
  List mod = fast_ols_cpp(X, y);
  Eigen::VectorXd beta = mod["b"];

  // Residuals and RSS
  Eigen::VectorXd residuals = y - X * beta;
  double rss = residuals.squaredNorm();     // == residuals.dot(residuals)
  double sigma2_hat = rss / (n - p);        // Estimated residual variance

  mod["s_b_2"] = sqrt(sigma2_hat * eigen_compute_single_entry_of_diagonal_matrix_cpp(mod["XtX"], 2));
  return mod;
}
