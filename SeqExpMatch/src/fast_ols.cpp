#include "_helper_functions.h"
using namespace Rcpp;

// [[Rcpp::export]]
List fast_ols_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int j = 2) {
  int n = X.rows();
  int p = X.cols();
  List mod = fast_ols_cpp(X, y);
  Eigen::VectorXd beta = mod["b"];

  // Residuals and RSS
  Eigen::VectorXd e = y - X * beta;
  double sse = e.squaredNorm();     // == e.dot(residuals)
  double sigma2_hat = sse / (n - p);        // Estimated residual variance

  mod["ssq_b_j"] = sigma2_hat * eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(mod["XtX"], j);
  return mod;
}
