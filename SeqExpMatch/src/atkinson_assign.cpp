#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
double atkinson_assign_weight_cpp(
  const NumericVector& w_prev,
  const NumericMatrix& X_prev,
  const NumericVector& xt_prev,
  int rank_prev,
  int t
) {
  int rows = w_prev.size();
  int cols = rank_prev + 2;
  if (rows == 0 || cols < 2) {
    return R::rbinom(1, 0.5);
  }

  Eigen::MatrixXd XprevWT(rows, cols);
  for (int i = 0; i < rows; ++i) {
    XprevWT(i, 0) = w_prev[i];
    XprevWT(i, 1) = 1.0;
    for (int j = 0; j < rank_prev; ++j) {
      XprevWT(i, j + 2) = X_prev(i, j);
    }
  }

  Eigen::MatrixXd XwtXw = XprevWT.transpose() * XprevWT;
  Eigen::FullPivLU<Eigen::MatrixXd> lu(XwtXw);
  if (!lu.isInvertible()) {
    return R::rbinom(1, 0.5);
  }

  Eigen::MatrixXd M = static_cast<double>(t - 1) * lu.inverse();
  Eigen::VectorXd row_segment(rank_prev + 1);
  for (int j = 0; j < rank_prev + 1; ++j) {
    row_segment(j) = M(0, j + 1);
  }

  Eigen::VectorXd xt(rank_prev + 1);
  xt(0) = 1.0;
  for (int j = 0; j < rank_prev; ++j) {
    xt(j + 1) = xt_prev[j];
  }
  double A = row_segment.dot(xt);
  if (A == 0 || !std::isfinite(A)) {
    return R::rbinom(1, 0.5);
  }

  double val = M(0, 0) / A + 1.0;
  double s_over_A_plus_one_sq = val * val;
  double prob = s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1.0);
  prob = std::max(0.0, std::min(1.0, prob));
  return R::rbinom(1, prob);
}
