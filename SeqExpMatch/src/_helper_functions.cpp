#include <RcppEigen.h>
using namespace Rcpp;


typedef Eigen::Map<Eigen::MatrixXd> MapMat;
typedef Eigen::Map<Eigen::VectorXd> MapVec;

// [[Rcpp::export]]
double eigen_compute_single_entry_of_diagonal_matrix_cpp(MapMat M, int j) {
  Eigen::VectorXd b;
  b.resize(M.rows());
  b.setZero();
  b(j - 1) = 1;

  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper> cg;
  cg.compute(M);
  Eigen::VectorXd x = cg.solve(b);

  return x(j - 1);
}

// [[Rcpp::export]]
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(MapMat X, MapVec w) {
  return X.transpose() * w.asDiagonal() * X;
}