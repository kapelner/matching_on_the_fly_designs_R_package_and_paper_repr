#include "_helper_functions.h"
using namespace Rcpp;


//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
//typedef Eigen::Map<Eigen::VectorXd> MapVec;

// [[Rcpp::export]]
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j) {
  Eigen::VectorXd b = Eigen::VectorXd::Unit(M.rows(), j - 1);

  // Use a direct LDLT decomposition instead of the iterative ConjugateGradient
  // solver. CG has a fixed iteration budget and fails for poorly-conditioned but
  // invertible matrices (e.g. scaled factor-dummy designs). LDLT is exact for
  // symmetric positive-definite matrices and has no convergence issues.
  Eigen::LDLT<Eigen::MatrixXd> ldlt(M);

  if (ldlt.info() != Eigen::Success) {
    return NA_REAL;
  }

  Eigen::VectorXd x = ldlt.solve(b);

  if (!x.allFinite()) {
    return NA_REAL;
  }

  return x(j - 1);
}

// [[Rcpp::export]]
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w) {
  return X.transpose() * w.asDiagonal() * X;
}

// [[Rcpp::export]]
double mean_cpp(const Eigen::VectorXd& x) {
  if (x.size() == 0) {
    return NA_REAL;
  }
  return x.mean();
}

// [[Rcpp::export]]
double var_cpp(const Eigen::VectorXd& x) {
  if (x.size() <= 1) {
    return NA_REAL;
  }
  return (x.array() - x.mean()).square().sum() / (x.size() - 1);
}
