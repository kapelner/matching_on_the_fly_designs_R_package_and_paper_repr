//#include <RcppEigen.h>
#include "_helper_functions.h"
using namespace Rcpp;


//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
//typedef Eigen::Map<Eigen::VectorXd> MapVec;

// [[Rcpp::export]]
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j) {
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
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w) {
  return X.transpose() * w.asDiagonal() * X;
}


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