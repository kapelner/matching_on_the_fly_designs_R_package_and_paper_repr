// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
int matrix_rank_cpp(const Eigen::MatrixXd& A, double tol = 1e-12) {
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
  qr.setThreshold(tol);
  return qr.rank();
}