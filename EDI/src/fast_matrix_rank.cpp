// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
int matrix_rank_cpp(const NumericMatrix& A_r, double tol = 1e-7) {
	Eigen::Map<const Eigen::MatrixXd> A(A_r.begin(), A_r.rows(), A_r.cols());
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
	qr.setThreshold(tol);
	return qr.rank();
}
