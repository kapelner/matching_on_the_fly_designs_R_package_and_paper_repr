#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

List compute_zhang_match_data_cpp(const IntegerVector& w,
                                  const IntegerVector& m_vec,
                                  const NumericVector& y,
                                  const NumericMatrix& X);

// [[Rcpp::export]]
List match_diffs_cpp(const Eigen::VectorXi& w,
				 		const Eigen::VectorXi& m_vec,
				 		const Eigen::VectorXd& y,
				 		const Eigen::MatrixXd& X,
				 		int m) {
	List match_data = compute_zhang_match_data_cpp(
		wrap(w),
		wrap(m_vec),
		wrap(y),
		wrap(X)
	);

	return List::create(
	_["yTs_matched"] = 		match_data["yTs_matched"],
	_["yCs_matched"] = 		match_data["yCs_matched"],
	_["X_matched_diffs"] = 	match_data["X_matched_diffs"]
	);
}
