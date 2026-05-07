#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

List compute_zhang_match_data_cpp(const NumericMatrix& X,
                                  const NumericVector& y,
                                  const IntegerVector& w,
                                  const IntegerVector& m_vec);

// [[Rcpp::export]]
List match_diffs_cpp(const Eigen::MatrixXd& X,
                        const Eigen::VectorXd& y,
                        const Eigen::VectorXi& w,
				 		const Eigen::VectorXi& m_vec,
				 		int m) {
	List match_data = compute_zhang_match_data_cpp(
		wrap(X),
		wrap(y),
		wrap(w),
		wrap(m_vec)
	);

	return List::create(
	_["yTs_matched"] = 		match_data["yTs_matched"],
	_["yCs_matched"] = 		match_data["yCs_matched"],
	_["X_matched_diffs"] = 	match_data["X_matched_diffs"]
	);
}
