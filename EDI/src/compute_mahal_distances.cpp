#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd compute_proportional_mahal_distances_cpp(
	SEXP xt_prev_sexp,              // fixed vector
	SEXP X_prev_sexp,               // full data matrix
	SEXP reservoir_indices_sexp,    // indices (1-based from R)
	SEXP S_xs_inv_sexp)             // inverse covariance matrix
{
	Rcpp::NumericVector xt_prev_r(xt_prev_sexp);
	Rcpp::NumericMatrix X_prev_r(X_prev_sexp);
	Rcpp::IntegerVector reservoir_indices_r(reservoir_indices_sexp);
	Rcpp::NumericMatrix S_xs_inv_r(S_xs_inv_sexp);
	Eigen::Map<const Eigen::VectorXd> xt_prev(xt_prev_r.begin(), xt_prev_r.size());
	Eigen::Map<const Eigen::MatrixXd> X_prev(X_prev_r.begin(), X_prev_r.nrow(), X_prev_r.ncol());
	Eigen::Map<const Eigen::VectorXi> reservoir_indices(reservoir_indices_r.begin(), reservoir_indices_r.size());
	Eigen::Map<const Eigen::MatrixXd> S_xs_inv(S_xs_inv_r.begin(), S_xs_inv_r.nrow(), S_xs_inv_r.ncol());
	const int n_R = reservoir_indices.size();
	const int p   = xt_prev.size();

	// Build difference matrix D (n_R x p) in one pass
	Eigen::MatrixXd D(n_R, p);
	for (int r = 0; r < n_R; r++)
		D.row(r) = xt_prev.transpose() - X_prev.row(reservoir_indices[r] - 1);

	// All quadratic forms d'*S_inv*d in one batched Eigen call:
	//   (D * S_inv) elementwise-product D, summed per row
	return (D * S_xs_inv).cwiseProduct(D).rowwise().sum();
}
