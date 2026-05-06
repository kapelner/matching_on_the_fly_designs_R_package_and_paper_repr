#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd compute_proportional_mahal_distances_cpp(
	const Eigen::VectorXd& xt_prev,              // fixed vector
	const Eigen::MatrixXd& X_prev,               // full data matrix
	const Eigen::VectorXi& reservoir_indices,    // indices (1-based from R)
	const Eigen::MatrixXd& S_xs_inv)             // inverse covariance matrix
{
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
