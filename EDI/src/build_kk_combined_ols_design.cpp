#include "_helper_functions.h"
using namespace Rcpp;

// Build the combined OLS design matrix and response vector for the KK
// combined-likelihood estimator (matched pairs + reservoir, both present).
//
// Column layout: [beta_0 | beta_T | beta_xs (p cols)]
//   Pair rows (rows 0..m-1, scaled by 1/sqrt2):
//     [0 | 1/sqrt2 | Xd_k / sqrt2]
//     y_comb[k] = yd[k] / sqrt2
//   Reservoir rows (rows m..m+nR-1):
//     [1 | w_r[i] | X_r[i, ]]
//     y_comb[m+i] = y_r[i]
//
// Scaling by 1/sqrt2 equalises all residual variances to sigma^2.
//
// [[Rcpp::export]]
List build_kk_combined_ols_design_cpp(
	const Eigen::VectorXd& yd,
	const Eigen::MatrixXd& Xd,
	const Eigen::VectorXd& y_r,
	const Eigen::VectorXd& w_r,
	const Eigen::MatrixXd& X_r
) {
	const int m   = (int)yd.size();
	const int nR  = (int)y_r.size();
	const int p   = (int)Xd.cols();
	const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

	Eigen::MatrixXd X_comb(m + nR, p + 2);
	Eigen::VectorXd y_comb(m + nR);

	// ---- Pair rows -------------------------------------------------------
	X_comb.block(0, 0, m, 1).setZero();
	X_comb.block(0, 1, m, 1).setConstant(inv_sqrt2);
	X_comb.block(0, 2, m, p) = Xd * inv_sqrt2;
	y_comb.head(m)            = yd * inv_sqrt2;

	// ---- Reservoir rows --------------------------------------------------
	X_comb.block(m, 0, nR, 1).setOnes();
	X_comb.block(m, 1, nR, 1) = w_r;
	X_comb.block(m, 2, nR, p) = X_r;
	y_comb.tail(nR)            = y_r;

	return List::create(
		Named("X_comb") = X_comb,
		Named("y_comb") = y_comb
	);
}
