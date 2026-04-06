#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;

namespace {

int atkinson_assign_weight_internal(
	const Eigen::VectorXd& w_prev,
	const Eigen::MatrixXd& X_prev,
	const Eigen::VectorXd& xt_prev,
	int t
) {
	int rows = w_prev.size();
	int p = X_prev.cols();
	int cols = p + 2;
	if (rows == 0 || cols < 2) {
	    return (R::unif_rand() < 0.5) ? 1 : 0;
	}

	Eigen::MatrixXd XprevWT(rows, cols);
    XprevWT.col(0) = w_prev;
    XprevWT.col(1).setOnes();
    XprevWT.rightCols(p) = X_prev;

	Eigen::MatrixXd XwtXw = XprevWT.transpose() * XprevWT;
	Eigen::FullPivLU<Eigen::MatrixXd> lu(XwtXw);
	if (!lu.isInvertible()) {
	    return (R::unif_rand() < 0.5) ? 1 : 0;
	}

	Eigen::MatrixXd M = static_cast<double>(t - 1) * lu.inverse();
	Eigen::VectorXd row_segment = M.row(0).segment(1, p + 1);

	Eigen::VectorXd xt(p + 1);
	xt(0) = 1.0;
    xt.tail(p) = xt_prev;

	double A = row_segment.dot(xt);
	if (A == 0 || !std::isfinite(A)) {
	    return (R::unif_rand() < 0.5) ? 1 : 0;
	}

	double val = M(0, 0) / A + 1.0;
	double s_over_A_plus_one_sq = val * val;
	double prob = s_over_A_plus_one_sq / (s_over_A_plus_one_sq + 1.0);
	prob = std::max(0.0, std::min(1.0, prob));
	return (R::unif_rand() < prob) ? 1 : 0;
}

} // namespace

// [[Rcpp::export]]
double atkinson_assign_weight_cpp(
	const NumericVector& w_prev,
	const NumericMatrix& X_prev,
	const NumericVector& xt_prev,
	int rank_prev,
	int t
) {
    Eigen::VectorXd w_e = as<Eigen::VectorXd>(w_prev);
    Eigen::MatrixXd X_e = as<Eigen::MatrixXd>(X_prev);
    Eigen::VectorXd xt_e = as<Eigen::VectorXd>(xt_prev);
    
    return static_cast<double>(atkinson_assign_weight_internal(w_e, X_e, xt_e, t));
}
