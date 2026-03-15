#include "_helper_functions.h"
using namespace Rcpp;


//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
//typedef Eigen::Map<Eigen::VectorXd> MapVec;

// [[Rcpp::export]]
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j) {
	Eigen::VectorXd b = Eigen::VectorXd::Unit(M.rows(), j - 1);

	// Prefer a direct LDLT decomposition for well-behaved symmetric systems.
	Eigen::LDLT<Eigen::MatrixXd> ldlt(M);

	if (ldlt.info() == Eigen::Success) {
		Eigen::VectorXd x = ldlt.solve(b);
		if (x.allFinite()) {
			return x(j - 1);
		}
	}

	// Fall back to a rank-revealing solve for near-singular or indefinite systems.
	Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(M);
	if (cod.rank() == 0) {
		return NA_REAL;
	}

	Eigen::VectorXd x = cod.solve(b);
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
