#include "_helper_functions.h"
using namespace Rcpp;


//typedef Eigen::Map<Eigen::MatrixXd> MapMat;
//typedef Eigen::Map<Eigen::VectorXd> MapVec;

double compute_diagonal_inverse_entry(const Eigen::MatrixXd& M, int j) {
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
double eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(Eigen::MatrixXd M, int j) {
    return compute_diagonal_inverse_entry(M, j);
}

// [[Rcpp::export]]
Eigen::MatrixXd eigen_Xt_times_diag_w_times_X_cpp(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::VectorXd> w) {
	return X.transpose() * w.asDiagonal() * X;
}

// [[Rcpp::export]]
Rcpp::List likelihood_ratio_test_from_negloglik_cpp(double unrestricted_neg_loglik,
                                                    double null_neg_loglik,
                                                    int df = 1) {
	return likelihood_ratio_test_from_negloglik(unrestricted_neg_loglik, null_neg_loglik, df);
}

// [[Rcpp::export]]
Rcpp::List score_test_from_score_information_cpp(const Eigen::VectorXd& score,
                                                 const Eigen::MatrixXd& information,
                                                 int tested_idx) {
	return score_test_from_score_information(score, information, tested_idx);
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
