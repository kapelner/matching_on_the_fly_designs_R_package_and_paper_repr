// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
List qr_reduce_full_rank_cpp(const Eigen::MatrixXd& X) {
	const int p = X.cols();
	const int n = X.rows();
	if (p == 0) {
		return List::create(
			Named("X_reduced") = X,
			Named("keep") = IntegerVector(0)
		);
	}

	ColPivHouseholderQR<MatrixXd> qr(X);
	const int rank = qr.rank();
	const auto& pivot = qr.colsPermutation().indices();
	std::vector<int> keep(rank);
	for (int i = 0; i < rank; ++i) {
		keep[i] = pivot(i);
	}
	std::sort(keep.begin(), keep.end());

	MatrixXd X_reduced(n, rank);
	IntegerVector keep_r(rank);
	for (int i = 0; i < rank; ++i) {
		X_reduced.col(i) = X.col(keep[i]);
		keep_r[i] = keep[i] + 1;
	}

	return List::create(
		Named("X_reduced") = X_reduced,
		Named("keep") = keep_r
	);
}

// [[Rcpp::export]]
List qr_reduce_preserve_cols_cpp(const Eigen::MatrixXd& X, IntegerVector required_cols) {
	const int p = X.cols();
	const int n = X.rows();
	if (p == 0) {
		return List::create(
			Named("X_reduced") = X,
			Named("keep") = IntegerVector(0)
		);
	}

	std::vector<int> required;
	required.reserve(required_cols.size());
	for (int j : required_cols) {
		if (j < 1 || j > p) stop("required column index is out of bounds");
		const int col = j - 1;
		bool seen = false;
		for (int kept : required) {
			if (kept == col) {
				seen = true;
				break;
			}
		}
		if (!seen) required.push_back(col);
	}

	std::vector<int> keep_required;
	keep_required.reserve(required.size());
	for (int j : required_cols) {
		const int col = j - 1;
		bool is_required = false;
		for (int r : required) {
			if (r == col) {
				is_required = true;
				break;
			}
		}
		if (!is_required) continue;

		std::vector<int> trial_keep = keep_required;
		trial_keep.push_back(col);
		MatrixXd X_trial(n, trial_keep.size());
		for (size_t i = 0; i < trial_keep.size(); ++i) {
			X_trial.col(i) = X.col(trial_keep[i]);
		}
		const int trial_rank = ColPivHouseholderQR<MatrixXd>(X_trial).rank();
		if (trial_rank > static_cast<int>(keep_required.size())) {
			keep_required.push_back(col);
		}
	}

	std::vector<int> remaining;
	remaining.reserve(p - keep_required.size());
	for (int j = 0; j < p; ++j) {
		bool already_kept = false;
		for (int kept : keep_required) {
			if (kept == j) {
				already_kept = true;
				break;
			}
		}
		if (!already_kept) remaining.push_back(j);
	}

	std::vector<int> keep = keep_required;
	keep.reserve(p);
	if (!remaining.empty()) {
		MatrixXd X_remaining(n, remaining.size());
		for (size_t i = 0; i < remaining.size(); ++i) {
			X_remaining.col(i) = X.col(remaining[i]);
		}

		MatrixXd X_residual = X_remaining;
		if (!keep_required.empty()) {
			MatrixXd X_required(n, keep_required.size());
			for (size_t i = 0; i < keep_required.size(); ++i) {
				X_required.col(i) = X.col(keep_required[i]);
			}
			ColPivHouseholderQR<MatrixXd> qr_required(X_required);
			X_residual -= X_required * qr_required.solve(X_remaining);
		}

		ColPivHouseholderQR<MatrixXd> qr_residual(X_residual);
		const int residual_rank = qr_residual.rank();
		const auto& residual_pivot = qr_residual.colsPermutation().indices();
		for (int i = 0; i < residual_rank; ++i) {
			const int candidate_col = remaining[residual_pivot(i)];
			std::vector<int> trial_keep = keep;
			trial_keep.push_back(candidate_col);
			MatrixXd X_trial(n, trial_keep.size());
			for (size_t j = 0; j < trial_keep.size(); ++j) {
				X_trial.col(j) = X.col(trial_keep[j]);
			}
			const int trial_rank = ColPivHouseholderQR<MatrixXd>(X_trial).rank();
			if (trial_rank > static_cast<int>(keep.size())) {
				keep.push_back(candidate_col);
			}
		}
	}

	std::sort(keep.begin(), keep.end());
	int final_p = keep.size();
	MatrixXd X_reduced(n, final_p);
	IntegerVector keep_r(final_p);
	for (int i = 0; i < final_p; ++i) {
		X_reduced.col(i) = X.col(keep[i]);
		keep_r[i] = keep[i] + 1;
	}

	return List::create(
		Named("X_reduced") = X_reduced,
		Named("keep") = keep_r
	);
}
