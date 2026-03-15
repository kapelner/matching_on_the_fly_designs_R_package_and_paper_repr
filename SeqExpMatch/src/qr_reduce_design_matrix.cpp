// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

static std::vector<int> sort_unique_ints(std::vector<int> x) {
	std::sort(x.begin(), x.end());
	x.erase(std::unique(x.begin(), x.end()), x.end());
	return x;
}

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
		required.push_back(j - 1);
	}
	required = sort_unique_ints(required);

	ColPivHouseholderQR<MatrixXd> qr_init(X);
	const int target_rank = qr_init.rank();
	const auto& pivot = qr_init.colsPermutation().indices();

	std::vector<int> candidate_order = required;
	candidate_order.reserve(p);
	for (int i = 0; i < p; ++i) {
		const int col = pivot(i);
		bool is_required = false;
		for (int r : required) if (r == col) { is_required = true; break; }
		if (!is_required) candidate_order.push_back(col);
	}

	std::vector<int> keep;
	keep.reserve(target_rank);
	for (int col : candidate_order) {
		std::vector<int> trial_keep = keep;
		trial_keep.push_back(col);
		
		MatrixXd X_trial(n, trial_keep.size());
		for (size_t i = 0; i < trial_keep.size(); ++i) {
			X_trial.col(i) = X.col(trial_keep[i]);
		}
		
		const int trial_rank = ColPivHouseholderQR<MatrixXd>(X_trial).rank();
		if (trial_rank > static_cast<int>(keep.size())) {
			keep.push_back(col);
		}
		if (static_cast<int>(keep.size()) >= target_rank) {
			break;
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
