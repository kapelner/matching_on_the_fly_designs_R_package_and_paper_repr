#include "_helper_functions.h"
#include <unordered_map>
using namespace Rcpp;

// Internal: iterate strata_m, locate each pair (exactly 2 members), check
// discordance, and fill pre-allocated output buffers.  Returns nd (number of
// discordant pairs written).  Caller pre-allocates y_01, t_diffs, X_diffs
// with at least n/2 rows.
static int find_discordant_pairs(
	const Eigen::VectorXd& y_m,
	const Eigen::VectorXd& w_m,
	const Eigen::MatrixXd& X_m,
	const Rcpp::IntegerVector& strata_m,
	Eigen::VectorXd& y_01,
	Eigen::VectorXd& t_diffs,
	Eigen::MatrixXd& X_diffs
) {
	const int n = (int)y_m.size();
	const int p = (int)X_m.cols();

	// Map stratum ID → index of first member seen; -1 once pair is complete.
	std::unordered_map<int,int> first_idx;
	first_idx.reserve(n);

	int nd = 0;
	for (int i = 0; i < n; i++) {
		int s = strata_m[i];
		auto it = first_idx.find(s);
		if (it == first_idx.end()) {
			first_idx[s] = i;                   // first member of this stratum
		} else if (it->second >= 0) {
			int i1 = it->second, i2 = i;
			it->second = -1;                    // mark stratum complete
			if (y_m[i1] == y_m[i2]) continue;  // concordant pair — skip
			double yd   = y_m[i1] - y_m[i2];
			y_01[nd]    = (yd < 0.0) ? 0.0 : 1.0;
			t_diffs[nd] = w_m[i1] - w_m[i2];
			if (p > 0) X_diffs.row(nd) = X_m.row(i1) - X_m.row(i2);
			nd++;
		}
		// 3rd+ member in a stratum: skip (should not occur in KK designs)
	}
	return nd;
}

// Collect discordant matched-pair differences for conditional logistic
// regression.  Used by clogit_helper (pairs only, no reservoir).
//
// Returns: y_01, t_diffs, X_diffs (each trimmed to nd rows), nd.
//
// [[Rcpp::export]]
List collect_discordant_pairs_cpp(
	const Eigen::VectorXd&     y_m,
	const Eigen::VectorXd&     w_m,
	const Eigen::MatrixXd&     X_m,
	const Rcpp::IntegerVector& strata_m
) {
	const int n_max = (int)y_m.size() / 2 + 1;
	const int p     = (int)X_m.cols();

	Eigen::VectorXd y_01(n_max);
	Eigen::VectorXd t_diffs(n_max);
	Eigen::MatrixXd X_diffs(n_max, p);

	int nd = find_discordant_pairs(y_m, w_m, X_m, strata_m, y_01, t_diffs, X_diffs);

	return List::create(
		Named("y_01")    = y_01.head(nd),
		Named("t_diffs") = t_diffs.head(nd),
		Named("X_diffs") = X_diffs.topRows(nd),
		Named("nd")      = nd
	);
}

// Build the combined clogit design matrix for the KK combined-likelihood
// estimator (matched pairs + reservoir, both present).
//
// Column layout: [beta_0 | beta_T | beta_xs (p cols)]
//   Pair rows (0..nd-1):      [0 | t_diff_k | X_diff_k]
//   Reservoir rows (nd..end): [1 | w_r[i]   | X_r[i,] ]
//
// beta_0 = 0 for pair rows encodes the conditional likelihood (intercept is
// eliminated by conditioning on the pair sum).
// Fitting logistic on the stacked dataset maximises L_cond(pairs)+L_marg(res).
//
// Returns: X_comb ((nd+nR) x (p+2)), y_comb (nd+nR), nd.
//
// [[Rcpp::export]]
List build_kk_combined_clogit_design_cpp(
	const Eigen::VectorXd&     y_m,
	const Eigen::VectorXd&     w_m,
	const Eigen::MatrixXd&     X_m,
	const Rcpp::IntegerVector& strata_m,
	const Eigen::VectorXd&     y_r,
	const Eigen::VectorXd&     w_r,
	const Eigen::MatrixXd&     X_r
) {
	const int nR    = (int)y_r.size();
	const int p     = (int)X_m.cols();
	const int n_max = (int)y_m.size() / 2 + 1;

	// Temporary pair-row buffers (worst-case allocation)
	Eigen::VectorXd y_01(n_max);
	Eigen::VectorXd t_diffs(n_max);
	Eigen::MatrixXd X_diffs(n_max, p);

	int nd = find_discordant_pairs(y_m, w_m, X_m, strata_m, y_01, t_diffs, X_diffs);

	// Build combined design [beta_0 | beta_T | beta_xs]
	Eigen::MatrixXd X_comb(nd + nR, p + 2);
	Eigen::VectorXd y_comb(nd + nR);

	// ---- Pair rows -------------------------------------------------------
	X_comb.block(0, 0, nd, 1).setZero();
	X_comb.block(0, 1, nd, 1) = t_diffs.head(nd);
	if (p > 0) X_comb.block(0, 2, nd, p) = X_diffs.topRows(nd);
	y_comb.head(nd) = y_01.head(nd);

	// ---- Reservoir rows --------------------------------------------------
	X_comb.block(nd, 0, nR, 1).setOnes();
	X_comb.block(nd, 1, nR, 1) = w_r;
	if (p > 0) X_comb.block(nd, 2, nR, p) = X_r;
	y_comb.tail(nR) = y_r;

	return List::create(
		Named("X_comb") = X_comb,
		Named("y_comb") = y_comb,
		Named("nd")     = nd
	);
}
