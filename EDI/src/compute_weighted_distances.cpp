#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_weighted_sqd_distances_cpp(
	NumericVector x_new,
	NumericMatrix X_all_scaled_col_subset,
	IntegerVector reservoir_indices,
	NumericVector covariate_weights) {

	int n = reservoir_indices.size();
	int d = x_new.size();
	NumericVector weighted_sqd_distances(n);

	for (int r = 0; r < n; ++r) {
		int row_idx = reservoir_indices[r] - 1; // 1-based to 0-based
		double sqd_weighted_sum = 0.0;
		for (int j = 0; j < d; ++j) {
			double diff = x_new[j] - X_all_scaled_col_subset(row_idx, j);
			sqd_weighted_sum += diff * diff * covariate_weights[j];
		}
		weighted_sqd_distances[r] = sqd_weighted_sum;
	}

	return weighted_sqd_distances;
}
