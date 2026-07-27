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

	const int n_rows = X_all_scaled_col_subset.nrow();
	const double* x_ptr = x_new.begin();
	const double* X_ptr = X_all_scaled_col_subset.begin();
	const int* idx_ptr = reservoir_indices.begin();
	const double* w_ptr = covariate_weights.begin();
	double* out_ptr = weighted_sqd_distances.begin();

	for (int j = 0; j < d; ++j) {
		const double xj = x_ptr[j];
		const double wj = w_ptr[j];
		const double* X_col = X_ptr + static_cast<size_t>(j) * n_rows;
#pragma omp simd
		for (int r = 0; r < n; ++r) {
			const int row_idx = idx_ptr[r] - 1; // 1-based to 0-based
			const double diff = xj - X_col[row_idx];
			out_ptr[r] += diff * diff * wj;
		}
	}

	return weighted_sqd_distances;
}
