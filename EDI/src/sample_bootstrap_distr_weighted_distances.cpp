#include <Rcpp.h>
#include <R_ext/Random.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_bootstrapped_weighted_sqd_distances_cpp(
	NumericMatrix X_all_scaled_col_subset,
	NumericVector covariate_weights,
	int t, // self$t
	int B) { // private$other_params$num_boot

	int d = covariate_weights.size();
	NumericVector bootstrapped_weighted_sqd_distances(B);

	GetRNGstate();

	for (int b = 0; b < B; ++b) {
		// sample two distinct indices from 0 to t-1
		int i1 = (int)(unif_rand() * t);
		int i2 = (int)(unif_rand() * t);
		if (t > 1) {
			while (i1 == i2) {
				i2 = (int)(unif_rand() * t);
			}
		}

		// compute delta and weighted sum
		double sqd_weighted_sum = 0.0;
		for (int j = 0; j < d; ++j) {
			double delta = X_all_scaled_col_subset(i1, j) - X_all_scaled_col_subset(i2, j);
			sqd_weighted_sum += delta * delta * covariate_weights[j];
		}

		bootstrapped_weighted_sqd_distances[b] = sqd_weighted_sum;
	}

	PutRNGstate();

	return bootstrapped_weighted_sqd_distances;
}
