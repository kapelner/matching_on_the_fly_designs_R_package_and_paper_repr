#include <Rcpp.h>
#include <random>
#include <utility>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_bootstrapped_weighted_sqd_distances_cpp(
    NumericMatrix X_all_scaled_col_subset,
    NumericVector covariate_weights,
    int t, // self$t
    int B) { // private$other_params$num_boot

  int d = covariate_weights.size();
  NumericVector bootstrapped_weighted_sqd_distances(B);

  for (int b = 0; b < B; ++b) {
    // sample two distinct indices from 0 to t-1
    IntegerVector idx = Rcpp::sample(t, 2, false) - 1; // adjust for C++ 0-indexing

    // extract the rows
    NumericVector x1 = X_all_scaled_col_subset(idx[0], _);
    NumericVector x2 = X_all_scaled_col_subset(idx[1], _);

    // compute delta
    double sqd_weighted_sum = 0.0;
    for (int j = 0; j < d; ++j) {
      double delta = x1[j] - x2[j];
      sqd_weighted_sum += delta * delta * covariate_weights[j];
    }

    bootstrapped_weighted_sqd_distances[b] = sqd_weighted_sum;
  }

  return bootstrapped_weighted_sqd_distances;
}
