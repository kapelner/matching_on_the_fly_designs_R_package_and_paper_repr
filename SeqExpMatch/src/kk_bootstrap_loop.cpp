#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Forward declaration
List match_stats_from_indices_cpp(
  const NumericVector& y,
  const NumericVector& w,
  const NumericMatrix& X,
  const IntegerVector& match_indic_b,
  const IntegerVector& i_b,
  int m
);

// [[Rcpp::export]]
NumericVector kk_bootstrap_loop_cpp(
  const IntegerMatrix& indices,
  const NumericVector& y,
  const NumericVector& w,
  const NumericMatrix& X,
  const IntegerVector& match_indic_b,
  int m,
  Function duplicate_inference_fn,
  Function compute_estimate_fn,
  int num_cores = 1
) {
  int B = indices.nrow();
  NumericVector estimates(B);

#ifdef _OPENMP
  if (num_cores > 1) {
    omp_set_num_threads(num_cores);
  }
#endif

#pragma omp parallel if(num_cores > 1)
  {
    // Each thread creates its own inference object copy for thread safety
    // This is only done once per thread, not once per bootstrap iteration
    SEXP thread_inf_obj = duplicate_inference_fn();

#pragma omp for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
      // Get bootstrap indices for this sample
      IntegerVector i_b = indices(b, _);

      // Compute match statistics for this bootstrap sample
      List kk_stats = match_stats_from_indices_cpp(y, w, X, match_indic_b, i_b, m);

      // Add thread-local inference object to the stats
      kk_stats["inf_obj"] = thread_inf_obj;

      // Call R function to compute estimate using these statistics
      // Each thread uses its own inference object
      NumericVector result = compute_estimate_fn(kk_stats);
      estimates[b] = result[0];
    }
  }

  return estimates;
}
