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
  const IntegerVector& original_match_indic, // Changed name
  const IntegerVector& i_b,
  int m
);

// [[Rcpp::export]]
NumericVector kk_bootstrap_loop_cpp(
  const IntegerMatrix& indices,
  const NumericVector& y,
  const NumericVector& w,
  const NumericMatrix& X,
  const IntegerVector& original_match_indic, // Added this argument
  int m,
  Function duplicate_inference_fn,
  Function compute_estimate_fn,
  int num_cores = 1
) {
  int B = indices.nrow();
  NumericVector estimates(B);

  // Create one inference object for serial execution
  SEXP thread_inf_obj_sexp = duplicate_inference_fn();
  bool thread_safe_inf_obj_valid = (thread_inf_obj_sexp != R_NilValue);

  // Run serially to avoid R callback issues in OpenMP
  for (int b = 0; b < B; ++b) {
    if (!thread_safe_inf_obj_valid) {
      estimates[b] = NA_REAL;
      continue; // Skip this iteration if the inference object is invalid
    }

    // Get bootstrap indices for this sample
    IntegerVector i_b = indices(b, _);

    // Compute match statistics for this bootstrap sample
    List kk_stats = match_stats_from_indices_cpp(y, w, X, original_match_indic, i_b, m); // Changed this line

    // Add thread-local inference object to the stats
    kk_stats["inf_obj"] = thread_inf_obj_sexp;
    // Also pass a flag indicating if the inf_obj is valid
    kk_stats["thread_safe_inf_obj_valid"] = thread_safe_inf_obj_valid;

    // Call R function to compute estimate using these statistics
    NumericVector result = compute_estimate_fn(kk_stats);

    if (result.length() > 0) {
      estimates[b] = result[0];
    } else {
      estimates[b] = NA_REAL; // Handle empty result from R function
    }
  }

  return estimates;
}
