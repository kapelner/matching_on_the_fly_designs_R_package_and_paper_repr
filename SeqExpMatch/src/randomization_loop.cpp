#include <RcppEigen.h>
// #ifdef _OPENMP
// #include <omp.h>
// #endif

// [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::plugins(openmp)]] // Temporarily disable OpenMP for debugging

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector randomization_loop_cpp(
  int nsim_exact_test,
  Function duplicate_design_fn,
  Function duplicate_inference_fn,
  Function run_randomization_iteration_fn,
  int num_cores = 1
) {
  NumericVector estimates(nsim_exact_test);

// #ifdef _OPENMP
//   if (num_cores > 1) {
//     omp_set_num_threads(num_cores);
//   }
// #endif

// #pragma omp parallel if(num_cores > 1) // Temporarily disable OpenMP for debugging
  {
    // Each thread creates its own design and inference object copies
    // This is done once per thread, not once per randomization iteration
    SEXP thread_des_obj = duplicate_design_fn();
    SEXP thread_inf_obj = duplicate_inference_fn();

    // Create a list to pass both objects to the callback
    List thread_objects = List::create(
      Named("design") = thread_des_obj,
      Named("inference") = thread_inf_obj
    );

// #pragma omp for schedule(dynamic) // Temporarily disable OpenMP for debugging
    for (int r = 0; r < nsim_exact_test; ++r) {
      // Call R function to run one randomization iteration
      // Each thread uses its own design and inference objects
      NumericVector result = run_randomization_iteration_fn(thread_objects);
      estimates[r] = result[0];
    }
  }

  return estimates;
}
