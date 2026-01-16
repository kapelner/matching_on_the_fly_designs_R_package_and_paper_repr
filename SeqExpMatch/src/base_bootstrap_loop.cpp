#include <RcppEigen.h>
// #ifdef _OPENMP
// #include <omp.h>
// #endif

// [[Rcpp::depends(RcppEigen)]]
// // [[Rcpp::plugins(openmp)]] // Temporarily disable OpenMP plugin for debugging

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector base_bootstrap_loop_cpp(
  const IntegerMatrix& indices,
  const NumericVector& y,
  const NumericVector& dead,
  const NumericMatrix& X,
  const NumericVector& w,
  Function duplicate_inference_fn,
  Function compute_estimate_fn,
  int num_cores = 1
) {
  int B = indices.nrow();
  int n = indices.ncol();
  int p = X.ncol();

  NumericVector estimates(B);

  Rcout << "DEBUG: Entering base_bootstrap_loop_cpp with B = " << B << ", num_cores = " << num_cores << std::endl;

// #ifdef _OPENMP
//   if (num_cores > 1) {
//     omp_set_num_threads(num_cores);
//   }
// #endif

// #pragma omp parallel if(num_cores > 1) // Temporarily disable OpenMP for debugging
  {
    SEXP thread_inf_obj_sexp;
    bool inf_obj_valid = false;
    try {
      thread_inf_obj_sexp = duplicate_inference_fn();
      if (thread_inf_obj_sexp == R_NilValue) {
        Rcout << "DEBUG: duplicate_inference_fn returned R_NilValue" << std::endl;
      } else if (TYPEOF(thread_inf_obj_sexp) != ENVSXP) {
        // R6 objects are environments (ENVSXP), not S4 objects
        Rcout << "DEBUG: duplicate_inference_fn returned non-environment object (type " << TYPEOF(thread_inf_obj_sexp) << ")" << std::endl;
      } else {
        Rcout << "DEBUG: duplicate_inference_fn returned valid R6 object (environment)" << std::endl;
        inf_obj_valid = true;
      }
    } catch (const Rcpp::exception& e) {
      Rcout << "DEBUG: Rcpp::exception in duplicate_inference_fn: " << e.what() << std::endl;
    } catch (...) {
      Rcout << "DEBUG: Unknown exception in duplicate_inference_fn" << std::endl;
    }

    // Use RObject instead of S4 since R6 objects are environments
    RObject thread_inf_obj = R_NilValue;
    if (inf_obj_valid) {
      thread_inf_obj = thread_inf_obj_sexp;
    }


// #pragma omp for schedule(dynamic) // Temporarily disable OpenMP for debugging
    for (int b = 0; b < B; ++b) {
      // Get bootstrap indices for this sample
      IntegerVector i_b = indices(b, _);

      // Create resampled vectors/matrices
      NumericVector y_b(n);
      NumericVector dead_b(n);
      NumericVector w_b(n);
      NumericMatrix X_b(n, p);

      for (int i = 0; i < n; ++i) {
        int idx = i_b[i] - 1;  // R indices are 1-based, C++ are 0-based
        y_b[i] = y[idx];
        dead_b[i] = dead[idx];
        w_b[i] = w[idx];
        for (int j = 0; j < p; ++j) {
          X_b(i, j) = X(idx, j);
        }
      }

      // Preserve column names from original X
      if (X.hasAttribute("dimnames")) {
        List dimnames = X.attr("dimnames");
        if (dimnames.size() >= 2 && !Rf_isNull(dimnames[1])) {
          List new_dimnames = List::create(R_NilValue, dimnames[1]);
          X_b.attr("dimnames") = new_dimnames;
        }
      }

      // Call R function to compute estimate using resampled data
      // Each thread uses its own inference object (thread_inf_obj)
      List bootstrap_data;
      if (inf_obj_valid) {
        bootstrap_data = List::create(
          Named("y") = y_b,
          Named("dead") = dead_b,
          Named("X") = X_b,
          Named("w") = w_b,
          Named("inf_obj") = thread_inf_obj
        );
      } else {
        bootstrap_data = List::create(
          Named("y") = y_b,
          Named("dead") = dead_b,
          Named("X") = X_b,
          Named("w") = w_b
        );
      }


      NumericVector result;
      try {
        Rcout << "DEBUG: Calling compute_estimate_fn for b = " << b << std::endl;
        result = compute_estimate_fn(bootstrap_data);
        Rcout << "DEBUG: compute_estimate_fn returned, result.length() = " << result.length() << std::endl;
        if (result.length() == 0) {
            Rcout << "DEBUG: compute_estimate_fn returned empty vector for b = " << b << std::endl;
            estimates[b] = NA_REAL;
        } else if (result.length() > 1) {
            Rcout << "DEBUG: compute_estimate_fn returned vector of length > 1 for b = " << b << std::endl;
            estimates[b] = result[0]; // Take the first element, or handle as appropriate
        } else {
            Rcout << "DEBUG: result[0] = " << result[0] << std::endl;
            estimates[b] = result[0];
        }
      } catch (const Rcpp::exception& e) {
        Rcout << "DEBUG: Rcpp::exception in compute_estimate_fn for b = " << b << ": " << e.what() << std::endl;
        estimates[b] = NA_REAL;
      } catch (...) {
        Rcout << "DEBUG: Unknown exception in compute_estimate_fn for b = " << b << std::endl;
        estimates[b] = NA_REAL;
      }
    }
  }

  Rcout << "DEBUG: Exiting base_bootstrap_loop_cpp." << std::endl;
  return estimates;
}
