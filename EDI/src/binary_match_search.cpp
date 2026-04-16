#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix draw_binary_match_assignments_cpp(IntegerMatrix indices_pairs, int n, int r, int num_cores) {
  int num_pairs = indices_pairs.nrow();
  NumericMatrix w_mat(n, r);
  
  // R's RNG is single-threaded. We draw all random values upfront serially
  // to remain thread-safe and maintain R's RNG state consistently.
  NumericVector rand_vals_all = Rcpp::runif(num_pairs * r);
  
  const int* pairs_ptr = indices_pairs.begin();
  const double* rand_ptr = rand_vals_all.begin();
  double* w_ptr = w_mat.begin();

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static)
  for (int j = 0; j < r; j++) {
    const double* rand_col = rand_ptr + (size_t)j * num_pairs;
    double* w_col = w_ptr + (size_t)j * n;
    
    for (int i = 0; i < num_pairs; i++) {
      // indices_pairs has 2 columns. indices_pairs(i, 0) and indices_pairs(i, 1)
      // indices_pairs is stored in column-major order.
      int idx1 = pairs_ptr[i] - 1;                // Col 0, Row i
      int idx2 = pairs_ptr[i + num_pairs] - 1;    // Col 1, Row i
      
      if (rand_col[i] < 0.5) {
        w_col[idx1] = 1.0;
        w_col[idx2] = 0.0;
      } else {
        w_col[idx1] = 0.0;
        w_col[idx2] = 1.0;
      }
    }
  }
  
  return w_mat;
}
