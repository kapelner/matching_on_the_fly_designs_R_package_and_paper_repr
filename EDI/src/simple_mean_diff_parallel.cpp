#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_simple_mean_diff_parallel_cpp(
	const NumericVector& y,
	const IntegerMatrix& w_mat,
	double delta,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = w_mat.rows();
	
	// 1. Use a standard C++ vector to avoid R API contention inside the loop
	std::vector<double> results_vec(nsim);
	
	// 2. Extract raw pointers for lock-free access
	const double* y_ptr = y.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static) if(use_parallel)
	for (int b = 0; b < nsim; ++b) {
		const int* w_col = w_ptr + (size_t)b * n;
		double sum_T = 0, sum_C = 0;
		int n_T = 0;

		for (int i = 0; i < n; ++i) {
			const int is_t = (w_col[i] == 1);
			sum_T += is_t * (y_ptr[i] + delta);
			sum_C += (1 - is_t) * y_ptr[i];
			n_T += is_t;
		}

		const int n_C = n - n_T;
		if (n_T == 0 || n_C == 0) {
			res_ptr[b] = NA_REAL;
		} else {
			res_ptr[b] = (sum_T / n_T) - (sum_C / n_C);
		}
	}

	// 3. Wrap back to R type only at the very end (main thread only)
	return wrap(results_vec);
}
