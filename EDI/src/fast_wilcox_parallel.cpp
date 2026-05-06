#include "_helper_functions.h"
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

void sort_order_by_values(std::vector<int>& order, const std::vector<double>& values) {
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (values[a] < values[b]) return true;
        if (values[b] < values[a]) return false;
        return a < b;
    });
}

double compute_single_wilcox_stat(int n, const double* y_ptr, const int* w_vec, double delta,
                                  std::vector<double>& values,
                                  std::vector<int>& order) {
    values.resize(n);
    order.resize(n);
    for (int i = 0; i < n; ++i) {
        values[i] = y_ptr[i] + (w_vec[i] == 1 ? delta : 0.0);
        order[i] = i;
    }
    sort_order_by_values(order, values);

    double sum_T = 0;
    int n_T = 0;
    double rank_sum_total = 0;
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && values[order[j]] == values[order[i]]) j++;
        double avg_rank = (i + j + 1.0) / 2.0;
        for (int k = i; k < j; ++k) {
            int idx = order[k];
            rank_sum_total += avg_rank;
            if (w_vec[idx] == 1) {
                sum_T += avg_rank;
                n_T++;
            }
        }
        i = j;
    }
    int n_C = n - n_T;
    if (n_T == 0 || n_C == 0) return NA_REAL;
    return (sum_T / n_T) - ((rank_sum_total - sum_T) / n_C);
}

} // namespace

// [[Rcpp::export]]
NumericVector compute_wilcox_distr_parallel_cpp(
    const NumericVector& y,
    const IntegerMatrix& w_mat,
    double delta,
    int num_cores) {

  int n = y.size();
  int nsim = w_mat.cols();
  std::vector<double> results_vec(nsim);
  const double* y_ptr = y.begin();
  const int* w_ptr = w_mat.begin();
  const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
  if (use_parallel) omp_set_num_threads(num_cores);
#endif

  if (delta == 0) {
    // Pre-rank once and reuse across all permutations
    std::vector<double> values(n);
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) {
        values[i] = y_ptr[i];
        order[i] = i;
    }
    sort_order_by_values(order, values);
    std::vector<double> ranks(n);
    double total_rank_sum = 0;
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && values[order[j]] == values[order[i]]) j++;
        double avg_rank = (i + j + 1.0) / 2.0;
        for (int k = i; k < j; ++k) {
            ranks[order[k]] = avg_rank;
            total_rank_sum += avg_rank;
        }
        i = j;
    }
    #pragma omp parallel for schedule(static) if(use_parallel)
    for (int b = 0; b < nsim; ++b) {
        const int* w_col = w_ptr + (size_t)b * n;
        double sum_T = 0;
        int n_T = 0;
        for (int i = 0; i < n; ++i) {
            const int is_t = (w_col[i] == 1);
            sum_T += is_t * ranks[i];
            n_T += is_t;
        }
        int n_C = n - n_T;
        results_vec[b] = (n_T == 0 || n_C == 0) ? NA_REAL :
            (sum_T / n_T) - ((total_rank_sum - sum_T) / n_C);
    }
  } else {
    #pragma omp parallel if(use_parallel)
    {
        std::vector<double> values;
        std::vector<int> order;
        values.reserve(n);
        order.reserve(n);

        #pragma omp for schedule(static)
        for (int b = 0; b < nsim; ++b) {
            const int* w_col = w_ptr + (size_t)b * n;
            results_vec[b] = compute_single_wilcox_stat(n, y_ptr, w_col, delta, values, order);
        }
    }
  }
  return wrap(results_vec);
}

// [[Rcpp::export]]
NumericVector compute_wilcox_distr_from_list_parallel_cpp(
    const NumericVector& y,
    const List& permutations,
    double delta,
    int num_cores) {

  int n = y.size();
  int nsim = permutations.size();
  std::vector<double> results_vec(nsim);
  const double* y_ptr = y.begin();
  const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
  if (use_parallel) omp_set_num_threads(num_cores);
#endif

  std::vector<const int*> w_ptrs(nsim);
  for (int b = 0; b < nsim; ++b) {
    List p = permutations[b];
    w_ptrs[b] = as<IntegerVector>(p["w"]).begin();
  }

  // Pre-rank if delta=0
  if (delta == 0) {
    std::vector<double> values(n);
    std::vector<int> order(n);
    for (int i = 0; i < n; ++i) {
        values[i] = y_ptr[i];
        order[i] = i;
    }
    sort_order_by_values(order, values);
    std::vector<double> ranks(n);
    double total_rank_sum = 0;
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && values[order[j]] == values[order[i]]) j++;
        double avg_rank = (i + j + 1.0) / 2.0;
        for (int k = i; k < j; ++k) {
            ranks[order[k]] = avg_rank;
            total_rank_sum += avg_rank;
        }
        i = j;
    }

    #pragma omp parallel for schedule(static) if(use_parallel)
    for (int b = 0; b < nsim; ++b) {
        double sum_T = 0;
        int n_T = 0;
        for (int i = 0; i < n; ++i) {
            const int is_t = (w_ptrs[b][i] == 1);
            sum_T += is_t * ranks[i];
            n_T += is_t;
        }
        int n_C = n - n_T;
        if (n_T == 0 || n_C == 0) results_vec[b] = NA_REAL;
        else results_vec[b] = (sum_T / n_T) - ((total_rank_sum - sum_T) / n_C);
    }
  } else {
    #pragma omp parallel if(use_parallel)
    {
        std::vector<double> values;
        std::vector<int> order;
        values.reserve(n);
        order.reserve(n);

        #pragma omp for schedule(static)
        for (int b = 0; b < nsim; ++b) {
            results_vec[b] = compute_single_wilcox_stat(n, y_ptr, w_ptrs[b], delta, values, order);
        }
    }
  }

  return wrap(results_vec);
}
