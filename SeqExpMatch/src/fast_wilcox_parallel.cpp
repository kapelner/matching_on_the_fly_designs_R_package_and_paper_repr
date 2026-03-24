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

double compute_single_wilcox_stat(int n, const double* y_ptr, const int* w_vec, double delta) {
    std::vector<std::pair<double, int>> val_idx(n);
    for (int i = 0; i < n; ++i) {
        val_idx[i] = std::make_pair(y_ptr[i] + (w_vec[i] == 1 ? delta : 0.0), i);
    }
    std::sort(val_idx.begin(), val_idx.end());

    double sum_T = 0;
    int n_T = 0;
    double rank_sum_total = 0;
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && val_idx[j].first == val_idx[i].first) j++;
        double avg_rank = (i + j + 1.0) / 2.0;
        for (int k = i; k < j; ++k) {
            int idx = val_idx[k].second;
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

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

  if (delta == 0) {
    // Pre-rank once and reuse across all permutations
    std::vector<std::pair<double, int>> val_idx(n);
    for (int i = 0; i < n; ++i) val_idx[i] = std::make_pair(y_ptr[i], i);
    std::sort(val_idx.begin(), val_idx.end());
    std::vector<double> ranks(n);
    double total_rank_sum = 0;
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && val_idx[j].first == val_idx[i].first) j++;
        double avg_rank = (i + j + 1.0) / 2.0;
        for (int k = i; k < j; ++k) {
            ranks[val_idx[k].second] = avg_rank;
            total_rank_sum += avg_rank;
        }
        i = j;
    }
    #pragma omp parallel for schedule(static)
    for (int b = 0; b < nsim; ++b) {
        const int* w_col = w_ptr + (size_t)b * n;
        double sum_T = 0;
        int n_T = 0;
        for (int i = 0; i < n; ++i) {
            if (w_col[i] == 1) { sum_T += ranks[i]; n_T++; }
        }
        int n_C = n - n_T;
        results_vec[b] = (n_T == 0 || n_C == 0) ? NA_REAL :
            (sum_T / n_T) - ((total_rank_sum - sum_T) / n_C);
    }
  } else {
    #pragma omp parallel for schedule(static)
    for (int b = 0; b < nsim; ++b) {
        const int* w_col = w_ptr + (size_t)b * n;
        results_vec[b] = compute_single_wilcox_stat(n, y_ptr, w_col, delta);
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

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

  std::vector<const int*> w_ptrs(nsim);
  for (int b = 0; b < nsim; ++b) {
    List p = permutations[b];
    w_ptrs[b] = as<IntegerVector>(p["w"]).begin();
  }

  // Pre-rank if delta=0
  if (delta == 0) {
    std::vector<std::pair<double, int>> val_idx(n);
    for (int i = 0; i < n; ++i) val_idx[i] = std::make_pair(y_ptr[i], i);
    std::sort(val_idx.begin(), val_idx.end());
    std::vector<double> ranks(n);
    double total_rank_sum = 0;
    for (int i = 0; i < n; ) {
        int j = i;
        while (j < n && val_idx[j].first == val_idx[i].first) j++;
        double avg_rank = (i + j + 1.0) / 2.0;
        for (int k = i; k < j; ++k) {
            ranks[val_idx[k].second] = avg_rank;
            total_rank_sum += avg_rank;
        }
        i = j;
    }

    #pragma omp parallel for schedule(static)
    for (int b = 0; b < nsim; ++b) {
        double sum_T = 0;
        int n_T = 0;
        for (int i = 0; i < n; ++i) {
            if (w_ptrs[b][i] == 1) {
                sum_T += ranks[i];
                n_T++;
            }
        }
        int n_C = n - n_T;
        if (n_T == 0 || n_C == 0) results_vec[b] = NA_REAL;
        else results_vec[b] = (sum_T / n_T) - ((total_rank_sum - sum_T) / n_C);
    }
  } else {
    #pragma omp parallel for schedule(static)
    for (int b = 0; b < nsim; ++b) {
        results_vec[b] = compute_single_wilcox_stat(n, y_ptr, w_ptrs[b], delta);
    }
  }

  return wrap(results_vec);
}
