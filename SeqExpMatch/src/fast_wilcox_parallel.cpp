#include <RcppEigen.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' Fast Wilcoxon Rank Sum Statistic for Multiple Permutations
//'
//' @param y The response vector (already shifted by delta if necessary).
//' @param w_mat Matrix of permutations (n x nsim).
//' @param num_cores Number of cores to use.
//'
//' @return A vector of rank-sum statistics (mean differences in ranks).
// [[Rcpp::export]]
NumericVector compute_wilcox_distr_parallel_cpp(const Eigen::VectorXd& y, const Eigen::MatrixXi& w_mat, int num_cores) {
  int n = y.size();
  int nsim = w_mat.cols();
  NumericVector results(nsim);

  // 1. Compute ranks once (O(n log n))
  std::vector<std::pair<double, int>> val_idx(n);
  for (int i = 0; i < n; ++i) {
    val_idx[i] = std::make_pair(y[i], i);
  }
  std::sort(val_idx.begin(), val_idx.end());

  Eigen::VectorXd ranks(n);
  for (int i = 0; i < n; ) {
    int j = i;
    while (j < n && val_idx[j].first == val_idx[i].first) j++;
    double avg_rank = (i + j + 1.0) / 2.0;
    for (int k = i; k < j; ++k) {
      ranks[val_idx[k].second] = avg_rank;
    }
    i = j;
  }

#ifdef _OPENMP
  if (num_cores > 1) {
    omp_set_num_threads(num_cores);
  }
#endif

#pragma omp parallel for schedule(static)
  for (int b = 0; b < nsim; ++b) {
    double sum_T = 0;
    int n_T = 0;
    for (int i = 0; i < n; ++i) {
      if (w_mat(i, b) == 1) {
        sum_T += ranks[i];
        n_T++;
      }
    }
    int n_C = n - n_T;
    if (n_T == 0 || n_C == 0) {
      results[b] = NA_REAL;
    } else {
      // Return difference in mean ranks
      results[b] = (sum_T / n_T) - ((ranks.sum() - sum_T) / n_C);
    }
  }

  return results;
}
