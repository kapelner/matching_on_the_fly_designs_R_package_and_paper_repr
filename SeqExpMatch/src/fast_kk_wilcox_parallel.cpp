#include <RcppEigen.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' Fast KK Wilcoxon Statistic for Multiple Permutations
// [[Rcpp::export]]
NumericVector compute_kk_wilcox_distr_parallel_cpp(
    const Eigen::VectorXd& y,
    const Eigen::MatrixXi& w_mat,
    const Eigen::MatrixXi& match_indic_mat,
    int num_cores) {

  int n = y.size();
  int nsim = w_mat.cols();
  NumericVector results(nsim);

#ifdef _OPENMP
  if (num_cores > 1) {
    omp_set_num_threads(num_cores);
  }
#endif

#pragma omp parallel for schedule(static)
  for (int b = 0; b < nsim; ++b) {
    // 1. Extract matched pairs and reservoir for this permutation
    std::vector<double> diffs;
    std::vector<double> y_r;
    std::vector<int> w_r;
    
    // First pass: find max match index
    int max_match = 0;
    for (int i = 0; i < n; ++i) {
      if (match_indic_mat(i, b) > max_match) max_match = match_indic_mat(i, b);
    }
    
    if (max_match > 0) {
      std::vector<double> match_T(max_match, 0.0);
      std::vector<double> match_C(max_match, 0.0);
      std::vector<bool> has_T(max_match, false);
      std::vector<bool> has_C(max_match, false);
      
      for (int i = 0; i < n; ++i) {
        int m = match_indic_mat(i, b);
        if (m > 0) {
          if (w_mat(i, b) == 1) { match_T[m-1] = y[i]; has_T[m-1] = true; }
          else { match_C[m-1] = y[i]; has_C[m-1] = true; }
        } else {
          y_r.push_back(y[i]);
          w_r.push_back(w_mat(i, b));
        }
      }
      
      for (int m = 0; m < max_match; ++m) {
        if (has_T[m] && has_C[m]) {
          diffs.push_back(match_T[m] - match_C[m]);
        }
      }
    } else {
      for (int i = 0; i < n; ++i) {
        y_r.push_back(y[i]);
        w_r.push_back(w_mat(i, b));
      }
    }

    double total_stat = 0;
    int n_components = 0;

    // 2. Signed-Rank for matched pairs
    int m_pairs = diffs.size();
    if (m_pairs > 0) {
      std::vector<double> abs_diffs(m_pairs);
      std::vector<int> signs(m_pairs);
      bool all_zero = true;
      for (int i = 0; i < m_pairs; ++i) {
        abs_diffs[i] = std::abs(diffs[i]);
        if (diffs[i] > 0) signs[i] = 1;
        else if (diffs[i] < 0) signs[i] = -1;
        else signs[i] = 0;
        if (diffs[i] != 0) all_zero = false;
      }
      
      if (!all_zero) {
        // Compute ranks of abs_diffs
        std::vector<std::pair<double, int>> val_idx(m_pairs);
        for (int i = 0; i < m_pairs; ++i) val_idx[i] = std::make_pair(abs_diffs[i], i);
        std::sort(val_idx.begin(), val_idx.end());
        
        double W_plus = 0;
        for (int i = 0; i < m_pairs; ) {
          int j = i;
          while (j < m_pairs && val_idx[j].first == val_idx[i].first) j++;
          double avg_rank = (i + j + 1.0) / 2.0;
          for (int k = i; k < j; ++k) {
            if (signs[val_idx[k].second] > 0) W_plus += avg_rank;
            else if (signs[val_idx[k].second] == 0) W_plus += 0.5 * avg_rank;
          }
          i = j;
        }
        
        double E_W = m_pairs * (m_pairs + 1.0) / 4.0;
        double V_W = m_pairs * (m_pairs + 1.0) * (2.0 * m_pairs + 1.0) / 24.0;
        total_stat += (W_plus - E_W) / std::sqrt(V_W);
        n_components++;
      }
    }

    // 3. Rank-Sum for reservoir
    int nr = y_r.size();
    if (nr > 0) {
      int nRT = 0;
      for (int w : w_r) if (w == 1) nRT++;
      int nRC = nr - nRT;
      
      if (nRT > 0 && nRC > 0) {
        std::vector<std::pair<double, int>> val_idx(nr);
        for (int i = 0; i < nr; ++i) val_idx[i] = std::make_pair(y_r[i], i);
        std::sort(val_idx.begin(), val_idx.end());
        
        double W_r = 0;
        for (int i = 0; i < nr; ) {
          int j = i;
          while (j < nr && val_idx[j].first == val_idx[i].first) j++;
          double avg_rank = (i + j + 1.0) / 2.0;
          for (int k = i; k < j; ++k) {
            if (w_r[val_idx[k].second] == 1) W_r += avg_rank;
          }
          i = j;
        }
        
        double E_W = nRT * (nr + 1.0) / 2.0; // Standard form for Wilcoxon rank sum
        // Note: wilcox.test uses U statistic (W - nRT*(nRT+1)/2). 
        // Standardized version is (W - E_W)/sqrt(V_W).
        double V_W = nRT * nRC * (nr + 1.0) / 12.0;
        total_stat += (W_r - E_W) / std::sqrt(V_W);
        n_components++;
      }
    }

    results[b] = (n_components > 0) ? total_stat : NA_REAL;
  }

  return results;
}
