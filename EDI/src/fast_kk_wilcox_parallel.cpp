#include <RcppEigen.h>
#include <algorithm>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

inline double logit_cpp(double p) {
    if (p <= 0.0) return -745.0; // Approx log(.Machine$double.eps)
    if (p >= 1.0) return 745.0;
    return std::log(p / (1.0 - p));
}

inline double inv_logit_cpp(double x) {
    if (x < -700) return 0.0;
    if (x > 700) return 1.0;
    return 1.0 / (1.0 + std::exp(-x));
}

inline double apply_shift(double y_val, double delta, int transform_code) {
    if (transform_code == 1) { // log
        return y_val * std::exp(delta);
    }
    if (transform_code == 2) { // logit
        return inv_logit_cpp(logit_cpp(y_val) + delta);
    }
    if (transform_code == 3) { // log1p
        return (y_val + 1.0) * std::exp(delta) - 1.0;
    }
    return y_val + delta;
}

} // namespace

//' Fast KK Wilcoxon Statistic for Multiple Permutations
//'
//' @param y Numeric response vector.
//' @param w_mat Integer matrix of permuted treatment assignments (n x r).
//' @param m_mat Integer matrix of match indicators (n x r).
//' @param delta Null treatment effect shift.
//' @param transform_code Integer code for response transformation.
//' @param is_fixed_matching Logical flag for fixed matching designs.
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of KK Wilcoxon statistics.
// [[Rcpp::export]]
NumericVector compute_kk_wilcox_distr_parallel_cpp(
    const NumericVector& y,
    const IntegerMatrix& w_mat,
    const IntegerMatrix& m_mat,
    double delta,
    int transform_code,
    bool is_fixed_matching,
    int num_cores) {

  int n = y.size();
  int nsim = w_mat.cols();
  std::vector<double> results_vec(nsim);
  
  const double* y_ptr = y.begin();
  const int* w_ptr = w_mat.begin();
  const int* m_ptr = m_mat.begin();
  double* res_ptr = results_vec.data();

  // Pre-calculate shifted responses for treatment group once
  std::vector<double> y_shifted(n);
  if (delta != 0.0) {
      for (int i = 0; i < n; ++i) {
          if (std::isfinite(y_ptr[i])) y_shifted[i] = apply_shift(y_ptr[i], delta, transform_code);
          else y_shifted[i] = NA_REAL;
      }
  }

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

  if (is_fixed_matching && nsim > 0) {
    // ULTRA-FAST PATH
    std::vector<int> res_idx;
    int max_match = 0;
    for (int i = 0; i < n; ++i) if (m_ptr[i] > max_match) max_match = m_ptr[i];
    
    std::vector<int> match_T_row(max_match, -1), match_C_row(max_match, -1);
    for (int i = 0; i < n; ++i) {
      int m = m_ptr[i];
      if (m > 0) {
        if (match_T_row[m-1] == -1) match_T_row[m-1] = i;
        else match_C_row[m-1] = i;
      } else {
        res_idx.push_back(i);
      }
    }
    
    std::vector<std::pair<int, int>> pairs;
    for (int m = 0; m < max_match; ++m) {
      if (match_T_row[m] != -1 && match_C_row[m] != -1) pairs.push_back({match_T_row[m], match_C_row[m]});
    }

    // Pre-rank reservoir if delta=0
    Eigen::VectorXd res_ranks(res_idx.size());
    if (delta == 0 && !res_idx.empty()) {
      int nr = res_idx.size();
      std::vector<std::pair<double, int>> val_idx(nr);
      for (int i = 0; i < nr; ++i) val_idx[i] = std::make_pair(y_ptr[res_idx[i]], i);
      std::sort(val_idx.begin(), val_idx.end());
      for (int i = 0; i < nr; ) {
        int j = i;
        while (j < nr && val_idx[j].first == val_idx[i].first) j++;
        double avg_rank = (i + j + 1.0) / 2.0;
        for (int k = i; k < j; ++k) res_ranks[val_idx[k].second] = avg_rank;
        i = j;
      }
    }

    #pragma omp parallel for schedule(static)
    for (int b = 0; b < nsim; ++b) {
      const int* w_col = w_ptr + (size_t)b * n;
      double total_stat = 0;
      int n_components = 0;

      if (!pairs.empty()) {
        int m_pairs = pairs.size();
        std::vector<double> abs_diffs(m_pairs);
        std::vector<int> signs(m_pairs);
        bool all_zero = true;
        for (int i = 0; i < m_pairs; ++i) {
          int idx1 = pairs[i].first;
          int idx2 = pairs[i].second;
          double y1 = (w_col[idx1] == 1 && delta != 0.0) ? y_shifted[idx1] : y_ptr[idx1];
          double y2 = (w_col[idx2] == 1 && delta != 0.0) ? y_shifted[idx2] : y_ptr[idx2];
          double diff = (w_col[idx1] == 1) ? (y1 - y2) : (y2 - y1);
          abs_diffs[i] = std::abs(diff);
          signs[i] = (diff > 0) ? 1 : (diff < 0 ? -1 : 0);
          if (diff != 0) all_zero = false;
        }
        if (!all_zero) {
          std::vector<std::pair<double, int>> d_val_idx(m_pairs);
          for (int i = 0; i < m_pairs; ++i) d_val_idx[i] = std::make_pair(abs_diffs[i], i);
          std::sort(d_val_idx.begin(), d_val_idx.end());
          double W_plus = 0;
          for (int i = 0; i < m_pairs; ) {
            int j = i;
            while (j < m_pairs && d_val_idx[j].first == d_val_idx[i].first) j++;
            double avg_rank = (i + j + 1.0) / 2.0;
            for (int k = i; k < j; ++k) if (signs[d_val_idx[k].second] > 0) W_plus += avg_rank;
            else if (signs[d_val_idx[k].second] == 0) W_plus += 0.5 * avg_rank;
            i = j;
          }
          total_stat += (W_plus - m_pairs*(m_pairs+1.0)/4.0) / std::sqrt(m_pairs*(m_pairs+1.0)*(2.0*m_pairs+1.0)/24.0);
          n_components++;
        }
      }

      if (!res_idx.empty()) {
        int nr = res_idx.size();
        int nRT = 0;
        for (int i : res_idx) if (w_col[i] == 1) nRT++;
        int nRC = nr - nRT;
        if (nRT > 0 && nRC > 0) {
          double W_r = 0;
          if (delta == 0) {
            for (int i = 0; i < nr; ++i) if (w_col[res_idx[i]] == 1) W_r += res_ranks[i];
          } else {
            std::vector<std::pair<double, int>> r_val_idx(nr);
            for (int i = 0; i < nr; ++i) {
              int idx = res_idx[i];
              double y_val = (w_col[idx] == 1) ? y_shifted[idx] : y_ptr[idx];
              r_val_idx[i] = std::make_pair(y_val, i);
            }
            std::sort(r_val_idx.begin(), r_val_idx.end());
            for (int i = 0; i < nr; ) {
              int j = i;
              while (j < nr && r_val_idx[j].first == r_val_idx[i].first) j++;
              double avg_rank = (i + j + 1.0) / 2.0;
              for (int k = i; k < j; ++k) if (w_col[res_idx[r_val_idx[k].second]] == 1) W_r += avg_rank;
              i = j;
            }
          }
          total_stat += (W_r - nRT*(nr+1.0)/2.0) / std::sqrt(nRT*nRC*(nr+1.0)/12.0);
          n_components++;
        }
      }
      res_ptr[b] = (n_components > 0) ? total_stat : NA_REAL;
    }
    return wrap(results_vec);
  }

#pragma omp parallel for schedule(static)
  for (int b = 0; b < nsim; ++b) {
    const int* w_col = w_ptr + (size_t)b * n;
    const int* m_col = m_ptr + (size_t)b * n;
    
    std::vector<double> diffs;
    std::vector<double> y_r;
    std::vector<int> w_r;
    int max_match = 0;
    for (int i = 0; i < n; ++i) if (m_col[i] > max_match) max_match = m_col[i];
    if (max_match > 0) {
      std::vector<double> match_T(max_match, 0.0), match_C(max_match, 0.0);
      std::vector<bool> has_T(max_match, false), has_C(max_match, false);
      for (int i = 0; i < n; ++i) {
        double y_val = (w_col[i] == 1 && delta != 0.0) ? y_shifted[i] : y_ptr[i];
        if (m_col[i] > 0) {
          int m = m_col[i];
          if (w_col[i] == 1) { match_T[m-1] = y_val; has_T[m-1] = true; }
          else { match_C[m-1] = y_val; has_C[m-1] = true; }
        } else {
          y_r.push_back(y_val);
          w_r.push_back(w_col[i]);
        }
      }
      for (int m = 0; m < max_match; ++m) if (has_T[m] && has_C[m]) diffs.push_back(match_T[m] - match_C[m]);
    } else {
      for (int i = 0; i < n; ++i) {
        double y_val = (w_col[i] == 1 && delta != 0.0) ? y_shifted[i] : y_ptr[i];
        y_r.push_back(y_val);
        w_r.push_back(w_col[i]);
      }
    }
    double total_stat = 0;
    int n_comp = 0;
    if (!diffs.empty()) {
      int m_pairs = diffs.size();
      std::vector<double> abs_diffs(m_pairs);
      std::vector<int> signs(m_pairs);
      bool all_zero = true;
      for (int i = 0; i < m_pairs; ++i) {
        abs_diffs[i] = std::abs(diffs[i]);
        signs[i] = (diffs[i] > 0) ? 1 : (diffs[i] < 0 ? -1 : 0);
        if (diffs[i] != 0) all_zero = false;
      }
      if (!all_zero) {
        std::vector<std::pair<double, int>> d_val_idx(m_pairs);
        for (int i = 0; i < m_pairs; ++i) d_val_idx[i] = std::make_pair(abs_diffs[i], i);
        std::sort(d_val_idx.begin(), d_val_idx.end());
        double W_plus = 0;
        for (int i = 0; i < m_pairs; ) {
          int j = i;
          while (j < m_pairs && d_val_idx[j].first == d_val_idx[i].first) j++;
          double avg_rank = (i + j + 1.0) / 2.0;
          for (int k = i; k < j; ++k) {
            if (signs[d_val_idx[k].second] > 0) W_plus += avg_rank;
            else if (signs[d_val_idx[k].second] == 0) W_plus += 0.5 * avg_rank;
          }
          i = j;
        }
        total_stat += (W_plus - m_pairs*(m_pairs+1.0)/4.0) / std::sqrt(m_pairs*(m_pairs+1.0)*(2.0*m_pairs+1.0)/24.0);
        n_comp++;
      }
    }
    if (!y_r.empty()) {
      int nr = y_r.size();
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
          for (int k = i; k < j; ++k) if (w_r[val_idx[k].second] == 1) W_r += avg_rank;
          i = j;
        }
        total_stat += (W_r - nRT*(nr+1.0)/2.0) / std::sqrt(nRT*nRC*(nr+1.0)/12.0);
        n_comp++;
      }
    }
    res_ptr[b] = (n_comp > 0) ? total_stat : NA_REAL;
  }
  return wrap(results_vec);
}
