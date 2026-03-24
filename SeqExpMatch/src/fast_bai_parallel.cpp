#include <RcppEigen.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

//' Fast Bai Adjusted T Statistic for Multiple Permutations
// [[Rcpp::export]]
NumericVector compute_bai_distr_parallel_cpp(
    const NumericVector& y,
    const IntegerMatrix& w_mat,
    const IntegerMatrix& m_mat,
    double delta,
    const IntegerMatrix& halves_idx,
    bool convex_flag,
    int num_cores) {

  int n = y.size();
  int nsim = w_mat.cols();
  int n_halves = halves_idx.rows();
  std::vector<double> results_vec(nsim);

  const double* y_ptr = y.begin();
  const int* w_ptr = w_mat.begin();
  const int* m_ptr = m_mat.begin();
  const int* h_ptr = halves_idx.begin();
  double* res_ptr = results_vec.data();

#ifdef _OPENMP
  omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static)
  for (int b = 0; b < nsim; ++b) {
    const int* w_col = w_ptr + (size_t)b * n;
    const int* m_col = m_ptr + (size_t)b * n;

    std::vector<double> d_i;
    std::vector<double> y_r;
    std::vector<int> w_r;
    
    int max_match = 0;
    for (int i = 0; i < n; ++i) if (m_col[i] > max_match) max_match = m_col[i];
    
    if (max_match > 0) {
      std::vector<double> match_T(max_match, 0.0), match_C(max_match, 0.0);
      std::vector<bool> has_T(max_match, false), has_C(max_match, false);
      for (int i = 0; i < n; ++i) {
        double y_val = y_ptr[i] + (w_col[i] == 1 ? delta : 0.0);
        int m = m_col[i];
        if (m > 0) {
          if (w_col[i] == 1) { match_T[m-1] = y_val; has_T[m-1] = true; }
          else { match_C[m-1] = y_val; has_C[m-1] = true; }
        } else {
          y_r.push_back(y_val);
          w_r.push_back(w_col[i]);
        }
      }
      for (int m = 0; m < max_match; ++m) if (has_T[m] && has_C[m]) d_i.push_back(match_T[m] - match_C[m]);
    } else {
      for (int i = 0; i < n; ++i) {
        y_r.push_back(y_ptr[i] + (w_col[i] == 1 ? delta : 0.0));
        w_r.push_back(w_col[i]);
      }
    }

    int m_size = d_i.size();
    if (m_size == 0 && y_r.empty()) { res_ptr[b] = NA_REAL; continue; }

    double d_bar = 0;
    if (m_size > 0) {
      for (double d : d_i) d_bar += d;
      d_bar /= m_size;
    }

    double r_bar = 0;
    double ssqR = 0;
    int nRT = 0, nRC = 0;
    if (!y_r.empty()) {
      double sumT = 0, sumC = 0;
      std::vector<double> yT, yC;
      for (size_t i = 0; i < y_r.size(); ++i) {
        if (w_r[i] == 1) { sumT += y_r[i]; nRT++; yT.push_back(y_r[i]); }
        else { sumC += y_r[i]; nRC++; yC.push_back(y_r[i]); }
      }
      if (nRT > 0 && nRC > 0) {
        r_bar = (sumT/nRT) - (sumC/nRC);
        if (nRT > 1 && nRC > 1) {
          double varT = 0, varC = 0;
          double meanT = sumT/nRT, meanC = sumC/nRC;
          for (double val : yT) varT += (val - meanT)*(val - meanT);
          for (double val : yC) varC += (val - meanC)*(val - meanC);
          ssqR = (varT/(nRT-1))/nRT + (varC/(nRC-1))/nRC;
        }
      }
    }

    double bai_var_d_bar = 0;
    if (m_size > 0) {
      double delta_sq = d_bar * d_bar;
      double tau_sq = 0;
      for (double d : d_i) tau_sq += d * d;
      tau_sq /= m_size;
      double lambda_squ = 0;
      if (n_halves > 0) {
        for (int i = 0; i < n_halves; ++i) {
          int idx1 = h_ptr[i] - 1; // halves_idx is IntegerMatrix, halves_idx(i, 0)
          int idx2 = h_ptr[i + n_halves] - 1; // halves_idx(i, 1)
          if (idx1 < m_size && idx2 < m_size) lambda_squ += d_i[idx1] * d_i[idx2];
        }
        lambda_squ /= n_halves;
      }
      bai_var_d_bar = std::max(1e-8, tau_sq - (lambda_squ + delta_sq) / 2.0) / m_size;
    }

    if (convex_flag && nRT > 1 && nRC > 1 && m_size > 0 && ssqR > 0) {
      double w_star = ssqR / (ssqR + bai_var_d_bar);
      res_ptr[b] = w_star * d_bar + (1 - w_star) * r_bar;
    } else if (m_size > 0) {
      res_ptr[b] = d_bar;
    } else {
      res_ptr[b] = r_bar;
    }
  }
  return wrap(results_vec);
}
