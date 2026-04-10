#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

namespace {

double median_in_place(std::vector<double>& values) {
    const size_t n = values.size();
    if (n == 0) {
        return NA_REAL;
    }

    const size_t mid = n / 2;
    std::nth_element(values.begin(), values.begin() + mid, values.end());
    const double upper = values[mid];

    if (n % 2 == 1) {
        return upper;
    }

    std::nth_element(values.begin(), values.begin() + mid - 1, values.begin() + mid);
    const double lower = values[mid - 1];
    return 0.5 * (lower + upper);
}

inline double logit_cpp(double x, double clamp) {
    if (x < clamp) x = clamp;
    if (x > 1.0 - clamp) x = 1.0 - clamp;
    return std::log(x / (1.0 - x));
}

inline double inv_logit_cpp(double x, double clamp) {
    double p;
    if (x >= 0.0) {
        const double z = std::exp(-x);
        p = 1.0 / (1.0 + z);
    } else {
        const double z = std::exp(x);
        p = z / (1.0 + z);
    }
    if (p < clamp) p = clamp;
    if (p > 1.0 - clamp) p = 1.0 - clamp;
    return p;
}

double hl_from_groups(const std::vector<double>& y_t, const std::vector<double>& y_c) {
    if (y_t.empty() || y_c.empty()) {
        return NA_REAL;
    }

    std::vector<double> diffs;
    diffs.reserve(y_t.size() * y_c.size());
    for (double yt : y_t) {
        for (double yc : y_c) {
            diffs.push_back(yt - yc);
        }
    }

    return median_in_place(diffs);
}

double hl_signed_rank(const std::vector<double>& pair_diffs) {
    if (pair_diffs.empty()) {
        return NA_REAL;
    }

    size_t m = pair_diffs.size();
    std::vector<double> walsh_avgs;
    walsh_avgs.reserve(m * (m + 1) / 2);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = i; j < m; ++j) {
            walsh_avgs.push_back(0.5 * (pair_diffs[i] + pair_diffs[j]));
        }
    }

    return median_in_place(walsh_avgs);
}

double estimate_hl_ssq_rank_sum(const std::vector<double>& y_t, const std::vector<double>& y_c) {
    if (y_t.size() < 2 || y_c.size() < 2) return NA_REAL;
    double sum = 0;
    double sum_sq = 0;
    int count = 0;
    for (double yt : y_t) {
        for (double yc : y_c) {
            double d = yt - yc;
            sum += d;
            sum_sq += d * d;
            count++;
        }
    }
    double var_diffs = (sum_sq - (sum * sum) / count) / (count - 1);
    return var_diffs / (y_t.size() + y_c.size());
}

double estimate_hl_ssq_signed_rank(const std::vector<double>& pair_diffs) {
    if (pair_diffs.size() < 2) return NA_REAL;
    double sum = 0;
    double sum_sq = 0;
    int count = 0;
    for (size_t i = 0; i < pair_diffs.size(); ++i) {
        for (size_t j = i; j < pair_diffs.size(); ++j) {
            double a = 0.5 * (pair_diffs[i] + pair_diffs[j]);
            sum += a;
            sum_sq += a * a;
            count++;
        }
    }
    double var_walsh = (sum_sq - (sum * sum) / count) / (count - 1);
    return var_walsh / pair_diffs.size();
}

double apply_shift(double y_val, double delta, int transform_code, double zero_one_logit_clamp) {
    if (transform_code == 1) {
        return y_val * std::exp(delta);
    }
    if (transform_code == 2) {
        return inv_logit_cpp(logit_cpp(y_val, zero_one_logit_clamp) + delta, zero_one_logit_clamp);
    }
    if (transform_code == 3) {
        return (y_val + 1.0) * std::exp(delta) - 1.0;
    }
    return y_val + delta;
}

} // namespace

// [[Rcpp::export]]
double wilcox_hl_point_estimate_cpp(const NumericVector& y, const IntegerVector& w) {
    std::vector<double> y_t;
    std::vector<double> y_c;
    y_t.reserve(y.size());
    y_c.reserve(y.size());

    const double* y_ptr = y.begin();
    const int* w_ptr = w.begin();
    int n = y.size();

    for (int i = 0; i < n; ++i) {
        if (!std::isfinite(y_ptr[i])) continue;
        if (w_ptr[i] == 1) y_t.push_back(y_ptr[i]);
        else if (w_ptr[i] == 0) y_c.push_back(y_ptr[i]);
    }

    return hl_from_groups(y_t, y_c);
}

// [[Rcpp::export]]
NumericVector compute_wilcox_hl_bootstrap_parallel_cpp(
    const NumericVector& y,
    const IntegerVector& w,
    const IntegerMatrix& indices_mat,
    int num_cores) {

    const int n = y.size();
    const int B = indices_mat.ncol();
    
    std::vector<double> results_vec(B, NA_REAL);
    const double* y_ptr = y.begin();
    const int* w_ptr = w.begin();
    const int* idx_ptr = indices_mat.begin();
    double* res_ptr = results_vec.data();

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
        const int* idx_col = idx_ptr + (size_t)b * n;
        if (idx_col[0] < 0) {
            res_ptr[b] = NA_REAL;
            continue;
        }

        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            const int idx = idx_col[i];
            if (idx < 0 || idx >= n) continue;
            
            const double y_val = y_ptr[idx];
            if (!std::isfinite(y_val)) continue;
            
            if (w_ptr[idx] == 1) y_t.push_back(y_val);
            else if (w_ptr[idx] == 0) y_c.push_back(y_val);
        }

        res_ptr[b] = hl_from_groups(y_t, y_c);
    }

    return wrap(results_vec);
}

// [[Rcpp::export]]
NumericVector compute_wilcox_hl_distr_parallel_cpp(
    const NumericVector& y,
    const IntegerMatrix& w_mat,
    double delta,
    int transform_code,
    double zero_one_logit_clamp,
    int num_cores) {

    const int n = y.size();
    const int nsim = w_mat.ncol();
    
    std::vector<double> results_vec(nsim, NA_REAL);
    const double* y_ptr = y.begin();
    const int* w_ptr = w_mat.begin();
    double* res_ptr = results_vec.data();

    std::vector<double> y_shifted(n);
    if (delta != 0.0) {
        for (int i = 0; i < n; ++i) {
            if (std::isfinite(y_ptr[i])) y_shifted[i] = apply_shift(y_ptr[i], delta, transform_code, zero_one_logit_clamp);
            else y_shifted[i] = NA_REAL;
        }
    }

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        const int* w_col = w_ptr + (size_t)b * n;
        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(y_ptr[i])) continue;

            if (w_col[i] == 1) {
                y_t.push_back(delta != 0.0 ? y_shifted[i] : y_ptr[i]);
            } else if (w_col[i] == 0) {
                y_c.push_back(y_ptr[i]);
            }
        }

        res_ptr[b] = hl_from_groups(y_t, y_c);
    }

    return wrap(results_vec);
}

// [[Rcpp::export]]
NumericVector compute_wilcox_kk_ivwc_bootstrap_parallel_cpp(
    const NumericVector& y,
    const IntegerVector& w,
    const IntegerVector& m_vec,
    const IntegerMatrix& indices_mat,
    const IntegerMatrix& m_mat,
    int num_cores) {

    int B = indices_mat.ncol();
    int n = y.size();
    
    std::vector<double> results_vec(B, NA_REAL);
    const double* y_ptr = y.begin();
    const int* w_ptr = w.begin();
    const int* idx_ptr = indices_mat.begin();
    const int* m_ptr = m_mat.begin();
    double* res_ptr = results_vec.data();

#ifdef _OPENMP
    omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
        const int* idx_col = idx_ptr + (size_t)b * n;
        const int* m_col = m_ptr + (size_t)b * n;
        
        std::vector<double> pair_diffs;
        std::vector<double> res_y_t;
        std::vector<double> res_y_c;

        std::map<int, std::vector<size_t>> pairs;
        for (int i = 0; i < n; ++i) {
            int mid = m_col[i];
            int idx = idx_col[i] - 1;
            if (idx < 0) continue;
            
            if (mid == 0) {
                if (w_ptr[idx] == 1) res_y_t.push_back(y_ptr[idx]);
                else res_y_c.push_back(y_ptr[idx]);
            } else {
                pairs[mid].push_back(idx);
            }
        }

        for (auto const& [mid, idxs] : pairs) {
            if (idxs.size() == 2) {
                int idx1 = idxs[0];
                int idx2 = idxs[1];
                if (w_ptr[idx1] == 1 && w_ptr[idx2] == 0) pair_diffs.push_back(y_ptr[idx1] - y_ptr[idx2]);
                else if (w_ptr[idx1] == 0 && w_ptr[idx2] == 1) pair_diffs.push_back(y_ptr[idx2] - y_ptr[idx1]);
            }
        }

        double beta_m = hl_signed_rank(pair_diffs);
        double ssq_m = estimate_hl_ssq_signed_rank(pair_diffs);
        double beta_r = hl_from_groups(res_y_t, res_y_c);
        double ssq_r = estimate_hl_ssq_rank_sum(res_y_t, res_y_c);

        bool m_ok = std::isfinite(beta_m) && std::isfinite(ssq_m) && ssq_m > 0;
        bool r_ok = std::isfinite(beta_r) && std::isfinite(ssq_r) && ssq_r > 0;

        if (m_ok && r_ok) {
            double w_star = ssq_r / (ssq_r + ssq_m);
            res_ptr[b] = w_star * beta_m + (1.0 - w_star) * beta_r;
        } else if (m_ok) {
            res_ptr[b] = beta_m;
        } else if (r_ok) {
            res_ptr[b] = beta_r;
        }
    }

    return wrap(results_vec);
}
