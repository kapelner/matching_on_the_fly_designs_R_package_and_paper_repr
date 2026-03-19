#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <cmath>
#include <vector>

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

inline double logit_cpp(double x) {
    return std::log(x / (1.0 - x));
}

inline double inv_logit_cpp(double x) {
    if (x >= 0.0) {
        const double z = std::exp(-x);
        return 1.0 / (1.0 + z);
    }
    const double z = std::exp(x);
    return z / (1.0 + z);
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

// Signed-rank HL estimate: median of Walsh averages of differences
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

// Simplified variance estimate for HL based on asymptotic normality of Wilcoxon stat
// This is a proxy for the SE back-calculated from CI width.
// For large n, SE(HL) is proportional to 1/sqrt(n)
// We'll use a robust sample-based approach if possible, but for IVWC weights, 
// even a consistent relative variance is helpful.
// Better: implement the actual Wilcoxon SE formula for the statistic and scale it.
double estimate_hl_ssq_rank_sum(const std::vector<double>& y_t, const std::vector<double>& y_c) {
    if (y_t.size() < 2 || y_c.size() < 2) return NA_REAL;
    // Use a simple heuristic: variance of pairwise differences / effective n
    // This is often a good proxy for the HL variance.
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

double apply_shift(double y_val, double delta, int transform_code) {
    if (transform_code == 1) {
        return y_val * std::exp(delta);
    }
    if (transform_code == 2) {
        return inv_logit_cpp(logit_cpp(y_val) + delta);
    }
    if (transform_code == 3) {
        return (y_val + 1.0) * std::exp(delta) - 1.0;
    }
    return y_val + delta;
}

} // namespace

// [[Rcpp::export]]
double wilcox_hl_point_estimate_cpp(const NumericVector& y, const IntegerVector& w) {
    if (y.size() != w.size()) {
        stop("y and w must have the same length");
    }

    std::vector<double> y_t;
    std::vector<double> y_c;
    y_t.reserve(y.size());
    y_c.reserve(y.size());

    for (int i = 0; i < y.size(); ++i) {
        if (!R_finite(y[i])) {
            continue;
        }
        if (w[i] == 1) {
            y_t.push_back(y[i]);
        } else if (w[i] == 0) {
            y_c.push_back(y[i]);
        }
    }

    if (y_t.empty() || y_c.empty()) {
        return NA_REAL;
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
    if (w.size() != n || indices_mat.nrow() != n) {
        stop("dimension mismatch in compute_wilcox_hl_bootstrap_parallel_cpp");
    }

    NumericVector results(B, NA_REAL);

#ifdef _OPENMP
    if (num_cores > 1) {
        omp_set_num_threads(num_cores);
    }
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
        if (indices_mat(0, b) < 0) {
            results[b] = NA_REAL;
            continue;
        }

        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            const int idx = indices_mat(i, b);
            if (idx < 0 || idx >= n) {
                y_t.clear();
                y_c.clear();
                break;
            }
            const double y_val = y[idx];
            if (!R_finite(y_val)) {
                continue;
            }
            if (w[idx] == 1) {
                y_t.push_back(y_val);
            } else if (w[idx] == 0) {
                y_c.push_back(y_val);
            }
        }

        results[b] = hl_from_groups(y_t, y_c);
    }

    return results;
}

// [[Rcpp::export]]
NumericVector compute_wilcox_hl_distr_parallel_cpp(
    const NumericVector& y,
    const IntegerMatrix& w_mat,
    double delta,
    int transform_code,
    int num_cores) {

    const int n = y.size();
    const int nsim = w_mat.ncol();
    if (w_mat.nrow() != n) {
        stop("dimension mismatch in compute_wilcox_hl_distr_parallel_cpp");
    }

    NumericVector results(nsim, NA_REAL);

    // Pre-calculate shifted responses for treatment group once
    std::vector<double> y_shifted(n);
    if (delta != 0.0) {
        for (int i = 0; i < n; ++i) {
            if (R_finite(y[i])) {
                y_shifted[i] = apply_shift(y[i], delta, transform_code);
            } else {
                y_shifted[i] = NA_REAL;
            }
        }
    }

#ifdef _OPENMP
    if (num_cores > 1) {
        omp_set_num_threads(num_cores);
    }
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        std::vector<double> y_t;
        std::vector<double> y_c;
        y_t.reserve(n);
        y_c.reserve(n);

        for (int i = 0; i < n; ++i) {
            const double y_orig = y[i];
            if (!R_finite(y_orig)) {
                continue;
            }

            const int w_val = w_mat(i, b);
            if (w_val == 1) {
                y_t.push_back(delta != 0.0 ? y_shifted[i] : y_orig);
            } else if (w_val == 0) {
                y_c.push_back(y_orig);
            }
        }

        results[b] = hl_from_groups(y_t, y_c);
    }

    return results;
}

// [[Rcpp::export]]
NumericVector compute_wilcox_kk_ivwc_bootstrap_parallel_cpp(
    const NumericVector& y,
    const IntegerVector& w,
    const IntegerVector& match_indic,
    const IntegerMatrix& indices_mat,
    const IntegerMatrix& match_indic_mat,
    int num_cores) {

    int B = indices_mat.ncol();
    int n = y.size();
    NumericVector results(B, NA_REAL);

#ifdef _OPENMP
    if (num_cores > 1) {
        omp_set_num_threads(num_cores);
    }
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
        std::vector<double> pair_diffs;
        std::vector<double> res_y_t;
        std::vector<double> res_y_c;

        // Process resampled data for this bootstrap iteration
        // We need to group by the new match_indic_mat
        std::map<int, std::vector<size_t>> pairs;
        for (int i = 0; i < n; ++i) {
            int mid = match_indic_mat(i, b);
            int idx = indices_mat(i, b) - 1;
            if (mid == 0) {
                if (w[idx] == 1) res_y_t.push_back(y[idx]);
                else res_y_c.push_back(y[idx]);
            } else {
                pairs[mid].push_back(idx);
            }
        }

        for (auto const& [mid, idxs] : pairs) {
            if (idxs.size() == 2) {
                int idx1 = idxs[0];
                int idx2 = idxs[1];
                if (w[idx1] == 1 && w[idx2] == 0) {
                    pair_diffs.push_back(y[idx1] - y[idx2]);
                } else if (w[idx1] == 0 && w[idx2] == 1) {
                    pair_diffs.push_back(y[idx2] - y[idx1]);
                }
            }
        }

        // Compute component estimates
        double beta_m = hl_signed_rank(pair_diffs);
        double ssq_m = estimate_hl_ssq_signed_rank(pair_diffs);
        
        double beta_r = hl_from_groups(res_y_t, res_y_c);
        double ssq_r = estimate_hl_ssq_rank_sum(res_y_t, res_y_c);

        bool m_ok = R_finite(beta_m) && R_finite(ssq_m) && ssq_m > 0;
        bool r_ok = R_finite(beta_r) && R_finite(ssq_r) && ssq_r > 0;

        if (m_ok && r_ok) {
            double w_star = ssq_r / (ssq_r + ssq_m);
            results[b] = w_star * beta_m + (1.0 - w_star) * beta_r;
        } else if (m_ok) {
            results[b] = beta_m;
        } else if (r_ok) {
            results[b] = beta_r;
        } else {
            results[b] = NA_REAL;
        }
    }

    return results;
}
