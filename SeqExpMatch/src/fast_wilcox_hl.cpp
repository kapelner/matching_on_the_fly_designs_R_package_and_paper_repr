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
            double y_val = y[i];
            if (!R_finite(y_val)) {
                continue;
            }

            const int w_val = w_mat(i, b);
            if (w_val == 1) {
                if (delta != 0.0) {
                    y_val = apply_shift(y_val, delta, transform_code);
                }
                y_t.push_back(y_val);
            } else if (w_val == 0) {
                y_c.push_back(y_val);
            }
        }

        results[b] = hl_from_groups(y_t, y_c);
    }

    return results;
}
