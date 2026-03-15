#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

//' Constrained MLE for Risk Difference (Miettinen-Nurminen)
//'
//' Solves the likelihood equations for p_C subject to p_T - p_C = delta.
//' Uses bisection on the score function (derivative of log-likelihood)
//' which is monotonic and well-behaved.
//' @keywords internal
// [[Rcpp::export]]
double mn_constrained_mle_pc_cpp(double x_t, double n_t, double x_c, double n_c, double delta) {
    double lower = std::max(0.0, -delta);
    double upper = std::min(1.0, 1.0 - delta);
    
    // Boundary guards
    if (upper <= lower) return (lower + upper) / 2.0;
    
    double eps = 1e-10;
    double l = lower + eps;
    double u = upper - eps;
    
    auto score_fn = [&](double pc) {
        double pt = pc + delta;
        double s = 0;
        if (pt > 0 && pt < 1) {
            s += x_t / pt - (n_t - x_t) / (1.0 - pt);
        } else if (pt <= 0) {
            s += 1e15; // Large positive
        } else {
            s -= 1e15; // Large negative
        }
        
        if (pc > 0 && pc < 1) {
            s += x_c / pc - (n_c - x_c) / (1.0 - pc);
        } else if (pc <= 0) {
            s += 1e15;
        } else {
            s -= 1e15;
        }
        return s;
    };
    
    // Bisection
    for (int i = 0; i < 60; ++i) {
        double mid = (l + u) / 2.0;
        double val = score_fn(mid);
        if (std::abs(val) < 1e-9) return mid;
        if (val > 0) {
            l = mid;
        } else {
            u = mid;
        }
    }
    return (l + u) / 2.0;
}

//' @keywords internal
// [[Rcpp::export]]
double mn_z_statistic_cpp(double x_t, double n_t, double x_c, double n_c, double delta, double p_t_obs, double p_c_obs) {
    if (n_t == 0 || n_c == 0) return NA_REAL;
    if (delta <= -1.0 || delta >= 1.0) return NA_REAL;
    
    double pc_tilde = mn_constrained_mle_pc_cpp(x_t, n_t, x_c, n_c, delta);
    double pt_tilde = pc_tilde + delta;
    
    double n_tot = n_t + n_c;
    double correction = n_tot / (n_tot - 1.0);
    
    double var_tilde = correction * (pt_tilde * (1.0 - pt_tilde) / n_t + pc_tilde * (1.0 - pc_tilde) / n_c);
    
    if (var_tilde <= 0) return NA_REAL;
    
    return ((p_t_obs - p_c_obs) - delta) / std::sqrt(var_tilde);
}

//' @keywords internal
// [[Rcpp::export]]
double mn_pvalue_cpp(double x_t, double n_t, double x_c, double n_c, double delta, double p_t_obs, double p_c_obs) {
    double z = mn_z_statistic_cpp(x_t, n_t, x_c, n_c, delta, p_t_obs, p_c_obs);
    if (!std::isfinite(z)) return NA_REAL;
    return 2.0 * R::pnorm(-std::abs(z), 0.0, 1.0, 1, 0);
}

//' @keywords internal
// [[Rcpp::export]]
NumericVector mn_ci_cpp(double x_t, double n_t, double x_c, double n_c, double p_t_obs, double p_c_obs, double alpha, double pval_epsilon) {
    double est = p_t_obs - p_c_obs;
    
    auto pval_fn = [&](double delta) {
        return mn_pvalue_cpp(x_t, n_t, x_c, n_c, delta, p_t_obs, p_c_obs);
    };
    
    auto find_bound = [&](double low_search, double high_search) {
        double l = low_search;
        double u = high_search;
        for (int i = 0; i < 50; ++i) {
            double mid = (l + u) / 2.0;
            double p = pval_fn(mid);
            if (!std::isfinite(p)) p = 0.0;
            if (std::abs(p - alpha) < pval_epsilon) return mid;
            if (p > alpha) {
                // Inside CI
                if (mid < est) u = mid; else l = mid;
            } else {
                // Outside CI
                if (mid < est) l = mid; else u = mid;
            }
        }
        return (l + u) / 2.0;
    };
    
    double lower = find_bound(-0.999999, est);
    double upper = find_bound(est, 0.999999);
    
    return NumericVector::create(lower, upper);
}
