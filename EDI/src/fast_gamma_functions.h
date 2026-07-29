#ifndef EDI_FAST_GAMMA_FUNCTIONS_H
#define EDI_FAST_GAMMA_FUNCTIONS_H

// Leaf header (no dependency on _helper_functions.h/optimization_starts.h) so
// both can use these without a circular include: _helper_functions.h already
// includes optimization_starts.h, so if optimization_starts.h needed these
// functions via _helper_functions.h, that would cycle.

#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <limits>

// Fast digamma via A&S 6.3.18 asymptotic expansion + recurrence shift.
// Accurate to ≤ 4e-12 relative error for x > 0; falls back to R::digamma for x <= 0.
inline double fast_digamma(double x) {
    if (x <= 0.0) return R::digamma(x);
    double r = 0.0;
    while (x < 8.0) { r -= 1.0 / x; x += 1.0; }
    const double ix = 1.0 / x, ix2 = ix * ix;
    r += std::log(x) - 0.5 * ix
       - ix2 * (1.0/12.0 - ix2 * (1.0/120.0 - ix2 * (1.0/252.0 - ix2 * (1.0/240.0 - ix2/132.0))));
    return r;
}

// Fast trigamma (psi') via A&S 6.4.12 asymptotic expansion + recurrence shift.
// Accurate to ≤ 3e-12 relative error for x > 0; falls back to R::trigamma for x <= 0.
inline double fast_trigamma(double x) {
    if (x <= 0.0) return R::trigamma(x);
    double r = 0.0;
    while (x < 8.0) { r += 1.0 / (x * x); x += 1.0; }
    const double ix = 1.0 / x, ix2 = ix * ix;
    r += ix + 0.5 * ix2
       + ix2 * ix * (1.0/6.0 - ix2 * (1.0/30.0 - ix2 * (1.0/42.0 - ix2 * (1.0/30.0 - ix2 * (5.0/66.0)))));
    return r;
}

// Vectorized fast_trigamma over an Eigen array. Requires x[i] > 0 for all i
// (callers with y+theta / theta arguments already guarantee this) — use the
// scalar fast_trigamma for values that may be <= 0. The recurrence shift is
// applied with a per-lane mask instead of a per-lane branch, and the shift
// count is bounded by the smallest element, so this stays a fixed small
// number of full-width Eigen array ops that the compiler can vectorize under
// -O3 -march=native, rather than n independent scalar loops.
inline Eigen::ArrayXd fast_trigamma_vec(const Eigen::ArrayXd& x) {
    const double min_x = x.size() > 0 ? x.minCoeff() : 8.0;
    const int n_shifts = min_x < 8.0 ? static_cast<int>(std::ceil(8.0 - min_x)) : 0;

    Eigen::ArrayXd xx = x;
    Eigen::ArrayXd r = Eigen::ArrayXd::Zero(x.size());
    // Reused across shift iterations instead of reallocated each pass — perf
    // showed malloc/free churn dominating this loop when active_d was a
    // fresh per-iteration temporary (see perf_experiments_final.md TODO-129).
    Eigen::ArrayXd active_d(x.size());
    for (int k = 0; k < n_shifts; ++k) {
        active_d = (xx < 8.0).cast<double>();
        r += active_d / (xx * xx);
        xx += active_d;
    }

    const Eigen::ArrayXd ix = xx.inverse();
    const Eigen::ArrayXd ix2 = ix * ix;
    r += ix + 0.5 * ix2
       + ix2 * ix * (1.0/6.0 - ix2 * (1.0/30.0 - ix2 * (1.0/42.0 - ix2 * (1.0/30.0 - ix2 * (5.0/66.0)))));
    return r;
}

inline double fast_lgamma_stirling(double x) {
    const double inv = 1.0 / x;
    const double inv2 = inv * inv;
    const double series = inv * (
        1.0 / 12.0 + inv2 * (
        -1.0 / 360.0 + inv2 * (
         1.0 / 1260.0 + inv2 * (
        -1.0 / 1680.0 + inv2 * (
         1.0 / 1188.0 + inv2 * (
        -691.0 / 360360.0 + inv2 * (1.0 / 156.0)))))));

    return (x - 0.5) * std::log(x) - x
        + 0.91893853320467274178032973640562 + series;
}

inline double fast_lgamma_lanczos(double x) {
    const double z = x - 1.0;
    double a = 0.99999999999980993;
    a += 676.5203681218851 / (z + 1.0);
    a += -1259.1392167224028 / (z + 2.0);
    a += 771.32342877765313 / (z + 3.0);
    a += -176.61502916214059 / (z + 4.0);
    a += 12.507343278686905 / (z + 5.0);
    a += -0.13857109526572012 / (z + 6.0);
    a += 9.9843695780195716e-6 / (z + 7.0);
    a += 1.5056327351493116e-7 / (z + 8.0);
    const double t = z + 7.5;
    return 0.91893853320467274178032973640562 + (z + 0.5) * std::log(t) - t + std::log(a);
}

// Fast log-gamma for positive arguments. Uses a Lanczos rational approximation
// for moderate x and Stirling for large x; falls back for nonpositive/nonfinite.
inline double fast_lgamma(double x) {
    if (x <= 0.0 || !std::isfinite(x)) return std::lgamma(x);
    if (x < 0.5) return fast_lgamma_lanczos(x + 1.0) - std::log(x);
    if (x < 8.0) return fast_lgamma_lanczos(x);
    return fast_lgamma_stirling(x);
}

// Fast log-beta: lbeta(a, b) = lgamma(a) + lgamma(b) - lgamma(a + b), via fast_lgamma.
inline double fast_lbeta(double a, double b) {
    return fast_lgamma(a) + fast_lgamma(b) - fast_lgamma(a + b);
}

// Fast negative-binomial PMF/log-PMF parameterized by mean mu (matches
// R::dnbinom_mu(x, size, mu, give_log)), via fast_lgamma in place of the
// three R::lgamma dispatches R's own dnbinom_mu implementation performs.
// Parameter/local names avoid "give_log"/"log_p": Rcpp's dpq macros.h
// #defines give_log to log_p, which would otherwise silently rewrite this
// signature and collide with a same-named local.
inline double fast_dnbinom_mu(double x, double size, double mu, bool return_log) {
    if (x < 0.0 || x != std::floor(x) || !std::isfinite(mu) || mu < 0.0 || size < 0.0) {
        return return_log ? -std::numeric_limits<double>::infinity() : 0.0;
    }
    if (mu == 0.0) {
        const double lp = (x == 0.0) ? 0.0 : -std::numeric_limits<double>::infinity();
        return return_log ? lp : std::exp(lp);
    }
    const double logp_val = fast_lgamma(x + size) - fast_lgamma(size) - fast_lgamma(x + 1.0)
        + size * (std::log(size) - std::log(size + mu))
        + x * (std::log(mu) - std::log(size + mu));
    return return_log ? logp_val : std::exp(logp_val);
}

// Regularized lower incomplete gamma P(a,x) via series expansion (A&S
// 6.5.29). Valid/convergent for x < a+1; used as a building block for
// fast_gammaq below, not called directly at any real call site.
inline double fast_gammap_series(double a, double x) {
    if (x <= 0.0) return 0.0;
    double ap = a;
    double sum = 1.0 / a;
    double del = sum;
    for (int n = 1; n <= 500; ++n) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (std::abs(del) < std::abs(sum) * 1e-16) break;
    }
    return sum * std::exp(-x + a * std::log(x) - fast_lgamma(a));
}

// Regularized upper incomplete gamma Q(a,x) via a continued fraction
// (modified Lentz's method, A&S 6.5.31). Valid/convergent for x >= a+1.
inline double fast_gammaq_cf(double a, double x) {
    constexpr double kMin = 1e-300;
    double b = x + 1.0 - a;
    double c = 1.0 / kMin;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i <= 500; ++i) {
        const double an = -static_cast<double>(i) * (static_cast<double>(i) - a);
        b += 2.0;
        d = an * d + b;
        if (std::abs(d) < kMin) d = kMin;
        c = b + an / c;
        if (std::abs(c) < kMin) c = kMin;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < 1e-16) break;
    }
    return std::exp(-x + a * std::log(x) - fast_lgamma(a)) * h;
}

// Regularized upper incomplete gamma Q(a,x) = 1 - P(a,x). Dispatches between
// the series (x < a+1) and continued-fraction (x >= a+1) forms, whichever
// converges quickly at that (a,x); both compute exp(-x + a*log(x) - lgamma(a))
// as their leading factor, so accuracy is governed by fast_lgamma's already
// -validated precision plus the truncated-series/CF residual (bounded by the
// 1e-16 stopping tolerance in each loop).
inline double fast_gammaq(double a, double x) {
    if (!(a > 0.0) || x < 0.0) return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return 1.0;
    return (x < a + 1.0) ? (1.0 - fast_gammap_series(a, x)) : fast_gammaq_cf(a, x);
}

// Upper-tail chi-squared p-value: P(X > statistic) for X ~ chi-squared(df).
// Matches R::pchisq(statistic, df, /*lower_tail=*/false, /*log_p=*/false),
// which is the only (lower_tail, log_p) combination any call site in this
// codebase uses -- via the standard chi-squared/gamma relation
// P(X > x) = Q(df/2, x/2).
inline double fast_pchisq_upper(double statistic, double df) {
    if (!std::isfinite(statistic) || !std::isfinite(df) || df <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (statistic <= 0.0) return 1.0;
    return fast_gammaq(df / 2.0, statistic / 2.0);
}

#endif
