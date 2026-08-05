#ifdef EDI_CORE_ONLY
#include <limits>
constexpr double NA_REAL = std::numeric_limits<double>::quiet_NaN();
#else
#include <Rcpp.h>
using namespace Rcpp;
#endif
#include <cmath>
#include <algorithm>
#include "fast_erfc.h"

struct WilsonCIBounds {
    double lower;
    double upper;
};

// Portable (EDI_CORE_ONLY-safe) core of wilson_score_interval_cpp below --
// identical logic, plain struct instead of Rcpp::NumericVector, so a
// separate Python binding translation unit (and newcombe_independent_ci_
// internal below, in the same TU) can call it directly.
WilsonCIBounds wilson_score_interval_internal(double x, double n, double alpha) {
    if (n <= 0) return WilsonCIBounds{NA_REAL, NA_REAL};
    double p = x / n;
    double z = fast_qnorm(1.0 - alpha / 2.0);
    double z2 = z * z;

    double denom = 2.0 * (n + z2);
    double term1 = 2.0 * n * p + z2;
    double term2 = z * std::sqrt(z2 + 4.0 * n * p * (1.0 - p));

    double lower = (term1 - term2) / denom;
    double upper = (term1 + term2) / denom;

    return WilsonCIBounds{std::max(0.0, lower), std::min(1.0, upper)};
}

//' Wilson Score Interval for a Single Proportion
//' @keywords internal
#ifndef EDI_CORE_ONLY
// [[Rcpp::export]]
NumericVector wilson_score_interval_cpp(double x, double n, double alpha) {
    WilsonCIBounds r = wilson_score_interval_internal(x, n, alpha);
    return NumericVector::create(r.lower, r.upper);
}
#endif

struct NewcombeCIBounds {
    double lower;
    double upper;
};

// Portable (EDI_CORE_ONLY-safe) core of newcombe_independent_ci_cpp below.
NewcombeCIBounds newcombe_independent_ci_internal(double x1, double n1, double x2, double n2, double alpha) {
    if (n1 <= 0 || n2 <= 0) return NewcombeCIBounds{NA_REAL, NA_REAL};

    double p1 = x1 / n1;
    double p2 = x2 / n2;
    double diff = p1 - p2;

    WilsonCIBounds ci1 = wilson_score_interval_internal(x1, n1, alpha);
    WilsonCIBounds ci2 = wilson_score_interval_internal(x2, n2, alpha);

    double l1 = ci1.lower, u1 = ci1.upper;
    double l2 = ci2.lower, u2 = ci2.upper;

    double lower = diff - std::sqrt(std::pow(p1 - l1, 2) + std::pow(u2 - p2, 2));
    double upper = diff + std::sqrt(std::pow(u1 - p1, 2) + std::pow(p2 - l2, 2));

    return NewcombeCIBounds{std::max(-1.0, lower), std::min(1.0, upper)};
}

//' Newcombe Hybrid Score Interval for Independent Proportions (Method 10)
//' @keywords internal
#ifndef EDI_CORE_ONLY
// [[Rcpp::export]]
NumericVector newcombe_independent_ci_cpp(double x1, double n1, double x2, double n2, double alpha) {
    NewcombeCIBounds r = newcombe_independent_ci_internal(x1, n1, x2, n2, alpha);
    return NumericVector::create(r.lower, r.upper);
}
#endif

#ifndef EDI_CORE_ONLY
//' Newcombe Hybrid Score Interval for Paired Proportions
//' @keywords internal
// [[Rcpp::export]]
NumericVector newcombe_paired_ci_cpp(double n11, double n10, double n01, double n00, double alpha) {
    double n = n11 + n10 + n01 + n00;
    if (n <= 0) return NumericVector::create(NA_REAL, NA_REAL);

    double p1 = (n11 + n10) / n;
    double p2 = (n11 + n01) / n;
    double diff = p1 - p2;

    NumericVector ci1 = wilson_score_interval_cpp(n11 + n10, n, alpha);
    NumericVector ci2 = wilson_score_interval_cpp(n11 + n01, n, alpha);

    double l1 = ci1[0], u1 = ci1[1];
    double l2 = ci2[0], u2 = ci2[1];

    // Pearson correlation for 2x2 paired table
    double denom_phi = std::sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00));
    double phi = (denom_phi > 0) ? (n11 * n00 - n10 * n01) / denom_phi : 0.0;

    // Newcombe's paired formula (Method 10 for paired data)
    double d_l = std::pow(p1 - l1, 2) - 2.0 * phi * (p1 - l1) * (u2 - p2) + std::pow(u2 - p2, 2);
    double d_u = std::pow(u1 - p1, 2) - 2.0 * phi * (u1 - p1) * (p2 - l2) + std::pow(p2 - l2, 2);

    double lower = diff - std::sqrt(std::max(0.0, d_l));
    double upper = diff + std::sqrt(std::max(0.0, d_u));

    return NumericVector::create(std::max(-1.0, lower), std::min(1.0, upper));
}
#endif
