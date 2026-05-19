#ifndef EDI_GLMM_LINKS_H
#define EDI_GLMM_LINKS_H

#include <cmath>
#include <Rmath.h>
#include "_helper_functions.h"

namespace glmm {

struct LogitLink {
    static inline double cdf(double x) { return plogis_safe(x); }
    static inline double pdf(double x) {
        double p = cdf(x);
        return p * (1.0 - p);
    }
    // pdf when the caller already has F = cdf(x): avoids a redundant cdf() call.
    static inline double pdf_from_cdf(double /*x*/, double F) { return F * (1.0 - F); }
    static inline double deriv_pdf(double x) {
        double p = cdf(x);
        return p * (1.0 - p) * (1.0 - 2.0 * p);
    }
};

struct ProbitLink {
    static inline double cdf(double x) { return R::pnorm(x, 0.0, 1.0, 1, 0); }
    static inline double pdf(double x) { return R::dnorm(x, 0.0, 1.0, 0); }
    static inline double pdf_from_cdf(double x, double /*F*/) { return R::dnorm(x, 0.0, 1.0, 0); }
    static inline double deriv_pdf(double x) { return -x * R::dnorm(x, 0.0, 1.0, 0); }
};

struct CauchitLink {
    static inline double cdf(double x) { return 0.5 + std::atan(x) / M_PI; }
    static inline double pdf(double x) { return 1.0 / (M_PI * (1.0 + x * x)); }
    static inline double pdf_from_cdf(double x, double /*F*/) { return 1.0 / (M_PI * (1.0 + x * x)); }
    static inline double deriv_pdf(double x) { return -2.0 * x / (M_PI * std::pow(1.0 + x * x, 2)); }
};

struct CloglogLink {
    static inline double cdf(double x) {
        if (x > 100.0) return 1.0;
        if (x < -700.0) return 0.0;
        return 1.0 - std::exp(-std::exp(x));
    }
    static inline double pdf(double x) {
        if (x > 700.0 || x < -700.0) return 0.0;
        const double ex = std::exp(x);
        return ex * std::exp(-ex);
    }
    // cloglog pdf does not simplify given F alone, so ignore F.
    static inline double pdf_from_cdf(double x, double /*F*/) {
        if (x > 700.0 || x < -700.0) return 0.0;
        const double ex = std::exp(x);
        return ex * std::exp(-ex);
    }
    static inline double deriv_pdf(double x) {
        if (x > 700.0 || x < -700.0) return 0.0;
        const double ex = std::exp(x);
        return ex * std::exp(-ex) * (1.0 - ex);
    }
};

} // namespace glmm

#endif
