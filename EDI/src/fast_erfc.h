#ifndef EDI_FAST_ERFC_H
#define EDI_FAST_ERFC_H

#include <cmath>
#include <limits>

// Cephes piecewise rational approximation for erfc. Falls back to libm for
// |x| > 5.6 to preserve extreme-tail behavior outside the probit fit range.
inline double fast_erfc_polevl(double x, const double* coef, int degree) {
    double result = coef[0];
    for (int i = 1; i <= degree; ++i) result = result * x + coef[i];
    return result;
}

inline double fast_erfc_p1evl(double x, const double* coef, int degree) {
    double result = x + coef[0];
    for (int i = 1; i < degree; ++i) result = result * x + coef[i];
    return result;
}

inline double fast_erfc(double x) {
    const double ax = std::abs(x);
    if (ax > 5.6) return std::erfc(x);

    if (ax < 1.0) {
        static constexpr double T[] = {
            9.60497373987051638749E0,
            9.00260197203842689217E1,
            2.23200534594684319226E3,
            7.00332514112805075473E3,
            5.55923013010394962768E4
        };
        static constexpr double U[] = {
            3.35617141647503099647E1,
            5.21357949780152679795E2,
            4.59432382970980127987E3,
            2.26290000613890934246E4,
            4.92673942608635921086E4
        };
        const double z = x * x;
        const double erf_x = x * fast_erfc_polevl(z, T, 4) /
            fast_erfc_p1evl(z, U, 5);
        return 1.0 - erf_x;
    }

    static constexpr double P[] = {
        2.46196981473530512524E-10,
        5.64189564831068821977E-1,
        7.46321056442269912687E0,
        4.86371970985681366614E1,
        1.96520832956077098289E2,
        5.26445194995477358631E2,
        9.34528527171957607540E2,
        1.02755188689515710272E3,
        5.57535335369399327526E2
    };
    static constexpr double Q[] = {
        1.32281951154744992508E1,
        8.67072140885989742329E1,
        3.54937778887819891062E2,
        9.75708501743205489753E2,
        1.82390916687909736289E3,
        2.24633760818710981792E3,
        1.65666309194161341348E3,
        5.57535340817727675546E2
    };
    const double y = std::exp(-ax * ax) *
        fast_erfc_polevl(ax, P, 8) / fast_erfc_p1evl(ax, Q, 8);
    return x < 0.0 ? 2.0 - y : y;
}

static constexpr double kSqrt1_2   = 0.7071067811865476;
static constexpr double k1_Sqrt2Pi = 0.3989422804014327;

// Standard normal CDF via fast_erfc. Clamped at ±8 to prevent log(0).
inline double pnorm_fast(double x) {
    if (x >= 8.0) return 1.0 - 6e-16;
    if (x <= -8.0) return 6e-16;
    return 0.5 * fast_erfc(-x * kSqrt1_2);
}

// Standard normal PDF: phi(x) = exp(-x^2/2) / sqrt(2*pi)
inline double dnorm_fast(double x) {
    return k1_Sqrt2Pi * std::exp(-0.5 * x * x);
}

// Log standard normal CDF: log(Phi(x)); avoids R::pnorm dispatch overhead.
inline double fast_log_pnorm(double x) {
    if (x >= 8.0)  return -6.661e-16;        // log(1 - 6e-16)
    if (x <= -8.0) return -35.0496;           // log(6e-16)
    return std::log(0.5 * fast_erfc(-x * kSqrt1_2));
}

// Log standard normal PDF: log(phi(x)) = -x^2/2 - log(sqrt(2*pi)).
inline double fast_log_dnorm(double x) {
    constexpr double kLnSqrt2Pi = 0.9189385332046728;
    return -0.5 * x * x - kLnSqrt2Pi;
}

// Standard normal quantile function (inverse CDF), i.e. R::qnorm5(p, 0, 1, 1, 0).
// Peter Acklam's rational approximation, accurate to ~1.2e-9 relative error
// across the full (0, 1) range including the extreme tails. A Halley
// refinement step (using pnorm_fast/dnorm_fast) was tried to push this
// closer to machine precision, but cost more than R::qnorm5's own AS
// 241-based implementation while both real call sites (cold-start/
// warm-start heuristics in fast_probit_regression.cpp and
// fast_ordinal_probit_regression.cpp, refined further by the optimizer
// regardless) have no need for better than ~1e-9 — so the refinement was
// dropped in favor of the ~2x faster unrefined approximation.
inline double fast_qnorm(double p) {
    if (!(p > 0.0) || !(p < 1.0)) {
        if (p == 0.0) return -HUGE_VAL;
        if (p == 1.0) return HUGE_VAL;
        return std::numeric_limits<double>::quiet_NaN();
    }

    static constexpr double a[] = {
        -3.969683028665376e+01,  2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00
    };
    static constexpr double b[] = {
        -5.447609879822406e+01,  1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01
    };
    static constexpr double c[] = {
        -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
        -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00
    };
    static constexpr double d[] = {
         7.784695709041462e-03,  3.224671290700398e-01,  2.445134137142996e+00,
         3.754408661907416e+00
    };
    constexpr double p_low = 0.02425;
    constexpr double p_high = 1.0 - p_low;

    double x;
    if (p < p_low) {
        const double q = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    } else if (p <= p_high) {
        const double q = p - 0.5, r = q * q;
        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
    } else {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    return x;
}

#endif // EDI_FAST_ERFC_H
