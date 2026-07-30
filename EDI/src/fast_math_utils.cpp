// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "fast_gamma_functions.h"
#include "fast_erfc.h"

using namespace Rcpp;

// Vectorized wrappers around EDI's internal scalar fast_* special-function
// kernels (used internally inside the NegBin/Beta/ZINB/Hurdle likelihoods,
// KK21 negative-binomial fitting, and probit cold-start heuristics). Exported
// as standalone utilities because benchmark/benchmark_model_fits.R's "Utility
// / Math Kernel Performance" table shows every one of them is consistently
// faster than its base R equivalent (see package_metadata/benchmark_model_fits.md).

//' @title Fast Digamma Function (C++)
//' @description Vectorized digamma (psi) function via an asymptotic expansion with recurrence shift, faster than base R's \code{digamma()}.
//' @param x Numeric vector of arguments.
//' @return A numeric vector of digamma values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_digamma_vec_cpp(NumericVector x) {
    const int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_digamma(x[i]);
    return out;
}

//' @title Fast Trigamma Function (C++)
//' @description Vectorized trigamma (psi') function via an asymptotic expansion with recurrence shift, faster than base R's \code{trigamma()}.
//' @param x Numeric vector of arguments.
//' @return A numeric vector of trigamma values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_trigamma_vec_cpp(NumericVector x) {
    const int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_trigamma(x[i]);
    return out;
}

//' @title Fast Log-Gamma Function (C++)
//' @description Vectorized log-gamma function via a Lanczos/Stirling approximation, faster than base R's \code{lgamma()}.
//' @param x Numeric vector of arguments.
//' @return A numeric vector of log-gamma values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_lgamma_vec_cpp(NumericVector x) {
    const int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_lgamma(x[i]);
    return out;
}

//' @title Fast Log-Beta Function (C++)
//' @description Vectorized log-beta function, \code{lgamma(a) + lgamma(b) - lgamma(a + b)}, computed via \code{fast_lgamma}. Faster than base R's \code{lbeta()}.
//' @param a Numeric vector of first shape arguments.
//' @param b Numeric vector of second shape arguments (recycled against \code{a} elementwise; must be the same length).
//' @return A numeric vector of log-beta values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_lbeta_vec_cpp(NumericVector a, NumericVector b) {
    const int n = a.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_lbeta(a[i], b[i]);
    return out;
}

//' @title Fast Mu-Parameterized Negative-Binomial Density (C++)
//' @description Vectorized negative-binomial density parameterized by mean \code{mu}, matching \code{R::dnbinom_mu(x, size, mu, give_log)} semantics via \code{fast_lgamma} in place of R's three \code{lgamma} dispatches. Faster than base R's \code{stats::dnbinom(x, size, mu = mu, log = ...)}.
//' @param x Numeric vector of non-negative integer counts.
//' @param size Dispersion parameter (single value).
//' @param mu Mean parameter (single value).
//' @param return_log Logical. If \code{TRUE}, return the log-density.
//' @return A numeric vector of (log-)density values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_dnbinom_mu_vec_cpp(NumericVector x, double size, double mu, bool return_log) {
    const int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_dnbinom_mu(x[i], size, mu, return_log);
    return out;
}

//' @title Fast Standard Normal Quantile Function (C++)
//' @description Vectorized standard normal inverse CDF via Peter Acklam's rational approximation (accurate to ~1.2e-9 relative error). Faster than base R's \code{stats::qnorm()}.
//' @param p Numeric vector of probabilities in \verb{(0, 1)}.
//' @return A numeric vector of quantiles.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_qnorm_vec_cpp(NumericVector p) {
    const int n = p.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_qnorm(p[i]);
    return out;
}

//' @title Fast Log Standard Normal CDF (C++)
//' @description Vectorized \code{log(Phi(x))} via \code{fast_erfc}, avoiding R's \code{pnorm} dispatch overhead. Faster than base R's \code{stats::pnorm(x, log.p = TRUE)}.
//' @param x Numeric vector of arguments.
//' @return A numeric vector of log CDF values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_log_pnorm_vec_cpp(NumericVector x) {
    const int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_log_pnorm(x[i]);
    return out;
}

//' @title Fast Log Standard Normal PDF (C++)
//' @description Vectorized \code{log(phi(x))} in closed form. Faster than base R's \code{stats::dnorm(x, log = TRUE)}.
//' @param x Numeric vector of arguments.
//' @return A numeric vector of log density values.
//' @export
//' @keywords internal
// [[Rcpp::export]]
NumericVector fast_log_dnorm_vec_cpp(NumericVector x) {
    const int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_log_dnorm(x[i]);
    return out;
}
