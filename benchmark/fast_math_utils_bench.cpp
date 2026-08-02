// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
// fast_gamma_functions.h/fast_erfc.h/ordinal_fixed_link_helpers.h are leaf
// headers (only depend on RcppEigen/Rmath/cmath/limits) so they compile
// standalone under Rcpp::sourceCpp() exactly like
// benchmark/fast_trigamma_speed_compare.cpp. fast_log1pexp now lives in the
// heavier _helper_functions.h (which pulls in LBFGSpp via
// "optimization/LBFGS.h") — resolved here via Rcpp::depends(RcppNumerical),
// the same precedent fast_trigamma_speed_compare.cpp already established,
// rather than copying the function out.
//
// Only the fast_* kernels NOT yet promoted to real EDI package exports live
// here (fast_digamma_vec_cpp, fast_trigamma_vec_cpp, fast_lgamma_vec_cpp,
// fast_lbeta_vec_cpp, fast_dnbinom_mu_vec_cpp, fast_qnorm_vec_cpp,
// fast_log_pnorm_vec_cpp, fast_log_dnorm_vec_cpp are now real EDI::
// exports — see EDI/src/fast_math_utils.cpp — and are called directly from
// the installed package in benchmark_model_fits.R, not sourced from here;
// duplicating them in this file would shadow the real exports once both are
// loaded into the same R session). These six are still prototype/benchmark-
// only kernels: promote them to EDI/src/fast_math_utils.cpp the same way,
// following the same "consistently faster in the R benchmark" precedent,
// once this file's own Utility table run confirms it.

#include "../EDI/src/fast_gamma_functions.h"
#include "../EDI/src/fast_erfc.h"
#include "../EDI/src/ordinal_fixed_link_helpers.h"
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
Rcpp::NumericVector fast_pchisq_upper_vec_cpp(Rcpp::NumericVector statistic, Rcpp::NumericVector df) {
    const int n = statistic.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_pchisq_upper(statistic[i], df[i]);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_erfc_vec_cpp(Rcpp::NumericVector x) {
    const int n = x.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_erfc(x[i]);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector pnorm_fast_vec_cpp(Rcpp::NumericVector x) {
    const int n = x.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = pnorm_fast(x[i]);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector dnorm_fast_vec_cpp(Rcpp::NumericVector x) {
    const int n = x.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = dnorm_fast(x[i]);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_atan_vec_cpp(Rcpp::NumericVector x) {
    const int n = x.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = edi_ordinal::fast_atan(x[i]);
    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_log1pexp_vec_cpp(Rcpp::NumericVector x) {
    const int n = x.size();
    Rcpp::NumericVector out(n);
    for (int i = 0; i < n; ++i) out[i] = fast_log1pexp(x[i]);
    return out;
}
