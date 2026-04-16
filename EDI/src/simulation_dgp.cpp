#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
List apply_treatment_and_noise_cpp(
    const NumericVector& y_base,
    const IntegerVector& w,
    const std::string&   response_type,
    double               betaT,
    double               sd_noise,
    double               prob_censoring,
    int                  n_ordinal_levels)
{
  const int n = y_base.size();
  if (w.size() != n)
    stop("apply_treatment_and_noise_cpp: y_base and w must have the same length.");

  NumericVector y(n);
  IntegerVector dead(n, 1);   // default: not censored

  // Pre-generate all Gaussian noise in one RNG call
  NumericVector eps = rnorm(n, 0.0, sd_noise);

  if (response_type == "continuous") {
    for (int i = 0; i < n; ++i) {
      const double bt = (w[i] == 1) ? betaT : 0.0;
      y[i] = y_base[i] + bt + eps[i];
    }

  } else if (response_type == "incidence") {
    // Pre-generate all uniform draws for Bernoulli in one call
    NumericVector u = runif(n);
    for (int i = 0; i < n; ++i) {
      const double bt  = (w[i] == 1) ? betaT : 0.0;
      double p_b = y_base[i];
      // clamp to (0.05, 0.95) — y_base is already plogis(scaled), stays in (0,1)
      if (p_b < 0.05) p_b = 0.05;
      if (p_b > 0.95) p_b = 0.95;
      // logit-shift: plogis(qlogis(p_b) + bt + eps)
      const double logit_p = std::log(p_b / (1.0 - p_b)) + bt + eps[i];
      const double p_new   = 1.0 / (1.0 + std::exp(-logit_p));
      y[i] = (u[i] < p_new) ? 1.0 : 0.0;
    }

  } else if (response_type == "proportion") {
    const double lo = 1e-9, hi = 1.0 - 1e-9;
    for (int i = 0; i < n; ++i) {
      const double bt = (w[i] == 1) ? betaT : 0.0;
      double v = y_base[i] + bt + eps[i];
      if (v < lo) v = lo;
      if (v > hi) v = hi;
      y[i] = v;
    }

  } else if (response_type == "count") {
    for (int i = 0; i < n; ++i) {
      const double bt  = (w[i] == 1) ? betaT : 0.0;
      const double lam = std::max(DBL_EPSILON, y_base[i] * std::exp(bt + eps[i]));
      y[i] = R::rpois(lam);
    }

  } else if (response_type == "survival") {
    NumericVector u_cens = runif(n);
    for (int i = 0; i < n; ++i) {
      const double bt = (w[i] == 1) ? betaT : 0.0;
      double t_i = std::max(DBL_EPSILON, y_base[i] * std::exp(bt + eps[i]));
      if (u_cens[i] < prob_censoring) {
        t_i    = R::runif(0.0, t_i);
        dead[i] = 0;
      }
      y[i] = t_i;
    }

  } else if (response_type == "ordinal") {
    for (int i = 0; i < n; ++i) {
      const double bt = (w[i] == 1) ? betaT : 0.0;
      int v = static_cast<int>(std::round(y_base[i] + bt + eps[i]));
      if (v < 1)               v = 1;
      if (v > n_ordinal_levels) v = n_ordinal_levels;
      y[i] = v;
    }

  } else {
    stop("apply_treatment_and_noise_cpp: unknown response_type '%s'.",
         response_type.c_str());
  }

  return List::create(_["y"] = y, _["dead"] = dead);
}
