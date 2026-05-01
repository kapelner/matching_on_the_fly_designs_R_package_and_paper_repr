#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <string>

using namespace Rcpp;

double expit(double u){
  return 1.0 / (1.0 + std::exp(-u));
}

// [[Rcpp::export]]
List apply_treatment_and_noise_cpp(
    const NumericVector& y_linear_model,
    const IntegerVector& w,
    const std::string&   response_type,
    double               betaT,
    double               sd_noise,
    double               prob_censoring,
    int                  n_ordinal_levels,
    double               phi_proportion,
    double               k_survival,
    double               incidence_clamp = 1e-9,
    double               proportion_clamp = 1e-9,
    double               count_clamp = 1e-9,
    double               survival_clamp = 1e-9)
{
  const int n = y_linear_model.size();
  if (w.size() != n)
    stop("apply_treatment_and_noise_cpp: y_linear_model and w must have the same length.");
  if (!std::isfinite(phi_proportion) || phi_proportion <= 0.0)
    stop("apply_treatment_and_noise_cpp: phi_proportion must be finite and > 0.");
  if (!std::isfinite(k_survival) || k_survival <= 0.0)
    stop("apply_treatment_and_noise_cpp: k_survival must be finite and > 0.");
  if (!std::isfinite(incidence_clamp) || incidence_clamp <= 0.0 || incidence_clamp >= 0.5)
    stop("apply_treatment_and_noise_cpp: incidence_clamp must be finite and in (0, 0.5).");
  if (!std::isfinite(proportion_clamp) || proportion_clamp <= 0.0 || proportion_clamp >= 0.5)
    stop("apply_treatment_and_noise_cpp: proportion_clamp must be finite and in (0, 0.5).");
  if (!std::isfinite(count_clamp) || count_clamp <= 0.0)
    stop("apply_treatment_and_noise_cpp: count_clamp must be finite and > 0.");
  if (!std::isfinite(survival_clamp) || survival_clamp <= 0.0)
    stop("apply_treatment_and_noise_cpp: survival_clamp must be finite and > 0.");
  NumericVector y(n);
  IntegerVector dead(n, 1);   // default: not censored

  if (response_type == "continuous") {
    // Pre-generate all Gaussian noise in one RNG call
    NumericVector eps = rnorm(n, 0.0, sd_noise);
    for (int i = 0; i < n; ++i) {
      const double bt_i = (w[i] == 1) ? betaT : 0.0;
      y[i] = y_linear_model[i] + bt_i + eps[i];
    }

  } else if (response_type == "incidence") {
    // Pre-generate all uniform draws for Bernoulli in one call
    NumericVector u = runif(n);
    for (int i = 0; i < n; ++i) {
      const double bt_i  = (w[i] == 1) ? betaT : 0.0;
      double p_i = expit(y_linear_model[i] + bt_i);
      p_i = std::min(1.0 - incidence_clamp, std::max(incidence_clamp, p_i));
      if (!std::isfinite(p_i) || p_i <= 0.0 || p_i >= 1.0)
        stop("apply_treatment_and_noise_cpp: incidence Bernoulli probability must be finite and in (0, 1).");
      y[i] = (u[i] < p_i) ? 1.0 : 0.0;
    }

  } else if (response_type == "proportion") {
    for (int i = 0; i < n; ++i) {
      const double bt_i = (w[i] == 1) ? betaT : 0.0;
      double mu_i = expit(y_linear_model[i] + bt_i);
      mu_i = std::min(1.0 - proportion_clamp, std::max(proportion_clamp, mu_i));
      if (!std::isfinite(mu_i) || mu_i <= 0.0 || mu_i >= 1.0)
        stop("apply_treatment_and_noise_cpp: proportion beta mean must be finite and in (0, 1).");
      y[i] = R::rbeta(mu_i * phi_proportion, (1.0 - mu_i) * phi_proportion);
    }

  } else if (response_type == "count") {
    for (int i = 0; i < n; ++i) {
      const double bt_i  = (w[i] == 1) ? betaT : 0.0;
      double mu_i = std::exp(y_linear_model[i] + bt_i);
      if (!std::isfinite(mu_i))
        stop("apply_treatment_and_noise_cpp: count Poisson mean must be finite.");
      mu_i = std::max(count_clamp, mu_i);
      y[i] = R::rpois(mu_i);
    }

  } else if (response_type == "survival") {
    NumericVector u_cens = runif(n);
    for (int i = 0; i < n; ++i) {
      const double bt_i = (w[i] == 1) ? betaT : 0.0;
      double mu_i = std::exp(y_linear_model[i] + bt_i);
      if (!std::isfinite(mu_i))
        stop("apply_treatment_and_noise_cpp: survival Weibull shape must be finite.");
      mu_i = std::max(survival_clamp, mu_i);
      double t_i = R::rweibull(mu_i, k_survival); //no epsilon necessary as we draw from the uniform to generate randomness
      if (u_cens[i] < prob_censoring) {
        t_i     = R::runif(0.0, t_i);
        dead[i] = 0;
      }
      y[i] = t_i;
    }

  } else if (response_type == "ordinal") {
    // Pre-generate all Gaussian noise in one RNG call
    NumericVector eps = rnorm(n, 0.0, sd_noise);
    for (int i = 0; i < n; ++i) {
      const double bt_i = (w[i] == 1) ? betaT : 0.0;
      int v = static_cast<int>(std::round(y_linear_model[i] + bt_i + eps[i]));
      if (v < 1)                v = 1;
      if (v > n_ordinal_levels) v = n_ordinal_levels;
      y[i] = v;
    }

  } else {
    stop("apply_treatment_and_noise_cpp: unknown response_type '%s'.",
         response_type.c_str());
  }

  return List::create(_["y"] = y, _["dead"] = dead);
}
