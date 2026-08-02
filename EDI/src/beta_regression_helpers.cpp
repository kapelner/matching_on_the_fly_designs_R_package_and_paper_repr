#ifdef EDI_CORE_ONLY
#include <Eigen/Dense>
#else
#include <RcppEigen.h>
#endif
#include <cmath>
#include <limits>
#include "fast_gamma_functions.h"

// [[Rcpp::depends(RcppEigen)]]

double beta_loglik_internal(const Eigen::Ref<const Eigen::VectorXd>& y,
                             const Eigen::Ref<const Eigen::VectorXd>& mu,
                             const double phi,
                             const Eigen::Ref<const Eigen::VectorXd>& wt) {
	const int n = y.size();
	double out = 0.0;

	for (int i = 0; i < n; i++) {
	const double mui = mu[i];
	const double yi  = y[i];
	const double wi  = wt[i];

	// clamp to avoid log(0)
	double mui_phi = mui * phi;
	double one_minus_mui_phi = (1.0 - mui) * phi;

	if (mui <= 0.0) mui_phi = 1e-12;
	if (mui >= 1.0) one_minus_mui_phi = 1e-12;

	double term =
		fast_lgamma(phi) -
		fast_lgamma(mui_phi) -
		fast_lgamma(one_minus_mui_phi) +
		(mui_phi - 1.0) * std::log(yi) +
		(one_minus_mui_phi - 1.0) * std::log1p(-yi);

	out += wi * term;
	}

	return out;
}

Eigen::VectorXd beta_dev_resids_internal(const Eigen::Ref<const Eigen::VectorXd>& y,
                                          const Eigen::Ref<const Eigen::VectorXd>& mu,
                                          const double phi,
                                          const Eigen::Ref<const Eigen::VectorXd>& wt) {
	const int n = y.size();
	Eigen::VectorXd out(n);

	for (int i = 0; i < n; i++) {
	const double mui = mu[i];
	const double yi  = y[i];
	const double wi  = wt[i];

	double mui_phi = mui * phi;
	double one_minus_mui_phi = (1.0 - mui) * phi;

	// numerical safeguard
	if (mui <= 0.0) mui_phi = 1e-12;
	if (mui >= 1.0) one_minus_mui_phi = 1e-12;
	if (yi <= 0.0) out[i] = std::numeric_limits<double>::quiet_NaN();  // avoid log(0)
	else if (yi >= 1.0) out[i] = std::numeric_limits<double>::quiet_NaN();
	else {
		double val = fast_lbeta(mui_phi, one_minus_mui_phi) -
		(mui_phi - 1.0) * std::log(yi) -
		(one_minus_mui_phi - 1.0) * std::log1p(-yi);

		out[i] = 2.0 * wi * val;
	}
	}

	return out;
}

double beta_aic_internal(const Eigen::Ref<const Eigen::VectorXd>& y,
                          const Eigen::Ref<const Eigen::VectorXd>& mu,
                          const double phi,
                          const Eigen::Ref<const Eigen::VectorXd>& wt) {
	const int n = y.size();
	double ll = 0.0;

	for (int i = 0; i < n; i++) {
	const double yi  = y[i];
	const double mui = mu[i];
	const double wi  = wt[i];

	double mui_phi = mui * phi;
	double one_minus_mui_phi = (1.0 - mui) * phi;

	// numerical safeguards
	if (mui <= 0.0) mui_phi = 1e-12;
	if (mui >= 1.0) one_minus_mui_phi = 1e-12;

	if (yi <= 0.0 || yi >= 1.0) continue; // skip invalid

	double term =
		fast_lgamma(phi) -
		fast_lgamma(mui_phi) -
		fast_lgamma(one_minus_mui_phi) +
		(mui_phi - 1.0) * std::log(yi) +
		(one_minus_mui_phi - 1.0) * std::log1p(-yi);

	ll += wi * term;
	}

	// AIC = -2 logLik + 2k
	return -2.0 * ll + 2.0 * (mu.size() + 1);
}

#ifndef EDI_CORE_ONLY
// [[Rcpp::export]]
double beta_loglik_cpp(SEXP y_sexp, SEXP mu_sexp, const double phi, SEXP wt_sexp) {
	Rcpp::NumericVector y_r(y_sexp);
	Rcpp::NumericVector mu_r(mu_sexp);
	Rcpp::NumericVector wt_r(wt_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> mu(mu_r.begin(), mu_r.size());
	Eigen::Map<const Eigen::VectorXd> wt(wt_r.begin(), wt_r.size());
	return beta_loglik_internal(y, mu, phi, wt);
}

// [[Rcpp::export]]
Eigen::VectorXd beta_dev_resids_cpp(SEXP y_sexp, SEXP mu_sexp, const double phi, SEXP wt_sexp) {
	Rcpp::NumericVector y_r(y_sexp);
	Rcpp::NumericVector mu_r(mu_sexp);
	Rcpp::NumericVector wt_r(wt_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> mu(mu_r.begin(), mu_r.size());
	Eigen::Map<const Eigen::VectorXd> wt(wt_r.begin(), wt_r.size());
	return beta_dev_resids_internal(y, mu, phi, wt);
}

// [[Rcpp::export]]
double beta_aic_cpp(SEXP y_sexp, SEXP mu_sexp, const double phi, SEXP wt_sexp) {
	Rcpp::NumericVector y_r(y_sexp);
	Rcpp::NumericVector mu_r(mu_sexp);
	Rcpp::NumericVector wt_r(wt_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> mu(mu_r.begin(), mu_r.size());
	Eigen::Map<const Eigen::VectorXd> wt(wt_r.begin(), wt_r.size());
	return beta_aic_internal(y, mu, phi, wt);
}
#endif // EDI_CORE_ONLY
