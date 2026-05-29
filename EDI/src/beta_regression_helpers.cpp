#include <RcppEigen.h>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double beta_loglik_cpp(SEXP y_sexp,
					 SEXP mu_sexp,
					 const double phi,
					 SEXP wt_sexp) {
	Rcpp::NumericVector y_r(y_sexp);
	Rcpp::NumericVector mu_r(mu_sexp);
	Rcpp::NumericVector wt_r(wt_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> mu(mu_r.begin(), mu_r.size());
	Eigen::Map<const Eigen::VectorXd> wt(wt_r.begin(), wt_r.size());
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
		R::lgammafn(phi) -
		R::lgammafn(mui_phi) -
		R::lgammafn(one_minus_mui_phi) +
		(mui_phi - 1.0) * std::log(yi) +
		(one_minus_mui_phi - 1.0) * std::log1p(-yi);

	out += wi * term;
	}

	return out;
}

// [[Rcpp::export]]
Eigen::VectorXd beta_dev_resids_cpp(SEXP y_sexp,
									SEXP mu_sexp,
									const double phi,
									SEXP wt_sexp) {
	Rcpp::NumericVector y_r(y_sexp);
	Rcpp::NumericVector mu_r(mu_sexp);
	Rcpp::NumericVector wt_r(wt_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> mu(mu_r.begin(), mu_r.size());
	Eigen::Map<const Eigen::VectorXd> wt(wt_r.begin(), wt_r.size());
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
	if (yi <= 0.0) out[i] = NA_REAL;  // avoid log(0)
	else if (yi >= 1.0) out[i] = NA_REAL;
	else {
		double val = R::lbeta(mui_phi, one_minus_mui_phi) -
		(mui_phi - 1.0) * std::log(yi) -
		(one_minus_mui_phi - 1.0) * std::log1p(-yi);

		out[i] = 2.0 * wi * val;
	}
	}

	return out;
}

// [[Rcpp::export]]
double beta_aic_cpp(SEXP y_sexp,
					SEXP mu_sexp,
					const double phi,
					SEXP wt_sexp) {
	Rcpp::NumericVector y_r(y_sexp);
	Rcpp::NumericVector mu_r(mu_sexp);
	Rcpp::NumericVector wt_r(wt_sexp);
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> mu(mu_r.begin(), mu_r.size());
	Eigen::Map<const Eigen::VectorXd> wt(wt_r.begin(), wt_r.size());
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
		R::lgammafn(phi) -
		R::lgammafn(mui_phi) -
		R::lgammafn(one_minus_mui_phi) +
		(mui_phi - 1.0) * std::log(yi) +
		(one_minus_mui_phi - 1.0) * std::log1p(-yi);

	ll += wi * term;
	}

	// AIC = -2 logLik + 2k
	return -2.0 * ll + 2.0 * (mu.size() + 1);
}
