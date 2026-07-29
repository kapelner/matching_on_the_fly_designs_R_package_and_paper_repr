// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(BH)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
Rcpp::DataFrame fast_qnorm_vs_boost_cpp(Rcpp::NumericVector p) {
	const int n = p.size();
	Rcpp::NumericVector fast(n), boost_ref(n), r_ref(n);
	boost::math::normal_distribution<double> std_normal(0.0, 1.0);
	for (int i = 0; i < n; ++i) {
		fast[i] = fast_qnorm(p[i]);
		boost_ref[i] = boost::math::quantile(std_normal, p[i]);
		r_ref[i] = R::qnorm5(p[i], 0.0, 1.0, 1, 0);
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("p") = p,
		Rcpp::Named("fast") = fast,
		Rcpp::Named("boost") = boost_ref,
		Rcpp::Named("r") = r_ref
	);
}

// [[Rcpp::export]]
Rcpp::DataFrame fast_log_pnorm_vs_boost_cpp(Rcpp::NumericVector x) {
	const int n = x.size();
	Rcpp::NumericVector fast(n), boost_ref(n), r_ref(n);
	boost::math::normal_distribution<double> std_normal(0.0, 1.0);
	for (int i = 0; i < n; ++i) {
		fast[i] = fast_log_pnorm(x[i]);
		boost_ref[i] = std::log(boost::math::cdf(std_normal, x[i]));
		r_ref[i] = R::pnorm5(x[i], 0.0, 1.0, 1, 1);
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("x") = x,
		Rcpp::Named("fast") = fast,
		Rcpp::Named("boost") = boost_ref,
		Rcpp::Named("r") = r_ref
	);
}

// [[Rcpp::export]]
Rcpp::DataFrame fast_lbeta_vs_boost_cpp(Rcpp::NumericVector a, Rcpp::NumericVector b) {
	const int n = a.size();
	Rcpp::NumericVector fast(n), boost_ref(n), r_ref(n);
	for (int i = 0; i < n; ++i) {
		fast[i] = fast_lbeta(a[i], b[i]);
		// log-space reference via boost::math::lgamma directly (not
		// log(boost::math::beta(...))): beta() computes the ratio of
		// (non-log) Gamma functions, which loses precision well before
		// double overflow for even moderately large a/b - lgamma-based
		// log-beta is the fair, numerically sound apples-to-apples
		// reference for what fast_lbeta itself computes.
		boost_ref[i] = boost::math::lgamma(a[i]) + boost::math::lgamma(b[i]) - boost::math::lgamma(a[i] + b[i]);
		r_ref[i] = R::lbeta(a[i], b[i]);
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("a") = a,
		Rcpp::Named("b") = b,
		Rcpp::Named("fast") = fast,
		Rcpp::Named("boost") = boost_ref,
		Rcpp::Named("r") = r_ref
	);
}

// [[Rcpp::export]]
Rcpp::DataFrame fast_dnbinom_mu_vs_r_cpp(Rcpp::NumericVector x, Rcpp::NumericVector size, Rcpp::NumericVector mu) {
	const int n = x.size();
	Rcpp::NumericVector fast_log(n), r_log(n), fast_p(n), r_p(n);
	for (int i = 0; i < n; ++i) {
		fast_log[i] = fast_dnbinom_mu(x[i], size[i], mu[i], true);
		r_log[i] = R::dnbinom_mu(x[i], size[i], mu[i], true);
		fast_p[i] = fast_dnbinom_mu(x[i], size[i], mu[i], false);
		r_p[i] = R::dnbinom_mu(x[i], size[i], mu[i], false);
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("x") = x,
		Rcpp::Named("size") = size,
		Rcpp::Named("mu") = mu,
		Rcpp::Named("fast_log") = fast_log,
		Rcpp::Named("r_log") = r_log,
		Rcpp::Named("fast_p") = fast_p,
		Rcpp::Named("r_p") = r_p
	);
}
