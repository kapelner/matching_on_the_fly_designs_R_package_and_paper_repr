// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(BH)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include <boost/math/special_functions/trigamma.hpp>
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
Rcpp::DataFrame fast_trigamma_vs_boost_cpp(Rcpp::NumericVector x) {
	const int n = x.size();
	Rcpp::NumericVector fast(n), boost_ref(n), r_ref(n);
	for (int i = 0; i < n; ++i) {
		fast[i] = fast_trigamma(x[i]);
		boost_ref[i] = boost::math::trigamma(x[i]);
		r_ref[i] = R::trigamma(x[i]);
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("x") = x,
		Rcpp::Named("fast") = fast,
		Rcpp::Named("boost") = boost_ref,
		Rcpp::Named("r") = r_ref
	);
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_trigamma_vec_vs_scalar_cpp(Rcpp::NumericVector x) {
	Eigen::Map<const Eigen::ArrayXd> xa(x.begin(), x.size());
	Eigen::ArrayXd out = fast_trigamma_vec(xa);
	return Rcpp::wrap(out);
}
