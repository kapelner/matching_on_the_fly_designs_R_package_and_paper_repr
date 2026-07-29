// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(BH)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include <boost/math/distributions/chi_squared.hpp>
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
Rcpp::DataFrame fast_pchisq_vs_boost_cpp(Rcpp::NumericVector statistic, Rcpp::NumericVector df) {
	const int n = statistic.size();
	Rcpp::NumericVector fast(n), boost_ref(n), r_ref(n);
	for (int i = 0; i < n; ++i) {
		fast[i] = fast_pchisq_upper(statistic[i], df[i]);
		boost::math::chi_squared_distribution<double> chisq(df[i]);
		boost_ref[i] = boost::math::cdf(boost::math::complement(chisq, statistic[i]));
		r_ref[i] = R::pchisq(statistic[i], df[i], false, false);
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("statistic") = statistic,
		Rcpp::Named("df") = df,
		Rcpp::Named("fast") = fast,
		Rcpp::Named("boost") = boost_ref,
		Rcpp::Named("r") = r_ref
	);
}
