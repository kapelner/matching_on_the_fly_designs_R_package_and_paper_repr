// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include <chrono>
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
double time_fast_pchisq_cpp(Rcpp::NumericVector statistic, Rcpp::NumericVector df, int reps) {
	const int n = statistic.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r)
		for (int i = 0; i < n; ++i)
			s += fast_pchisq_upper(statistic[i], df[i]);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_r_pchisq_cpp(Rcpp::NumericVector statistic, Rcpp::NumericVector df, int reps) {
	const int n = statistic.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r)
		for (int i = 0; i < n; ++i)
			s += R::pchisq(statistic[i], df[i], false, false);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}
