// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include <chrono>
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
double time_fast_trigamma_cpp(Rcpp::NumericVector x, int reps) {
	const int n = x.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) {
		for (int i = 0; i < n; ++i) {
			s += fast_trigamma(x[i]);
		}
	}
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");  // keep s live, never true
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_r_trigamma_cpp(Rcpp::NumericVector x, int reps) {
	const int n = x.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) {
		for (int i = 0; i < n; ++i) {
			s += R::trigamma(x[i]);
		}
	}
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}
