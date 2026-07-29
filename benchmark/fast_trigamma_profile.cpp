// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
double bench_fast_trigamma_scalar_loop_cpp(Rcpp::NumericVector x, int reps) {
	const int n = x.size();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) {
		for (int i = 0; i < n; ++i) {
			s += fast_trigamma(x[i] + 0.0000001 * r);
		}
	}
	return s;
}

// [[Rcpp::export]]
double bench_fast_trigamma_vec_loop_cpp(Rcpp::NumericVector x, int reps) {
	Eigen::Map<const Eigen::ArrayXd> xa(x.begin(), x.size());
	double s = 0.0;
	for (int r = 0; r < reps; ++r) {
		Eigen::ArrayXd shifted = xa + 0.0000001 * r;
		s += fast_trigamma_vec(shifted).sum();
	}
	return s;
}
