// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <RcppEigen.h>
#include <chrono>
#include "../EDI/src/_helper_functions.h"

// [[Rcpp::export]]
double time_fast_qnorm_cpp(Rcpp::NumericVector p, int reps) {
	const int n = p.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += fast_qnorm(p[i]);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_r_qnorm_cpp(Rcpp::NumericVector p, int reps) {
	const int n = p.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += R::qnorm5(p[i], 0.0, 1.0, 1, 0);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_fast_log_pnorm_cpp(Rcpp::NumericVector x, int reps) {
	const int n = x.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += fast_log_pnorm(x[i]);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_r_pnorm_cpp(Rcpp::NumericVector x, int reps) {
	const int n = x.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += R::pnorm5(x[i], 0.0, 1.0, 1, 1);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_fast_lbeta_cpp(Rcpp::NumericVector a, Rcpp::NumericVector b, int reps) {
	const int n = a.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += fast_lbeta(a[i], b[i]);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_r_lbeta_cpp(Rcpp::NumericVector a, Rcpp::NumericVector b, int reps) {
	const int n = a.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += R::lbeta(a[i], b[i]);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_fast_dnbinom_mu_cpp(Rcpp::NumericVector x, Rcpp::NumericVector size, Rcpp::NumericVector mu, int reps) {
	const int n = x.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += fast_dnbinom_mu(x[i], size[i], mu[i], true);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}

// [[Rcpp::export]]
double time_r_dnbinom_mu_cpp(Rcpp::NumericVector x, Rcpp::NumericVector size, Rcpp::NumericVector mu, int reps) {
	const int n = x.size();
	auto start = std::chrono::steady_clock::now();
	double s = 0.0;
	for (int r = 0; r < reps; ++r) for (int i = 0; i < n; ++i) s += R::dnbinom_mu(x[i], size[i], mu[i], true);
	auto end = std::chrono::steady_clock::now();
	if (s != s) Rcpp::stop("nan");
	return std::chrono::duration<double, std::milli>(end - start).count();
}
