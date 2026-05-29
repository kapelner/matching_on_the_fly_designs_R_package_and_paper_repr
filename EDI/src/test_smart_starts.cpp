#include <RcppEigen.h>
#include "_helper_functions.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector test_ols_smart_cold_start_beta_cpp(SEXP X_sexp, SEXP y_sexp) {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	return wrap(ols_smart_cold_start_beta(X, y));
}

// [[Rcpp::export]]
NumericVector test_ols_smart_cold_start_beta_on_log1p_cpp(SEXP X_sexp, SEXP y_sexp) {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	return wrap(ols_smart_cold_start_beta_on_log1p(X, y));
}

// [[Rcpp::export]]
NumericVector test_finalize_warm_start_beta_cpp(SEXP smart_cold_start_sexp,
                                           SEXP legacy_start_sexp,
                                           bool use_smart = true,
                                           Nullable<IntegerVector> fixed_idx = R_NilValue,
                                           Nullable<NumericVector> fixed_values = R_NilValue) {
	NumericVector smart_cold_start_r(smart_cold_start_sexp); NumericVector legacy_start_r(legacy_start_sexp);
	Eigen::Map<const Eigen::VectorXd> smart_cold_start(smart_cold_start_r.begin(), smart_cold_start_r.size());
	Eigen::Map<const Eigen::VectorXd> legacy_start(legacy_start_r.begin(), legacy_start_r.size());
	FixedParamSpec fixed_spec = make_fixed_param_spec(legacy_start.size(), fixed_idx, fixed_values);
	return wrap(finalize_warm_start_beta(smart_cold_start, legacy_start, fixed_spec, use_smart));
}

// [[Rcpp::export]]
List test_weibull_aft_start_cpp(SEXP X_sexp, SEXP y_sexp, SEXP dead_sexp) {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp); NumericVector dead_r(dead_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXd> dead(dead_r.begin(), dead_r.size());
	WeibullStart warm_start_params = weibull_aft_start(X, y, dead);
	return List::create(
		_["beta"] = wrap(warm_start_params.beta),
		_["log_sigma"] = warm_start_params.log_sigma,
		_["params"] = wrap(weibull_start_to_params(warm_start_params))
	);
}

// [[Rcpp::export]]
List test_ordinal_start_cpp(SEXP X_sexp, SEXP y_sexp,
                            std::string link = "logit") {
	NumericMatrix X_r(X_sexp); NumericVector y_r(y_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	edi_ordinal::Link ordinal_link = edi_ordinal::Link::Logit;
	if (link == "probit") ordinal_link = edi_ordinal::Link::Probit;
	else if (link == "cloglog") ordinal_link = edi_ordinal::Link::Cloglog;
	else if (link == "cauchit") ordinal_link = edi_ordinal::Link::Cauchit;
	OrdinalStart warm_start_params = ordinal_smart_cold_start(X, y, ordinal_link);
	return List::create(
		_["alpha"] = wrap(warm_start_params.alpha),
		_["beta"] = wrap(warm_start_params.beta),
		_["params"] = wrap(ordinal_start_to_params(warm_start_params))
	);
}
