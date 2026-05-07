#include <RcppEigen.h>
#include "_helper_functions.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector test_ols_start_beta_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
	return wrap(ols_start_beta(X, y));
}

// [[Rcpp::export]]
NumericVector test_ols_start_beta_on_log1p_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
	return wrap(ols_start_beta_on_log1p(X, y));
}

// [[Rcpp::export]]
NumericVector test_finalize_start_beta_cpp(const Eigen::VectorXd& smart_start,
                                           const Eigen::VectorXd& legacy_start,
                                           bool use_smart = true,
                                           Nullable<IntegerVector> fixed_idx = R_NilValue,
                                           Nullable<NumericVector> fixed_values = R_NilValue) {
	FixedParamSpec fixed_spec = make_fixed_param_spec(legacy_start.size(), fixed_idx, fixed_values);
	return wrap(finalize_start_beta(smart_start, legacy_start, fixed_spec, use_smart));
}

// [[Rcpp::export]]
List test_weibull_aft_start_cpp(const Eigen::MatrixXd& X,
                                const Eigen::VectorXd& y,
                                const Eigen::VectorXd& dead) {
	WeibullStart start = weibull_aft_start(X, y, dead);
	return List::create(
		_["beta"] = wrap(start.beta),
		_["log_sigma"] = start.log_sigma,
		_["params"] = wrap(weibull_start_to_params(start))
	);
}

// [[Rcpp::export]]
List test_ordinal_start_cpp(const Eigen::MatrixXd& X,
                            const Eigen::VectorXd& y,
                            std::string link = "logit") {
	edi_ordinal::Link ordinal_link = edi_ordinal::Link::Logit;
	if (link == "probit") ordinal_link = edi_ordinal::Link::Probit;
	else if (link == "cloglog") ordinal_link = edi_ordinal::Link::Cloglog;
	else if (link == "cauchit") ordinal_link = edi_ordinal::Link::Cauchit;
	OrdinalStart start = ordinal_start_from_ols(X, y, ordinal_link);
	return List::create(
		_["alpha"] = wrap(start.alpha),
		_["beta"] = wrap(start.beta),
		_["params"] = wrap(ordinal_start_to_params(start))
	);
}
