#ifdef EDI_CORE_ONLY
#include "_helper_functions_core.h"
#include "result_map.h"
#else
#include "_helper_functions.h"
#include "result_map_rcpp.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#endif
#include "ordinal_fixed_link_helpers.h"
#include <algorithm>
#include <vector>
#include <optional>
#include <string>
#include <limits>
#include <stdexcept>

// [[Rcpp::depends(RcppEigen)]]

#ifndef EDI_CORE_ONLY
using namespace Rcpp;
#endif
using namespace Eigen;

class OrdinalCLLRegression {
private:
    edi_ordinal::FixedOrdinalRegression m_model;

public:
    OrdinalCLLRegression(const Eigen::Ref<const Eigen::MatrixXd>& X, const Eigen::Ref<const Eigen::VectorXd>& y) :
        m_model(X, y, edi_ordinal::Link::Cloglog, 1.0) {}

    static std::vector<double> init_levels(const Eigen::Ref<const Eigen::VectorXd>& y) {
        return edi_ordinal::init_levels(y);
    }

    double neg_log_likelihood(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        return m_model.neg_log_likelihood(params);
    }

    double operator()(const Eigen::Ref<const Eigen::VectorXd>& params, Eigen::Ref<Eigen::VectorXd> grad) const {
        return m_model(params, grad);
    }

    MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        return m_model.hessian(params);
    }

    MatrixXd expected_hessian(const Eigen::Ref<const Eigen::VectorXd>& params) const {
        return m_model.expected_hessian(params);
    }
};

#ifndef EDI_CORE_ONLY
// [[Rcpp::export]]
List fast_cloglog_link_eval_cpp(const NumericVector& x) {
    const R_xlen_t n = x.size();
    Eigen::VectorXd cdf(n), pdf(n), pdf_derivative(n);
    for (R_xlen_t i = 0; i < n; ++i) {
        double F = NA_REAL, f = NA_REAL, fp = NA_REAL;
        edi_ordinal::cdf_pdf_fpdf(edi_ordinal::Link::Cloglog, x[i], F, f, fp);
        cdf[i] = F;
        pdf[i] = f;
        pdf_derivative[i] = fp;
    }
    return edi::to_rcpp_list(edi::ResultMap()
        .set("cdf", cdf)
        .set("pdf", pdf)
        .set("pdf_derivative", pdf_derivative));
}

// [[Rcpp::export]]
SEXP get_ordinal_cloglog_regression_score_cpp(const Rcpp::NumericMatrix& X,
														 const Rcpp::NumericVector& y,
														 const Rcpp::NumericVector& params,
														 Nullable<IntegerVector> fixed_idx = R_NilValue,
														 Nullable<NumericVector> fixed_values = R_NilValue) {
	Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
	Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
	Eigen::Map<const Eigen::VectorXd> map_params(params.begin(), params.size());

	OrdinalCLLRegression model(map_X, map_y);
	FixedParamSpec fixed_spec = make_fixed_param_spec(
		map_params.size(),
		nullable_to_optional<Eigen::VectorXi>(fixed_idx),
		nullable_to_optional<Eigen::VectorXd>(fixed_values));
	Eigen::VectorXd par = apply_fixed_values(map_params, fixed_spec);
	Eigen::VectorXd grad(par.size());
	model(par, grad);
	return wrap(-grad);
}

// [[Rcpp::export]]
SEXP get_ordinal_cloglog_regression_hessian_cpp(const Rcpp::NumericMatrix& X,
														   const Rcpp::NumericVector& y,
														   const Rcpp::NumericVector& params,
														   Nullable<IntegerVector> fixed_idx = R_NilValue,
														   Nullable<NumericVector> fixed_values = R_NilValue) {
	Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
	Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());
	Eigen::Map<const Eigen::VectorXd> map_params(params.begin(), params.size());

	OrdinalCLLRegression model(map_X, map_y);
	FixedParamSpec fixed_spec = make_fixed_param_spec(
		map_params.size(),
		nullable_to_optional<Eigen::VectorXi>(fixed_idx),
		nullable_to_optional<Eigen::VectorXd>(fixed_values));
	Eigen::VectorXd par = apply_fixed_values(map_params, fixed_spec);
	return wrap(-model.hessian(par));
}
#endif // EDI_CORE_ONLY

// Portable core: see fast_ordinal_cauchit_regression.cpp for the analogous
// pattern this mirrors (fit-only when estimate_only, else fit + vcov/ssq_b_j
// folded in from the with_var variant below). K < 2 (degenerate outcome)
// returns an empty ResultMap, matching the original fit-only R export.
edi::ResultMap fast_ordinal_cloglog_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    std::optional<Eigen::VectorXd> warm_start_params = std::nullopt,
    bool smart_cold_start = true,
    int maxit = 100,
    double tol = 1e-6,
    std::string optimization_alg = "lbfgs",
    std::optional<Eigen::VectorXi> fixed_idx = std::nullopt,
    std::optional<Eigen::VectorXd> fixed_values = std::nullopt,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info = std::nullopt,
    bool estimate_only = false
) {
    OrdinalCLLRegression model(X, y);
    int p = (int)X.cols();
    int K = (int)OrdinalCLLRegression::init_levels(y).size();
    if (K < 2) return edi::ResultMap();
    int n_alpha = K - 1;
    int n_params = n_alpha + p;

    VectorXd params(n_params);
    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    if (warm_start_params.has_value()) {
        params = *warm_start_params;
        if (params.size() != n_params) throw std::invalid_argument("warm_start_params must have length equal to the number of model parameters");
    } else {
        OrdinalStart legacy_start;
        legacy_start.alpha = VectorXd(n_alpha);
        for (int k = 0; k < n_alpha; ++k) {
            legacy_start.alpha[k] = -1.0 + 2.0 * (k + 1) / K;
        }
        legacy_start.beta = VectorXd::Zero(p);
        params = ordinal_start_to_params(
            smart_cold_start ? ordinal_smart_cold_start_or_legacy(X, y, edi_ordinal::Link::Cloglog, legacy_start, fixed_spec)
                        : legacy_start
        );
    }

    params = apply_fixed_values(params, fixed_spec);

    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.has_value()) {
        info_start = *warm_start_fisher_info;
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(model, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs", 0, info_start_ptr);
    params = fit.params;

    if (estimate_only) {
        return edi::ResultMap()
            .set("b", params.tail(p))
            .set("alpha", params.head(n_alpha))
            .set("n_params", n_params)
            .set("params", params)
            .set("converged", fit.converged)
            .set("iterations", fit.niter);
    }

    MatrixXd H_full = model.hessian(params);

    double ssq_b_2 = std::numeric_limits<double>::quiet_NaN();
    MatrixXd vcov;
    bool have_vcov = false;
    if (fit.converged) {
        FixedParameterFunctor<OrdinalCLLRegression> fixed_obj(model, fixed_spec, params);
        VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
        MatrixXd H_free = fixed_obj.hessian(params_free);
        int n_alpha2 = n_params - p;
        int free_j = -1;
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == n_alpha2) { free_j = jj + 1; break; }
        if (p >= 1 && free_j > 0) ssq_b_2 = compute_diagonal_inverse_entry(H_free, free_j);
        MatrixXd cov_free = covariance_from_information(H_free);
        vcov = expand_free_covariance(n_params, fixed_spec, cov_free, true);
        have_vcov = true;
    }

    edi::ResultMap rm = edi::ResultMap()
        .set("b", params.tail(p))
        .set("alpha", params.head(n_alpha))
        .set("n_params", n_params)
        .set("params", params)
        .set("neg_loglik", fit.value)
        .set("converged", fit.converged)
        .set("iterations", fit.niter)
        .set("ssq_b_j", ssq_b_2)
        .set("observed_information", H_full)
        .set("fisher_information", H_full)
        .set("information", H_full)
        .set("information_type", std::string("observed"));
    if (have_vcov) {
        rm.set("vcov", vcov);
    } else {
        rm.set("vcov", std::monostate{});
    }
    return rm;
}

#ifndef EDI_CORE_ONLY
//' @title Fast Ordinal Cloglog Regression (C++)
//' @description High-performance ordinal cloglog regression fitting using Newton-Raphson.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param warm_start_params Optional starting values for [alpha, beta]. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param optimization_alg Optimization algorithm.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, thresholds, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ordinal_cloglog_regression_cpp(const Rcpp::NumericMatrix& X,
                                          const Rcpp::NumericVector& y,
                                          Nullable<NumericVector> warm_start_params = R_NilValue,
                                          bool smart_cold_start = true,
                                          int maxit = 100,
                                          double tol = 1e-6,
                                          std::string optimization_alg = "lbfgs",
                                          Nullable<IntegerVector> fixed_idx = R_NilValue,
                                          Nullable<NumericVector> fixed_values = R_NilValue,
                                          Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                          bool estimate_only = false) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());

    edi::ResultMap res = fast_ordinal_cloglog_regression_internal(
        map_X, map_y,
        nullable_to_optional<Eigen::VectorXd>(warm_start_params),
        smart_cold_start, maxit, tol, optimization_alg,
        nullable_to_optional<Eigen::VectorXi>(fixed_idx),
        nullable_to_optional<Eigen::VectorXd>(fixed_values),
        nullable_to_optional<Eigen::MatrixXd>(warm_start_fisher_info),
        estimate_only);
    return edi::to_rcpp_list(res);
}

//' @title Fast Ordinal Cloglog Regression with Variance (C++)
//' @description Ordinal cloglog regression fitting with full variance-covariance matrix.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param warm_start_params Optional starting values for [alpha, beta]. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param optimization_alg Optimization algorithm.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, thresholds, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_ordinal_cloglog_regression_with_var_cpp(const Rcpp::NumericMatrix& X,
                                                   const Rcpp::NumericVector& y,
                                                   Nullable<NumericVector> warm_start_params = R_NilValue,
                                                   bool smart_cold_start = true,
                                                   std::string optimization_alg = "lbfgs",
                                                   Nullable<IntegerVector> fixed_idx = R_NilValue,
                                                   Nullable<NumericVector> fixed_values = R_NilValue,
                                                   Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    Eigen::Map<const Eigen::MatrixXd> map_X(X.begin(), X.rows(), X.cols());
    Eigen::Map<const Eigen::VectorXd> map_y(y.begin(), y.size());

    edi::ResultMap res = fast_ordinal_cloglog_regression_internal(
        map_X, map_y,
        nullable_to_optional<Eigen::VectorXd>(warm_start_params),
        smart_cold_start, 100, 1e-6, optimization_alg,
        nullable_to_optional<Eigen::VectorXi>(fixed_idx),
        nullable_to_optional<Eigen::VectorXd>(fixed_values),
        nullable_to_optional<Eigen::MatrixXd>(warm_start_fisher_info),
        false);
    if (res.entries().empty()) return edi::to_rcpp_list(edi::ResultMap().set("b", NA_REAL).set("ssq_b_2", NA_REAL));
    return edi::to_rcpp_list(res);
}
#endif // EDI_CORE_ONLY
