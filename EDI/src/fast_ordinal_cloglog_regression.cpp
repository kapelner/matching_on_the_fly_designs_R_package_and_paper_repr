#include "_helper_functions.h"
#include "ordinal_fixed_link_helpers.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class OrdinalCLLRegression {
private:
    edi_ordinal::FixedOrdinalRegression m_model;

public:
    OrdinalCLLRegression(const MatrixXd& X, const VectorXd& y) :
        m_model(X, y, edi_ordinal::Link::Cloglog, 1.0) {}

    static std::vector<double> init_levels(const VectorXd& y) {
        return edi_ordinal::init_levels(y);
    }

    double neg_log_likelihood(const VectorXd& params) const {
        return m_model.neg_log_likelihood(params);
    }

    double operator()(const VectorXd& params, VectorXd& grad) const {
        return m_model(params, grad);
    }

    MatrixXd hessian(const VectorXd& params) const {
        return m_model.hessian(params);
    }
};

// [[Rcpp::export]]
Eigen::VectorXd get_ordinal_cloglog_regression_score_cpp(const Eigen::MatrixXd& X,
														 const Eigen::VectorXd& y,
														 const Eigen::VectorXd& params,
														 Nullable<IntegerVector> fixed_idx = R_NilValue,
														 Nullable<NumericVector> fixed_values = R_NilValue) {
	OrdinalCLLRegression model(X, y);
	FixedParamSpec fixed_spec = make_fixed_param_spec(params.size(), fixed_idx, fixed_values);
	Eigen::VectorXd par = apply_fixed_values(params, fixed_spec);
	Eigen::VectorXd grad(par.size());
	model(par, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_ordinal_cloglog_regression_hessian_cpp(const Eigen::MatrixXd& X,
														   const Eigen::VectorXd& y,
														   const Eigen::VectorXd& params,
														   Nullable<IntegerVector> fixed_idx = R_NilValue,
														   Nullable<NumericVector> fixed_values = R_NilValue) {
	OrdinalCLLRegression model(X, y);
	FixedParamSpec fixed_spec = make_fixed_param_spec(params.size(), fixed_idx, fixed_values);
	Eigen::VectorXd par = apply_fixed_values(params, fixed_spec);
	return -model.hessian(par);
}

// [[Rcpp::export]]
List fast_ordinal_cloglog_regression_cpp(const Eigen::MatrixXd& X, 
                                          const Eigen::VectorXd& y, 
                                          int maxit = 100, 
                                          double tol = 1e-6, 
                                          std::string optimization_alg = "newton_raphson",
                                          Nullable<IntegerVector> fixed_idx = R_NilValue,
                                          Nullable<NumericVector> fixed_values = R_NilValue) {
    OrdinalCLLRegression model(X, y);
    int p = X.cols();
    int K = OrdinalCLLRegression::init_levels(y).size();
    if (K < 2) return List::create();
    int n_alpha = K - 1;
    int n_params = n_alpha + p;

    VectorXd params(n_params);
    // Initialize alpha
    for (int k = 0; k < n_alpha; ++k) {
        params[k] = -1.0 + 2.0 * (k + 1) / K;
    }
    // Initialize beta
    params.tail(p).setZero();

    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    LikelihoodFitResult fit = optimize_fixed_likelihood(model, params, fixed_spec, maxit, tol, optimization_alg, "newton_raphson");
    params = fit.params;

    return List::create(
        Named("b") = params.tail(p),
        Named("alpha") = params.head(n_alpha),
        Named("n_params") = n_params,
        Named("params") = params,
        Named("neg_loglik") = fit.value,
        Named("converged") = fit.converged
    );
}

// [[Rcpp::export]]
List fast_ordinal_cloglog_regression_with_var_cpp(const Eigen::MatrixXd& X, 
                                                   const Eigen::VectorXd& y, 
                                                   std::string optimization_alg = "newton_raphson",
                                                   Nullable<IntegerVector> fixed_idx = R_NilValue,
                                                   Nullable<NumericVector> fixed_values = R_NilValue) {
    List res = fast_ordinal_cloglog_regression_cpp(X, y, 100, 1e-6, optimization_alg, fixed_idx, fixed_values);
    if (res.size() == 0) return List::create(Named("b") = NumericVector::create(NA_REAL), Named("ssq_b_2") = NA_REAL);
    
    VectorXd params = res["params"];
    bool converged = res["converged"];
    OrdinalCLLRegression model(X, y);
    int n_params = params.size();
    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    
    double ssq_b_2 = NA_REAL;
    if (converged) {
        FixedParameterFunctor<OrdinalCLLRegression> fixed_obj(model, fixed_spec, params);
        VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
        MatrixXd H_free = fixed_obj.hessian(params_free);
        
        FullPivLU<MatrixXd> lu(H_free);
        if (lu.isInvertible()) {
            MatrixXd inv_free = lu.inverse();
            MatrixXd cov_mat = expand_free_covariance(n_params, fixed_spec, inv_free, true);
            int p = X.cols();
            int n_alpha = n_params - p;
            if (p >= 1) ssq_b_2 = cov_mat(n_alpha, n_alpha);
        }
    }

    return List::create(
        Named("b") = res["b"],
        Named("alpha") = res["alpha"],
        Named("params") = params,
        Named("neg_loglik") = res["neg_loglik"],
        Named("ssq_b_2") = ssq_b_2,
        Named("converged") = converged
    );
}
