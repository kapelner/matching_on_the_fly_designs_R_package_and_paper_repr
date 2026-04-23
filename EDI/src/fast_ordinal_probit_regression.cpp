#include "_helper_functions.h"
#include "ordinal_fixed_link_helpers.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class OrdinalProbitRegression {
private:
    edi_ordinal::FixedOrdinalRegression m_model;

public:
    OrdinalProbitRegression(const MatrixXd& X, const VectorXd& y) :
        m_model(X, y, edi_ordinal::Link::Probit, -1.0) {}

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

static MatrixXd pseudo_inverse_symmetric_probit(const MatrixXd& A, double tol = 1e-8) {
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    VectorXd sing = svd.singularValues();
    MatrixXd Dinv = MatrixXd::Zero(sing.size(), sing.size());
    for (int i = 0; i < sing.size(); ++i) {
        if (sing[i] > tol) {
            Dinv(i, i) = 1.0 / sing[i];
        }
    }
    return svd.matrixV() * Dinv * svd.matrixU().transpose();
}

// [[Rcpp::export]]
Eigen::VectorXd get_ordinal_probit_regression_score_cpp(const Eigen::MatrixXd& X,
														const Eigen::VectorXd& y,
														const Eigen::VectorXd& params,
														Nullable<IntegerVector> fixed_idx = R_NilValue,
														Nullable<NumericVector> fixed_values = R_NilValue) {
	OrdinalProbitRegression model(X, y);
	FixedParamSpec fixed_spec = make_fixed_param_spec(params.size(), fixed_idx, fixed_values);
	Eigen::VectorXd par = apply_fixed_values(params, fixed_spec);
	Eigen::VectorXd grad(par.size());
	model(par, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_ordinal_probit_regression_hessian_cpp(const Eigen::MatrixXd& X,
														  const Eigen::VectorXd& y,
														  const Eigen::VectorXd& params,
														  Nullable<IntegerVector> fixed_idx = R_NilValue,
														  Nullable<NumericVector> fixed_values = R_NilValue) {
	OrdinalProbitRegression model(X, y);
	FixedParamSpec fixed_spec = make_fixed_param_spec(params.size(), fixed_idx, fixed_values);
	Eigen::VectorXd par = apply_fixed_values(params, fixed_spec);
	return -model.hessian(par);
}

// [[Rcpp::export]]
List fast_ordinal_probit_regression_cpp(const Eigen::MatrixXd& X, 
                                         const Eigen::VectorXd& y, 
                                         int maxit = 100, 
                                         double tol = 1e-6, 
                                         std::string optimization_alg = "newton_raphson",
                                         Nullable<IntegerVector> fixed_idx = R_NilValue,
                                         Nullable<NumericVector> fixed_values = R_NilValue) {
    OrdinalProbitRegression model(X, y);
    const int p = X.cols();
    const int K = OrdinalProbitRegression::init_levels(y).size();
    const int n_alpha = K - 1;
    const int n_params = n_alpha + p;

    VectorXd params(n_params);
    for (int k = 0; k < n_alpha; ++k) {
        double q = static_cast<double>(k + 1) / static_cast<double>(K);
        params[k] = R::qnorm5(q, 0.0, 1.0, 1, 0);
    }
    params.tail(p).setZero();

    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    LikelihoodFitResult fit = optimize_fixed_likelihood(model, params, fixed_spec, maxit, tol, optimization_alg, "newton_raphson");
    params = fit.params;

    return List::create(
        Named("b") = params.tail(p),
        Named("alpha") = params.head(n_alpha),
        Named("n_params") = n_params,
        Named("params") = params,
        Named("converged") = fit.converged
    );
}

// [[Rcpp::export]]
List fast_ordinal_probit_regression_with_var_cpp(const Eigen::MatrixXd& X, 
                                                  const Eigen::VectorXd& y, 
                                                  std::string optimization_alg = "newton_raphson",
                                                  Nullable<IntegerVector> fixed_idx = R_NilValue,
                                                  Nullable<NumericVector> fixed_values = R_NilValue) {
    List res = fast_ordinal_probit_regression_cpp(X, y, 100, 1e-6, optimization_alg, fixed_idx, fixed_values);
    VectorXd params = res["params"];
    bool converged = res["converged"];
    OrdinalProbitRegression model(X, y);
    int n_params = params.size();
    FixedParamSpec fixed_spec = make_fixed_param_spec(n_params, fixed_idx, fixed_values);
    
    double ssq_b_2 = NA_REAL;
    if (converged) {
        FixedParameterFunctor<OrdinalProbitRegression> fixed_obj(model, fixed_spec, params);
        VectorXd params_free = subset_vector(params, fixed_spec.free_idx);
        MatrixXd H_free = fixed_obj.hessian(params_free);
        MatrixXd inv_free;

        FullPivLU<MatrixXd> lu(H_free);
        if (lu.isInvertible()) {
            inv_free = lu.inverse();
        } else {
            inv_free = pseudo_inverse_symmetric_probit(H_free);
        }

        MatrixXd cov_mat = expand_free_covariance(n_params, fixed_spec, inv_free, true);
        const int p = X.cols();
        const int n_alpha = n_params - p;
        ssq_b_2 = (p >= 1) ? cov_mat(n_alpha, n_alpha) : NA_REAL;
        if (!R_finite(ssq_b_2) || ssq_b_2 <= 0) {
            ssq_b_2 = NA_REAL;
        }
    }

    return List::create(
        Named("b") = res["b"],
        Named("ssq_b_2") = ssq_b_2,
        Named("converged") = converged
    );
}
