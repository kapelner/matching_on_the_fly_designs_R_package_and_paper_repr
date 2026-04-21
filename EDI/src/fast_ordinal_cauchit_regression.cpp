#include "_helper_functions.h"
#include "ordinal_fixed_link_helpers.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

class OrdinalCauchitRegression {
private:
    edi_ordinal::FixedOrdinalRegression m_model;

public:
    OrdinalCauchitRegression(const MatrixXd& X, const VectorXd& y) :
        m_model(X, y, edi_ordinal::Link::Cauchit, -1.0) {}

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
Eigen::VectorXd get_ordinal_cauchit_regression_score_cpp(const Eigen::MatrixXd& X,
														 const Eigen::VectorXd& y,
														 const Eigen::VectorXd& params) {
	OrdinalCauchitRegression model(X, y);
	Eigen::VectorXd grad(params.size());
	model(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_ordinal_cauchit_regression_hessian_cpp(const Eigen::MatrixXd& X,
														   const Eigen::VectorXd& y,
														   const Eigen::VectorXd& params) {
	OrdinalCauchitRegression model(X, y);
	return -model.hessian(params);
}

// [[Rcpp::export]]
List fast_ordinal_cauchit_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-6, std::string optimization_alg = "newton_raphson") {
    OrdinalCauchitRegression model(X, y);
    int p = X.cols();
    int K = OrdinalCauchitRegression::init_levels(y).size();
    if (K < 2) return List::create();
    int n_alpha = K - 1;
    int n_params = n_alpha + p;

    VectorXd params(n_params);
    // Initialize alpha using cauchit quantiles
    for (int k = 0; k < n_alpha; ++k) {
        double prob = static_cast<double>(k + 1) / static_cast<double>(K);
        params[k] = std::tan(M_PI * (prob - 0.5));
    }
    // Initialize beta
    params.tail(p).setZero();

    LikelihoodFitResult fit = optimize_likelihood(model, params, maxit, tol, optimization_alg, "newton_raphson");
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
List fast_ordinal_cauchit_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, std::string optimization_alg = "newton_raphson") {
    List res = fast_ordinal_cauchit_regression_cpp(X, y, 100, 1e-6, optimization_alg);
    if (res.size() == 0) return List::create(Named("b") = NumericVector::create(NA_REAL), Named("ssq_b_2") = NA_REAL);
    
    VectorXd params = res["params"];
    bool converged = res["converged"];
    OrdinalCauchitRegression model(X, y);
    MatrixXd H = model.hessian(params);
    
    FullPivLU<MatrixXd> lu(H);
    if (!lu.isInvertible()) {
        return List::create(Named("b") = res["b"], Named("ssq_b_2") = NA_REAL, Named("converged") = converged);
    }
    
    MatrixXd cov_mat = lu.inverse();
    int p = X.cols();
    int n_alpha = params.size() - p;
    
    double ssq_b_2 = (p >= 1) ? cov_mat(n_alpha, n_alpha) : NA_REAL;

    return List::create(
        Named("b") = res["b"],
        Named("ssq_b_2") = ssq_b_2,
        Named("converged") = converged
    );
}
