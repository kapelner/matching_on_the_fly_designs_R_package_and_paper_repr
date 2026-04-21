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
														const Eigen::VectorXd& params) {
	OrdinalProbitRegression model(X, y);
	Eigen::VectorXd grad(params.size());
	model(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_ordinal_probit_regression_hessian_cpp(const Eigen::MatrixXd& X,
														  const Eigen::VectorXd& y,
														  const Eigen::VectorXd& params) {
	OrdinalProbitRegression model(X, y);
	return -model.hessian(params);
}

// [[Rcpp::export]]
List fast_ordinal_probit_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-6, std::string optimization_alg = "newton_raphson") {
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
List fast_ordinal_probit_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, std::string optimization_alg = "newton_raphson") {
    List res = fast_ordinal_probit_regression_cpp(X, y, 100, 1e-6, optimization_alg);
    VectorXd params = res["params"];
    bool converged = res["converged"];
    OrdinalProbitRegression model(X, y);
    MatrixXd H = model.hessian(params);
    MatrixXd cov_mat;

    FullPivLU<MatrixXd> lu(H);
    if (lu.isInvertible()) {
        cov_mat = lu.inverse();
    } else {
        cov_mat = pseudo_inverse_symmetric_probit(H);
    }

    const int p = X.cols();
    const int n_alpha = params.size() - p;
    double ssq_b_2 = (p >= 1) ? cov_mat(n_alpha, n_alpha) : NA_REAL;
    if (!R_finite(ssq_b_2) || ssq_b_2 <= 0) {
        ssq_b_2 = NA_REAL;
    }

    return List::create(
        Named("b") = res["b"],
        Named("ssq_b_2") = ssq_b_2,
        Named("converged") = converged
    );
}
