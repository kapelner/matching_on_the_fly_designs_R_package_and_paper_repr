#include "_helper_functions.h"
#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace LBFGSpp;

namespace {

struct DigammaFunctor {
	double operator()(double x) const {
		return R::digamma(x);
	}
};

class BetaRegression {
private:
	const Eigen::VectorXd m_y;
	const Eigen::MatrixXd m_X;
	const int m_n;
	const int m_p;
	const Eigen::VectorXd m_log_y;
	const Eigen::VectorXd m_log1_y;

public:
	BetaRegression(const Eigen::VectorXd& y, const Eigen::MatrixXd& X) :
		m_y(y), m_X(X), m_n(X.rows()), m_p(X.cols()),
		m_log_y(y.array().log().matrix()),
		m_log1_y((1.0 - y.array()).log().matrix()) {}

	double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad) {
		Eigen::VectorXd beta = params.head(m_p);
		double log_phi = params[m_p];
		double phi = std::exp(log_phi);

		Eigen::VectorXd eta = m_X * beta;
		Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();

		double epsilon = 1e-8;
		for(int i=0; i<m_n; ++i) {
			if (mu[i] < epsilon) mu[i] = epsilon;
			if (mu[i] > 1.0 - epsilon) mu[i] = 1.0 - epsilon;
		}

		double neg_ll = - (
			(m_n * R::lgammafn(phi)) -
			(mu.array() * phi).unaryExpr([](double x){ return R::lgammafn(x); }).sum() -
			((1.0 - mu.array()) * phi).unaryExpr([](double x){ return R::lgammafn(x); }).sum() +
			((mu.array() * phi - 1.0) * m_log_y.array()).sum() +
			(((1.0 - mu.array()) * phi - 1.0) * m_log1_y.array()).sum()
		);

		grad.resize(m_p + 1);
		Eigen::VectorXd mu_phi = mu.array() * phi;
		Eigen::VectorXd one_minus_mu_phi = (1.0 - mu.array()) * phi;
		Eigen::VectorXd d_mu_d_eta = mu.array() * (1.0 - mu.array());

		Eigen::VectorXd d_neg_ll_d_mu = (
			mu_phi.unaryExpr(DigammaFunctor()).array() -
			one_minus_mu_phi.unaryExpr(DigammaFunctor()).array() -
			m_log_y.array() + m_log1_y.array()
		) * phi;

		grad.head(m_p) = m_X.transpose() * (d_neg_ll_d_mu.array() * d_mu_d_eta.array()).matrix();

		double d_neg_ll_d_phi = (
			-m_n * R::digamma(phi) +
			(mu.array() * mu_phi.unaryExpr(DigammaFunctor()).array()).sum() +
			((1.0 - mu.array()) * one_minus_mu_phi.unaryExpr(DigammaFunctor()).array()).sum() -
			(mu.array() * m_log_y.array()).sum() -
			((1.0 - mu.array()) * m_log1_y.array()).sum()
		);
		grad[m_p] = d_neg_ll_d_phi * phi;

		return neg_ll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& params) {
		int total_p = m_p + 1;
		Eigen::MatrixXd H(total_p, total_p);
		H.setZero();
		double h = 1e-6;
		Eigen::VectorXd grad_at_params(total_p);
		operator()(params, grad_at_params);

		for (int i = 0; i < total_p; ++i) {
			Eigen::VectorXd p_plus = params;
			p_plus[i] += h;
			Eigen::VectorXd g_plus(total_p);
			operator()(p_plus, g_plus);
			H.col(i) = (g_plus - grad_at_params) / h;
		}
		H = (H + H.transpose()) / 2.0;
		return H;
	}
};

ModelResult fast_beta_regression_internal(const Eigen::MatrixXd& X,
                                        const Eigen::VectorXd& y,
                                        const Eigen::VectorXd* start_beta = nullptr,
                                        double start_phi = 10.0) {
    int p = X.cols();
    ModelResult res;
    Eigen::VectorXd params(p + 1);
    if (start_beta) params.head(p) = *start_beta;
    else params.head(p).setZero();
    params[p] = std::log(start_phi);

    BetaRegression fun(y, X);
    LBFGSParam<double> lbfgs_params;
    lbfgs_params.epsilon = 1e-6;
    lbfgs_params.max_iterations = 1000;

    LBFGSSolver<double> solver(lbfgs_params);
    double neg_ll;
    int niter = solver.minimize(fun, params, neg_ll);

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // phi
    res.XtWX = fun.hessian(params); // This is actually the Hessian H
    res.converged = (niter < 1000);
    return res;
}

} // namespace

// [[Rcpp::export]]
List fast_beta_regression_cpp(const Eigen::MatrixXd& X,
								const NumericVector& y,
								Nullable<NumericVector> start_beta = R_NilValue,
								double start_phi = 10.0,
								bool compute_std_errs = false) {

	Eigen::VectorXd y_eigen = as<Eigen::VectorXd>(y);
    Eigen::VectorXd sb;
    Eigen::VectorXd* sb_ptr = nullptr;
    if (start_beta.isNotNull()) {
        sb = as<Eigen::VectorXd>(start_beta);
        sb_ptr = &sb;
    }

    ModelResult fit = fast_beta_regression_internal(X, y_eigen, sb_ptr, start_phi);

	return List::create(
		Named("coefficients") = fit.b,
		Named("phi") = fit.dispersion,
		Named("converged") = fit.converged
	);
}

// [[Rcpp::export]]
List fast_beta_regression_with_var_cpp(const Eigen::MatrixXd& X,
									 const NumericVector& y,
									 Nullable<NumericVector> start_beta = R_NilValue,
									 double start_phi = 10.0,
									 bool compute_std_errs = true) {

	Eigen::VectorXd y_eigen = as<Eigen::VectorXd>(y);
    Eigen::VectorXd sb;
    Eigen::VectorXd* sb_ptr = nullptr;
    if (start_beta.isNotNull()) {
        sb = as<Eigen::VectorXd>(start_beta);
        sb_ptr = &sb;
    }

    ModelResult fit = fast_beta_regression_internal(X, y_eigen, sb_ptr, start_phi);
	Eigen::MatrixXd cov_mat = fit.XtWX.inverse();

	return List::create(
		Named("coefficients") = fit.b,
		Named("phi") = fit.dispersion,
		Named("vcov") = cov_mat,
		Named("std_errs") = cov_mat.diagonal().array().sqrt(),
        Named("converged") = fit.converged
	);
}
