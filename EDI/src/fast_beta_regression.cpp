#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

namespace {

struct DigammaFunctor {
	double operator()(double x) const {
		return R::digamma(x);
	}
};

struct TrigammaFunctor {
	double operator()(double x) const {
		return R::trigamma(x);
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
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total_p, total_p);
		Eigen::VectorXd beta = params.head(m_p);
		double phi = std::exp(params[m_p]);
		Eigen::VectorXd eta = m_X * beta;
		Eigen::VectorXd mu = (1.0 / (1.0 + (-eta).array().exp())).matrix();
		double epsilon = 1e-8;
		for (int i = 0; i < m_n; ++i) {
			if (mu[i] < epsilon) mu[i] = epsilon;
			if (mu[i] > 1.0 - epsilon) mu[i] = 1.0 - epsilon;
		}

		Eigen::VectorXd a = mu.array() * phi;
		Eigen::VectorXd b = (1.0 - mu.array()) * phi;
		Eigen::VectorXd dig_a = a.unaryExpr(DigammaFunctor());
		Eigen::VectorXd dig_b = b.unaryExpr(DigammaFunctor());
		Eigen::VectorXd tri_a = a.unaryExpr(TrigammaFunctor());
		Eigen::VectorXd tri_b = b.unaryExpr(TrigammaFunctor());

		for (int i = 0; i < m_n; ++i) {
			double mui = mu[i];
			double dmu = mui * (1.0 - mui);
			double d2mu = dmu * (1.0 - 2.0 * mui);
			double C = dig_a[i] - dig_b[i] - m_log_y[i] + m_log1_y[i];
			double B = phi * C;
			double B_mu = phi * phi * (tri_a[i] + tri_b[i]);
			double w_beta = B_mu * dmu * dmu + B * d2mu;
			Eigen::VectorXd x = m_X.row(i).transpose();

			H.topLeftCorner(m_p, m_p).noalias() += w_beta * (x * x.transpose());

			double B_log_phi = phi * (C + a[i] * tri_a[i] - b[i] * tri_b[i]);
			H.topRightCorner(m_p, 1).noalias() += B_log_phi * dmu * x;
		}

		double D = -m_n * R::digamma(phi);
		double D_phi = -m_n * R::trigamma(phi);
		for (int i = 0; i < m_n; ++i) {
			double mui = mu[i];
			D += mui * dig_a[i] + (1.0 - mui) * dig_b[i] -
				mui * m_log_y[i] - (1.0 - mui) * m_log1_y[i];
			D_phi += mui * mui * tri_a[i] + (1.0 - mui) * (1.0 - mui) * tri_b[i];
		}
		H(m_p, m_p) = phi * D + phi * phi * D_phi;
		H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();
		return H;
	}
};

ModelResult fast_beta_regression_internal(const Eigen::MatrixXd& X,
                                        const Eigen::VectorXd& y,
                                        const Eigen::VectorXd* start_beta = nullptr,
                                        double start_phi = 10.0,
                                        Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                        Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                        std::string optimization_alg = "newton_raphson",
                                        Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
    int p = X.cols();
    ModelResult res;
    Eigen::VectorXd params(p + 1);
    if (start_beta) params.head(p) = *start_beta;
    else params.head(p).setZero();
    params[p] = std::log(start_phi);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);

    BetaRegression fun(y, X);
    
    Eigen::MatrixXd H_start;
    const Eigen::MatrixXd* h_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        H_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        h_ptr = &H_start;
    }
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, 1000, 1e-6, optimization_alg, "newton_raphson", 0, h_ptr);
    params = fit.params;

    res.b = params.head(p);
    res.dispersion = std::exp(params[p]); // phi
    res.XtWX = fun.hessian(params); // This is actually the Hessian H
    res.converged = fit.converged;
    return res;
}

} // namespace

//' @title Compute Beta Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a beta regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (in (0, 1)).
//' @param params A numeric vector of parameters [beta, log_phi].
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::VectorXd get_beta_regression_score_cpp(const Eigen::MatrixXd& X,
                                              const Eigen::VectorXd& y,
                                              const Eigen::VectorXd& params) {
    BetaRegression fun(y, X);
    Eigen::VectorXd grad(params.size());
    fun(params, grad);
    return -grad; // Return the actual score (gradient of log-likelihood)
}

//' @title Compute Beta Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a beta regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses.
//' @param params A numeric vector of parameters [beta, log_phi].
//' @return A numeric matrix representing the Hessian.
//' @export
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd get_beta_regression_hessian_cpp(const Eigen::MatrixXd& X,
                                                const Eigen::VectorXd& y,
                                                const Eigen::VectorXd& params) {
    BetaRegression fun(y, X);
    return -fun.hessian(params); // Return the actual Hessian of log-likelihood (Fisher Information is -Hessian)
}

//' @title Fast Beta Regression (C++)
//' @description High-performance beta regression fitting using Newton-Raphson or L-BFGS.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (in (0, 1)).
//' @param start_beta Optional starting values for coefficients.
//' @param start_phi Optional starting value for precision parameter phi.
//' @param compute_std_errs Deprecated.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, phi, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = runif(10)
//' fast_beta_regression_cpp(X, y)
// [[Rcpp::export]]
List fast_beta_regression_cpp(const Eigen::MatrixXd& X,
								const NumericVector& y,
								Nullable<NumericVector> start_beta = R_NilValue,
								double start_phi = 10.0,
                                bool compute_std_errs = false,
                                Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                std::string optimization_alg = "newton_raphson",
                                Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {

	Eigen::VectorXd y_eigen = as<Eigen::VectorXd>(y);
    Eigen::VectorXd sb;
    Eigen::VectorXd* sb_ptr = nullptr;
    if (start_beta.isNotNull()) {
        sb = as<Eigen::VectorXd>(start_beta);
        sb_ptr = &sb;
    }

    ModelResult fit = fast_beta_regression_internal(X, y_eigen, sb_ptr, start_phi, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);

    Eigen::VectorXd params_full(fit.b.size() + 1);
    params_full.head(fit.b.size()) = fit.b;
    params_full[fit.b.size()] = std::log(fit.dispersion);
    BetaRegression fun_neg_ll(y_eigen, X);
    Eigen::VectorXd dummy_grad(params_full.size());
    double neg_loglik = fun_neg_ll(params_full, dummy_grad);

	return List::create(
		Named("coefficients") = fit.b,
		Named("phi") = fit.dispersion,
		Named("neg_loglik") = neg_loglik,
		Named("converged") = fit.converged
	);
}

//' @title Fast Beta Regression with Variance (C++)
//' @description Beta regression with full variance-covariance matrix and standard error estimation.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (in (0, 1)).
//' @param start_beta Optional starting values for coefficients.
//' @param start_phi Optional starting value for precision parameter phi.
//' @param compute_std_errs Deprecated.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, phi, vcov, standard errors, and convergence status.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = runif(10)
//' fast_beta_regression_with_var_cpp(X, y)
// [[Rcpp::export]]
List fast_beta_regression_with_var_cpp(const Eigen::MatrixXd& X,
									 const NumericVector& y,
									 Nullable<NumericVector> start_beta = R_NilValue,
									 double start_phi = 10.0,
                                     bool compute_std_errs = true,
                                     Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "newton_raphson") {

	Eigen::VectorXd y_eigen = as<Eigen::VectorXd>(y);
    Eigen::VectorXd sb;
    Eigen::VectorXd* sb_ptr = nullptr;
    if (start_beta.isNotNull()) {
        sb = as<Eigen::VectorXd>(start_beta);
        sb_ptr = &sb;
    }

    ModelResult fit = fast_beta_regression_internal(X, y_eigen, sb_ptr, start_phi, fixed_idx, fixed_values, optimization_alg);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols() + 1, fixed_idx, fixed_values);
    Eigen::MatrixXd H_free = subset_matrix(fit.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
	Eigen::MatrixXd cov_free = H_free.inverse();
    Eigen::MatrixXd cov_mat = expand_free_covariance(X.cols() + 1, fixed_spec, cov_free, true);
    Eigen::VectorXd se = cov_mat.diagonal().array().sqrt();

    Eigen::VectorXd params_full(fit.b.size() + 1);
    params_full.head(fit.b.size()) = fit.b;
    params_full[fit.b.size()] = std::log(fit.dispersion);
    BetaRegression fun_neg_ll(y_eigen, X);
    Eigen::VectorXd dummy_grad(params_full.size());
    double neg_loglik = fun_neg_ll(params_full, dummy_grad);

	return List::create(
		Named("coefficients") = fit.b,
		Named("phi") = fit.dispersion,
		Named("neg_loglik") = neg_loglik,
		Named("vcov") = cov_mat,
		Named("std_errs") = se,
        Named("converged") = fit.converged
	);
}
