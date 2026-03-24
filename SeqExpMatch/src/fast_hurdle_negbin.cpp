#include "_helper_functions.h"
#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace Eigen;
using namespace LBFGSpp;

// Forward declaration from fast_logistic_regression.cpp
ModelResult fast_logistic_regression_internal(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8);

namespace {

double clamp_exp_arg_hnb(double eta) {
	return std::min(eta, 700.0);
}

class TruncatedNegBinCount {
private:
	const MatrixXd m_X;
	const VectorXi m_y;
	const int m_n;
	const int m_p;

public:
	TruncatedNegBinCount(const MatrixXd& X, const VectorXi& y)
		: m_X(X), m_y(y), m_n(X.rows()), m_p(X.cols()) {}

	double operator()(const VectorXd& params, VectorXd& grad) {
		const VectorXd beta = params.head(m_p);
		const double log_theta = params[m_p];
		const double theta = std::exp(log_theta);

		VectorXd eta = m_X * beta;
		for (int i = 0; i < m_n; ++i) {
			eta[i] = clamp_exp_arg_hnb(eta[i]);
		}
		VectorXd mu = eta.array().exp().matrix();
		mu = mu.array().max(1e-10);

		double neg_ll = 0.0;
		VectorXd score_beta = VectorXd::Zero(m_p);
		double score_log_theta = 0.0;

		for (int i = 0; i < m_n; ++i) {
			const double mu_i = mu[i];
			const double r = theta;
			const double yi = static_cast<double>(m_y[i]);

			double log_p0 = r * (std::log(r) - std::log(r + mu_i));
			double p0 = std::exp(log_p0);
			p0 = std::min(std::max(p0, 1e-12), 1.0 - 1e-12);
			const double trunc_denom = 1.0 - p0;

			neg_ll -= R::dnbinom_mu(yi, r, mu_i, true) - std::log(trunc_denom);

			const double standard_eta_score = yi - mu_i * (yi + r) / (r + mu_i);
			const double trunc_eta_corr = -mu_i * p0 * r / (trunc_denom * (r + mu_i));
			const double eta_score = standard_eta_score + trunc_eta_corr;
			score_beta.noalias() += eta_score * m_X.row(i).transpose();

			const double dlogf_dr =
				R::digamma(yi + r) - R::digamma(r) +
				std::log(r) - std::log(r + mu_i) +
				1.0 - (yi + r) / (r + mu_i);
			const double dlogp0_dr =
				std::log(r) - std::log(r + mu_i) +
				1.0 - r / (r + mu_i);
			score_log_theta += r * (dlogf_dr + (p0 * dlogp0_dr) / trunc_denom);
		}

		grad.resize(m_p + 1);
		grad.head(m_p) = -score_beta;
		grad[m_p] = -score_log_theta;

		return neg_ll;
	}

	MatrixXd hessian(const VectorXd& params) {
		MatrixXd H(m_p + 1, m_p + 1);
		H.setZero();

		const double h = 1e-5;
		VectorXd grad_at_params(m_p + 1);
		operator()(params, grad_at_params);

		for (int i = 0; i < m_p + 1; ++i) {
			VectorXd params_plus_h = params;
			params_plus_h[i] += h;
			VectorXd grad_plus_h(m_p + 1);
			operator()(params_plus_h, grad_plus_h);
			H.col(i) = (grad_plus_h - grad_at_params) / h;
		}

		H = (H + H.transpose()) / 2.0;
		return H;
	}
};

}

// [[Rcpp::export]]
List fast_hurdle_negbin_with_var_cpp(const Eigen::MatrixXd& Xmm,
									 const Eigen::VectorXd& y,
									 int j = 2,
									 int maxit = 1000,
									 double tol = 1e-8) {
	const int n = Xmm.rows();
	const int p = Xmm.cols();

	VectorXd y_pos_ind = (y.array() > 0.0).cast<double>();
	Eigen::VectorXd hurdle_b = Eigen::VectorXd::Constant(p, NA_REAL);
	double hurdle_ssq_b_j = NA_REAL;
	double hurdle_ssq_b_2 = NA_REAL;
	bool hurdle_converged = false;
    
	if (y_pos_ind.minCoeff() < y_pos_ind.maxCoeff()) {
        ModelResult hurdle_res = fast_logistic_regression_internal(Xmm, y_pos_ind);
		hurdle_b = hurdle_res.b;
		hurdle_ssq_b_j = compute_diagonal_inverse_entry(hurdle_res.XtWX, j);
		if (p >= 2) hurdle_ssq_b_2 = compute_diagonal_inverse_entry(hurdle_res.XtWX, 2);
		hurdle_converged = hurdle_res.converged;
	}

	std::vector<int> pos_rows;
	for (int i = 0; i < n; ++i) {
		if (y[i] > 0.0) pos_rows.push_back(i);
	}

	if (static_cast<int>(pos_rows.size()) <= p) {
		return List::create(
			Named("b") = NumericVector(p, NA_REAL),
			Named("ssq_b_j") = NA_REAL,
			Named("ssq_b_2") = NA_REAL,
			Named("theta_hat") = NA_REAL,
			Named("converged") = false,
			Named("hurdle_b") = hurdle_b,
			Named("hurdle_ssq_b_j") = hurdle_ssq_b_j,
			Named("hurdle_ssq_b_2") = hurdle_ssq_b_2,
			Named("hurdle_converged") = hurdle_converged
		);
	}

	MatrixXd X_pos(pos_rows.size(), p);
	VectorXi y_pos(pos_rows.size());
	for (size_t k = 0; k < pos_rows.size(); ++k) {
		const int i = pos_rows[k];
		X_pos.row(k) = Xmm.row(i);
		y_pos[k] = static_cast<int>(y[i]);
	}

	VectorXd params = VectorXd::Zero(p + 1);
	params[0] = std::log(std::max(1.0, y_pos.cast<double>().mean()));

	TruncatedNegBinCount fun(X_pos, y_pos);
	LBFGSParam<double> lbfgs_params;
	lbfgs_params.epsilon = tol;
	lbfgs_params.max_iterations = maxit;

	LBFGSSolver<double> solver(lbfgs_params);
	double neg_ll = NA_REAL;
	bool converged = false;
	try {
		int niter = solver.minimize(fun, params, neg_ll);
        converged = (niter < maxit);
	} catch (const std::exception& e) {
		Rcpp::warning(e.what());
	}

	VectorXd beta = params.head(p);
	double theta_hat = std::exp(params[p]);
	MatrixXd H = fun.hessian(params);
	double ssq_b_j = NA_REAL;
	double ssq_b_2 = NA_REAL;
	if (H.allFinite()) {
		ssq_b_j = compute_diagonal_inverse_entry(H, j);
		if (p >= 2) ssq_b_2 = compute_diagonal_inverse_entry(H, 2);
	}

	return List::create(
		Named("b") = beta,
		Named("ssq_b_j") = ssq_b_j,
		Named("ssq_b_2") = ssq_b_2,
		Named("theta_hat") = theta_hat,
		Named("converged") = converged,
		Named("hurdle_b") = hurdle_b,
		Named("hurdle_ssq_b_j") = hurdle_ssq_b_j,
		Named("hurdle_ssq_b_2") = hurdle_ssq_b_2,
		Named("hurdle_converged") = hurdle_converged
	);
}
