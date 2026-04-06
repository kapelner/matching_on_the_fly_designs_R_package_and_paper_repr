#include <RcppEigen.h>
#include <optimization/LBFGS.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace LBFGSpp;

struct DigammaFunctor {
	double operator()(double x) const {
		return R::digamma(x);
	}
};

struct LgammaFunctor {
	double operator()(double x) const {
		return R::lgammafn(x);
	}
};

inline Eigen::VectorXd logistic(const Eigen::VectorXd& x){
	return (1.0 / (1.0 + (-x).array().exp())).matrix();
}

inline void clamp_probs(Eigen::VectorXd& v, double lower = 1e-8, double upper = 1.0 - 1e-8){
	for (int i = 0; i < v.size(); ++i){
		if (v[i] < lower){
			v[i] = lower;
		} else if (v[i] > upper){
			v[i] = upper;
		}
	}
}

class ZeroOneInflatedBeta {
public:
	ZeroOneInflatedBeta(const Eigen::VectorXd& y, const Eigen::MatrixXd& Xfull):
		m_y(y),
		m_Xfull(Xfull),
		m_n(Xfull.rows()),
		m_p(Xfull.cols())
	{
		m_n_zero = 0;
		m_n_one = 0;
		std::vector<int> beta_idx;
		for (int i = 0; i < m_n; ++i){
			if (m_y[i] <= 0){
				++m_n_zero;
			} else if (m_y[i] >= 1){
				++m_n_one;
			} else {
				beta_idx.push_back(i);
			}
		}

		m_n_beta = beta_idx.size();
		m_X_beta.resize(m_n_beta, m_p);
		m_y_beta.resize(m_n_beta);
		for (int j = 0; j < m_n_beta; ++j){
			int row = beta_idx[j];
			m_X_beta.row(j) = m_Xfull.row(row);
			m_y_beta[j] = m_y[row];
		}
		m_log_y_beta = m_y_beta.array().log().matrix();
		m_log1_y_beta = (1.0 - m_y_beta.array()).log().matrix();
	}

	double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad){
		Eigen::VectorXd beta = params.head(m_p);
		double log_phi = params[m_p];
		double phi = std::exp(log_phi);

		double alpha0 = params[m_p + 1];
		double alpha1 = params[m_p + 2];

		Eigen::VectorXd eta = m_Xfull * beta;
		Eigen::VectorXd mu = logistic(eta);
		clamp_probs(mu);

		Eigen::VectorXd eta_beta = m_X_beta * beta;
		Eigen::VectorXd mu_beta = logistic(eta_beta);
		clamp_probs(mu_beta);

		Eigen::Vector3d logits;
		logits << alpha0, alpha1, 0.0;
		double max_logit = logits.maxCoeff();
		Eigen::Vector3d exp_logits = (logits.array() - max_logit).exp();
		double sum_exp = exp_logits.sum();
		Eigen::Vector3d pis = exp_logits / sum_exp;
		double pi0 = pis[0];
		double pi1 = pis[1];
		double pib = pis[2];

		double mixture_loglik = m_n_zero * std::log(pi0) +
			m_n_one * std::log(pi1) +
			m_n_beta * std::log(pib);

		Eigen::VectorXd mu_beta_phi = mu_beta.array() * phi;
		Eigen::VectorXd one_minus_mu_beta_phi = (1.0 - mu_beta.array()) * phi;
		double sum_loggamma_mu_phi = mu_beta_phi.unaryExpr(LgammaFunctor()).sum();
		double sum_loggamma_one_minus_mu_phi = one_minus_mu_beta_phi.unaryExpr(LgammaFunctor()).sum();
		double sum_term3 = ((mu_beta_phi.array() - 1.0) * m_log_y_beta.array()).sum();
		double sum_term4 = ((one_minus_mu_beta_phi.array() - 1.0) * m_log1_y_beta.array()).sum();
		double neg_ll_beta = -(
			m_n_beta * R::lgammafn(phi) -
			sum_loggamma_mu_phi -
			sum_loggamma_one_minus_mu_phi +
			sum_term3 +
			sum_term4
		);

		double neg_ll = -mixture_loglik + neg_ll_beta;

		grad.resize(m_p + 3);

		Eigen::VectorXd d_mu_d_eta_beta = mu_beta.array() * (1.0 - mu_beta.array());
		Eigen::VectorXd digamma_mu_phi = mu_beta_phi.unaryExpr(DigammaFunctor());
		Eigen::VectorXd digamma_one_minus = one_minus_mu_beta_phi.unaryExpr(DigammaFunctor());

		Eigen::VectorXd d_neg_ll_d_mu_beta = (
			(digamma_mu_phi - digamma_one_minus - m_log_y_beta + m_log1_y_beta).array()
		) * phi;

		grad.head(m_p) = m_X_beta.transpose() *
			(d_neg_ll_d_mu_beta.array() * d_mu_d_eta_beta.array()).matrix();

		double sum_mu_digamma = (mu_beta.array() * digamma_mu_phi.array()).sum();
		double sum_one_minus_digamma = ((1.0 - mu_beta.array()) * digamma_one_minus.array()).sum();
		double sum_mu_logy = (mu_beta.array() * m_log_y_beta.array()).sum();
		double sum_one_minus_log1y = ((1.0 - mu_beta.array()) * m_log1_y_beta.array()).sum();

		double d_neg_ll_d_phi = -m_n_beta * R::digamma(phi) +
			sum_mu_digamma +
			sum_one_minus_digamma -
			sum_mu_logy -
			sum_one_minus_log1y;
		grad[m_p] = d_neg_ll_d_phi * phi;

		grad[m_p + 1] = pi0 * m_n - m_n_zero;
		grad[m_p + 2] = pi1 * m_n - m_n_one;

		return neg_ll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& params){
		Eigen::MatrixXd H(m_p + 3, m_p + 3);
		H.setZero();

		Eigen::VectorXd grad_at_params(m_p + 3);
		operator()(params, grad_at_params);

		double h = 1e-6;
		for (int i = 0; i < m_p + 3; ++i){
			Eigen::VectorXd params_plus = params;
			params_plus[i] += h;
			Eigen::VectorXd grad_plus(m_p + 3);
			operator()(params_plus, grad_plus);
			H.col(i) = (grad_plus - grad_at_params) / h;
		}
		H = (H + H.transpose()) * 0.5;
		return H;
	}

private:
	const Eigen::VectorXd m_y;
	const Eigen::MatrixXd m_Xfull;
	Eigen::MatrixXd m_X_beta;
	Eigen::VectorXd m_y_beta;
	Eigen::VectorXd m_log_y_beta;
	Eigen::VectorXd m_log1_y_beta;
	int m_n;
	int m_p;
	int m_n_beta;
	int m_n_zero;
	int m_n_one;
};

// [[Rcpp::export]]
List fast_zero_one_inflated_beta_cpp(Eigen::MatrixXd Xfull, NumericVector y, NumericVector init){
	Eigen::VectorXd y_eigen = Rcpp::as<Eigen::VectorXd>(y);
	Eigen::VectorXd params = Rcpp::as<Eigen::VectorXd>(init);

	ZeroOneInflatedBeta fun(y_eigen, Xfull);
	LBFGSParam<double> lbfgs_params;
	lbfgs_params.epsilon = 1e-6;
	lbfgs_params.max_iterations = 1500;
	LBFGSSolver<double> solver(lbfgs_params);

	double neg_ll;
	solver.minimize(fun, params, neg_ll);

	Eigen::MatrixXd H = fun.hessian(params);
	int dim = params.size();
	NumericMatrix vcov_mat(dim, dim);
	bool has_vcov = false;
	if (H.rows() == dim && H.cols() == dim && H.allFinite()){
		Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
		if (lu.isInvertible()){
			Eigen::MatrixXd inv = lu.inverse();
			for (int i = 0; i < dim; ++i){
				for (int j = 0; j < dim; ++j){
					vcov_mat(i, j) = inv(i, j);
				}
			}
			has_vcov = true;
		}
	}
	if (!has_vcov){
		for (int i = 0; i < dim; ++i){
			for (int j = 0; j < dim; ++j){
				vcov_mat(i, j) = NA_REAL;
			}
		}
	}

	NumericVector coeff = wrap(params);
	return List::create(
		Named("coefficients") = coeff,
		Named("vcov") = vcov_mat,
		Named("neg_loglik") = neg_ll
	);
}
