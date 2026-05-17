#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;

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
	ZeroOneInflatedBeta(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& X_zero_one):
		m_y(y),
		m_X(X),
		m_X_zero_one(X_zero_one),
		m_n(X.rows()),
		m_p(X.cols()),
		m_p_zero_one(X_zero_one.cols())
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
			m_X_beta.row(j) = m_X.row(row);
			m_y_beta[j] = m_y[row];
		}
		m_log_y_beta = m_y_beta.array().log().matrix();
		m_log1_y_beta = (1.0 - m_y_beta.array()).log().matrix();
	}

	double operator()(const Eigen::VectorXd& params, Eigen::VectorXd& grad){
		Eigen::VectorXd beta = params.head(m_p);
		double log_phi = params[m_p];
		double phi = std::exp(log_phi);

		Eigen::VectorXd gamma0 = params.segment(m_p + 1, m_p_zero_one);
		Eigen::VectorXd gamma1 = params.tail(m_p_zero_one);

		Eigen::VectorXd eta = m_X * beta;
		Eigen::VectorXd mu = logistic(eta);
		clamp_probs(mu);

		Eigen::VectorXd eta_beta = m_X_beta * beta;
		Eigen::VectorXd mu_beta = logistic(eta_beta);
		clamp_probs(mu_beta);

		Eigen::VectorXd eta0 = m_X_zero_one * gamma0;
		Eigen::VectorXd eta1 = m_X_zero_one * gamma1;
		Eigen::VectorXd pi0(m_n);
		Eigen::VectorXd pi1(m_n);
		Eigen::VectorXd pib(m_n);
		double mixture_loglik = 0.0;
		Eigen::VectorXd grad_gamma0 = Eigen::VectorXd::Zero(m_p_zero_one);
		Eigen::VectorXd grad_gamma1 = Eigen::VectorXd::Zero(m_p_zero_one);
		for (int i = 0; i < m_n; ++i){
			double max_logit = std::max(std::max(eta0[i], eta1[i]), 0.0);
			double e0 = std::exp(eta0[i] - max_logit);
			double e1 = std::exp(eta1[i] - max_logit);
			double eb = std::exp(-max_logit);
			double denom = e0 + e1 + eb;
			pi0[i] = e0 / denom;
			pi1[i] = e1 / denom;
			pib[i] = eb / denom;
			const Eigen::VectorXd xzi = m_X_zero_one.row(i).transpose();
			if (m_y[i] <= 0){
				mixture_loglik += std::log(pi0[i]);
				grad_gamma0.noalias() += xzi * (1.0 - pi0[i]);
				grad_gamma1.noalias() -= xzi * pi1[i];
			} else if (m_y[i] >= 1){
				mixture_loglik += std::log(pi1[i]);
				grad_gamma0.noalias() -= xzi * pi0[i];
				grad_gamma1.noalias() += xzi * (1.0 - pi1[i]);
			} else {
				mixture_loglik += std::log(pib[i]);
				grad_gamma0.noalias() -= xzi * pi0[i];
				grad_gamma1.noalias() -= xzi * pi1[i];
			}
		}

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

		grad.resize(m_p + 1 + 2 * m_p_zero_one);

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
		grad.segment(m_p + 1, m_p_zero_one) = -grad_gamma0;
		grad.tail(m_p_zero_one) = -grad_gamma1;

		return neg_ll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& params){
		Eigen::MatrixXd H(m_p + 1 + 2 * m_p_zero_one, m_p + 1 + 2 * m_p_zero_one);
		H.setZero();

		Eigen::VectorXd beta = params.head(m_p);
		double log_phi = params[m_p];
		double phi = std::exp(log_phi);

		if (m_n_beta > 0){
			Eigen::VectorXd eta_beta = m_X_beta * beta;
			Eigen::VectorXd mu_beta = logistic(eta_beta);
			clamp_probs(mu_beta);

			double d_neg_ll_d_phi = -m_n_beta * R::digamma(phi);
			double d2_neg_ll_d_phi2 = -m_n_beta * R::trigamma(phi);

			for (int i = 0; i < m_n_beta; ++i){
				const double mu = mu_beta[i];
				const double one_minus_mu = 1.0 - mu;
				const double dmu_deta = mu * one_minus_mu;
				const double d2mu_deta2 = dmu_deta * (1.0 - 2.0 * mu);
				const double a = mu * phi;
				const double b = one_minus_mu * phi;
				const double digamma_a = R::digamma(a);
				const double digamma_b = R::digamma(b);
				const double trigamma_a = R::trigamma(a);
				const double trigamma_b = R::trigamma(b);
				const double c = digamma_a - digamma_b - m_log_y_beta[i] + m_log1_y_beta[i];
				const double dc_dmu = phi * (trigamma_a + trigamma_b);
				const double dc_dphi = mu * trigamma_a - one_minus_mu * trigamma_b;
				const double dscore_eta_deta =
					phi * (dc_dmu * dmu_deta * dmu_deta + c * d2mu_deta2);
				const double dscore_eta_dlogphi =
					phi * dmu_deta * (c + phi * dc_dphi);
				const Eigen::VectorXd x = m_X_beta.row(i).transpose();

				H.topLeftCorner(m_p, m_p).noalias() += dscore_eta_deta * (x * x.transpose());
				H.topRightCorner(m_p, 1).noalias() += dscore_eta_dlogphi * x;

				d_neg_ll_d_phi +=
					mu * digamma_a +
					one_minus_mu * digamma_b -
					mu * m_log_y_beta[i] -
					one_minus_mu * m_log1_y_beta[i];
				d2_neg_ll_d_phi2 +=
					mu * mu * trigamma_a +
					one_minus_mu * one_minus_mu * trigamma_b;
			}

			H(m_p, m_p) = phi * d_neg_ll_d_phi + phi * phi * d2_neg_ll_d_phi2;
		}

		Eigen::VectorXd gamma0 = params.segment(m_p + 1, m_p_zero_one);
		Eigen::VectorXd gamma1 = params.tail(m_p_zero_one);
		Eigen::MatrixXd H00 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::MatrixXd H11 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::MatrixXd H01 = Eigen::MatrixXd::Zero(m_p_zero_one, m_p_zero_one);
		Eigen::VectorXd eta0 = m_X_zero_one * gamma0;
		Eigen::VectorXd eta1 = m_X_zero_one * gamma1;
		for (int i = 0; i < m_n; ++i){
			double max_logit = std::max(std::max(eta0[i], eta1[i]), 0.0);
			double e0 = std::exp(eta0[i] - max_logit);
			double e1 = std::exp(eta1[i] - max_logit);
			double eb = std::exp(-max_logit);
			double denom = e0 + e1 + eb;
			double pi0 = e0 / denom;
			double pi1 = e1 / denom;
			const Eigen::VectorXd xzi = m_X_zero_one.row(i).transpose();
			H00.noalias() += pi0 * (1.0 - pi0) * (xzi * xzi.transpose());
			H11.noalias() += pi1 * (1.0 - pi1) * (xzi * xzi.transpose());
			H01.noalias() += (-pi0 * pi1) * (xzi * xzi.transpose());
		}
		const int g0_start = m_p + 1;
		const int g1_start = m_p + 1 + m_p_zero_one;
		H.block(g0_start, g0_start, m_p_zero_one, m_p_zero_one) = H00;
		H.block(g1_start, g1_start, m_p_zero_one, m_p_zero_one) = H11;
		H.block(g0_start, g1_start, m_p_zero_one, m_p_zero_one) = H01;
		H.block(g1_start, g0_start, m_p_zero_one, m_p_zero_one) = H01.transpose();
		return H;
	}

private:
	const Eigen::VectorXd m_y;
	const Eigen::MatrixXd m_X;
	const Eigen::MatrixXd m_X_zero_one;
	Eigen::MatrixXd m_X_beta;
	Eigen::VectorXd m_y_beta;
	Eigen::VectorXd m_log_y_beta;
	Eigen::VectorXd m_log1_y_beta;
	int m_n;
	int m_p;
	int m_p_zero_one;
	int m_n_beta;
	int m_n_zero;
	int m_n_one;
};

// [[Rcpp::export]]
Eigen::VectorXd get_zero_one_inflated_beta_score_cpp(Eigen::MatrixXd X,
													 Eigen::MatrixXd X_zero_one,
													 NumericVector y,
													 NumericVector params){
	Eigen::VectorXd y_eigen = Rcpp::as<Eigen::VectorXd>(y);
	Eigen::VectorXd par = Rcpp::as<Eigen::VectorXd>(params);
	ZeroOneInflatedBeta fun(y_eigen, X, X_zero_one);
	Eigen::VectorXd grad(par.size());
	fun(par, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_zero_one_inflated_beta_hessian_cpp(Eigen::MatrixXd X,
													   Eigen::MatrixXd X_zero_one,
													   NumericVector y,
													   NumericVector params){
	Eigen::VectorXd y_eigen = Rcpp::as<Eigen::VectorXd>(y);
	Eigen::VectorXd par = Rcpp::as<Eigen::VectorXd>(params);
	ZeroOneInflatedBeta fun(y_eigen, X, X_zero_one);
	return -fun.hessian(par);
}

//' @title Fast Zero/One-Inflated Beta Regression (C++)
//' @description High-performance zero/one-inflated beta regression fitting using Newton-Raphson or L-BFGS.
//' @param X Matrix of predictors for the beta component.
//' @param X_zero_one Matrix of predictors for the zero and one inflation components.
//' @param y Vector of responses in [0, 1].
//' @param warm_start_params Optional starting values for all parameters. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first iteration.
//' @return A list containing coefficients, vcov, and convergence status.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_zero_one_inflated_beta_cpp(Eigen::MatrixXd X,
									 Eigen::MatrixXd X_zero_one,
									 NumericVector y,
									 Nullable<NumericVector> warm_start_params = R_NilValue,
									 bool smart_cold_start = true,
									 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
									 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
									 std::string optimization_alg = "lbfgs",
                                     Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue){
	Eigen::VectorXd y_eigen = Rcpp::as<Eigen::VectorXd>(y);
	int p = X.cols();
	int p_zo = X_zero_one.cols();
	int total = p + 1 + 2 * p_zo;
	Eigen::VectorXd params(total);
	
	if (warm_start_params.isNotNull()) {
		params = Rcpp::as<Eigen::VectorXd>(warm_start_params);
	} else if (smart_cold_start) {
		// Beta component: OLS on logit(y) for entries in (0, 1)
		std::vector<int> idx_mid;
		for(int i=0; i<y_eigen.size(); ++i) if(y_eigen[i] > 0 && y_eigen[i] < 1) idx_mid.push_back(i);
		if (idx_mid.size() > (size_t)p) {
			Eigen::MatrixXd X_mid(idx_mid.size(), p);
			Eigen::VectorXd y_logit(idx_mid.size());
			for(size_t i=0; i<idx_mid.size(); ++i) {
				X_mid.row(i) = X.row(idx_mid[i]);
				double yi = y_eigen[idx_mid[i]];
				y_logit[i] = std::log(yi / (1.0 - yi));
			}
			params.head(p) = safe_ols_solve(X_mid, y_logit);
		} else {
			params.head(p).setZero();
		}
		params[p] = 2.0; // log_phi warm_start_params
		
		// Zero/One components: OLS on indicators
		Eigen::VectorXd y_is_zero = (y_eigen.array() == 0.0).cast<double>();
		Eigen::VectorXd y_is_one  = (y_eigen.array() == 1.0).cast<double>();
		params.segment(p + 1, p_zo) = ols_smart_cold_start_beta(X_zero_one, y_is_zero);
		params.tail(p_zo)           = ols_smart_cold_start_beta(X_zero_one, y_is_one);
	} else {
		params.setZero();
		params[p] = 2.0;
	}
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	ZeroOneInflatedBeta fun(y_eigen, X, X_zero_one);
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }
    
	LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, fixed_spec, 1500, 1e-6, optimization_alg, "lbfgs", 0, info_start_ptr);
	params = fit.params;

	Eigen::MatrixXd H = fun.hessian(params);
	int dim = params.size();
	NumericMatrix vcov_mat(dim, dim);
	bool has_vcov = false;
	if (H.rows() == dim && H.cols() == dim && H.allFinite()){
		Eigen::MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::FullPivLU<Eigen::MatrixXd> lu(H_free);
		if (lu.isInvertible()){
			Eigen::MatrixXd inv_free = lu.inverse();
			Eigen::MatrixXd inv = expand_free_covariance(dim, fixed_spec, inv_free, true);
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

	int p_zero_one = X_zero_one.cols();
	return List::create(
		Named("b") = params.head(p),
		Named("log_phi") = params[p],
		Named("zero_one_b0") = params.segment(p + 1, p_zero_one),
		Named("zero_one_b1") = params.tail(p_zero_one),
		Named("params") = params,
		Named("vcov") = vcov_mat,
		Named("neg_loglik") = fit.value,
		Named("fisher_information") = H
	);
}
