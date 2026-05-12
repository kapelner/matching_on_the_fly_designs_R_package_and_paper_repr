#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace Eigen;

// Forward declaration from fast_logistic_regression.cpp
ModelResult fast_logistic_regression_internal(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& weights = Eigen::VectorXd(),
	Rcpp::Nullable<Rcpp::NumericVector> start_beta = R_NilValue,
	bool smart_start = true,
	int maxit = 100,
	double tol = 1e-8,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs");

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

		VectorXd eta = (m_X * beta).array().min(700.0).matrix();
		VectorXd mu = eta.array().exp().max(1e-10).matrix();

		double neg_ll = 0.0;
		VectorXd eta_score = VectorXd::Zero(m_n);
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
			eta_score[i] = standard_eta_score + trunc_eta_corr;

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
		grad.head(m_p) = -(m_X.transpose() * eta_score);
		grad[m_p] = -score_log_theta;

		return neg_ll;
	}

	MatrixXd hessian(const VectorXd& params) {
		const int total_p = m_p + 1;
		MatrixXd H = MatrixXd::Zero(total_p, total_p);
		const VectorXd beta = params.head(m_p);
		const double r = std::exp(params[m_p]);
		VectorXd eta = (m_X * beta).array().min(700.0).matrix();
		VectorXd beta_weights = VectorXd::Zero(m_n);
		VectorXd cross_weights = VectorXd::Zero(m_n);

		for (int i = 0; i < m_n; ++i) {
			const double mu_i = std::max(std::exp(eta[i]), 1e-10);
			const double yi = static_cast<double>(m_y[i]);
			const double denom = r + mu_i;
			const double denom_sq = denom * denom;

			double log_p0 = r * (std::log(r) - std::log(denom));
			double p0 = std::exp(log_p0);
			p0 = std::min(std::max(p0, 1e-12), 1.0 - 1e-12);
			const double q0 = 1.0 - p0;

			const double c_eta = r * mu_i / denom;
			const double dlogp0_dr =
				std::log(r) - std::log(denom) + 1.0 - r / denom;

			const double d_score_eta_d_eta =
				- (yi + r) * r * mu_i / denom_sq
				- r * r * mu_i * p0 / (denom_sq * q0)
				+ c_eta * c_eta * p0 / (q0 * q0);
			beta_weights[i] = -d_score_eta_d_eta;

			const double d_score_eta_d_log_r =
				r * mu_i * (yi - mu_i) / denom_sq
				- r * mu_i * mu_i * p0 / (denom_sq * q0)
				- r * r * mu_i * p0 * dlogp0_dr / (denom * q0 * q0);
			cross_weights[i] = -d_score_eta_d_log_r;

			const double dlogf_dr =
				R::digamma(yi + r) - R::digamma(r) +
				std::log(r) - std::log(denom) +
				1.0 - (yi + r) / denom;
			const double d2logf_dr2 =
				R::trigamma(yi + r) - R::trigamma(r) +
				1.0 / r - 1.0 / denom + (yi - mu_i) / denom_sq;
			const double d2logp0_dr2 =
				1.0 / r - 1.0 / denom - mu_i / denom_sq;
			const double p0_over_q0 = p0 / q0;
			const double d_score_log_r_d_log_r =
				r * (dlogf_dr + p0_over_q0 * dlogp0_dr) +
				r * r * (
					d2logf_dr2 +
					p0 * dlogp0_dr * dlogp0_dr / (q0 * q0) +
					p0_over_q0 * d2logp0_dr2
				);
			H(m_p, m_p) -= d_score_log_r_d_log_r;
		}

		H.topLeftCorner(m_p, m_p).noalias() = weighted_crossprod(m_X, beta_weights);
		H.topRightCorner(m_p, 1).noalias() = m_X.transpose() * cross_weights;
		H.bottomLeftCorner(1, m_p) = H.topRightCorner(m_p, 1).transpose();
		return H;
	}
};

}

static List build_positive_hurdle_negbin_data(const Eigen::MatrixXd& X,
											  const Eigen::VectorXd& y) {
	std::vector<int> pos_rows;
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_rows.push_back(i);
	}

	const int p = X.cols();
	MatrixXd X_pos(pos_rows.size(), p);
	VectorXi y_pos(pos_rows.size());
	for (size_t k = 0; k < pos_rows.size(); ++k) {
		const int i = pos_rows[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = static_cast<int>(y[i]);
	}
	return List::create(Named("X_pos") = X_pos, Named("y_pos") = y_pos);
}

// [[Rcpp::export]]
Eigen::VectorXd get_hurdle_negbin_count_score_cpp(const Eigen::MatrixXd& X,
												  const Eigen::VectorXd& y,
												  const Eigen::VectorXd& params) {
	List pos = build_positive_hurdle_negbin_data(X, y);
	MatrixXd X_pos = pos["X_pos"];
	VectorXi y_pos = pos["y_pos"];
	TruncatedNegBinCount fun(X_pos, y_pos);
	VectorXd grad(params.size());
	fun(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_hurdle_negbin_count_hessian_cpp(const Eigen::MatrixXd& X,
													const Eigen::VectorXd& y,
													const Eigen::VectorXd& params) {
	List pos = build_positive_hurdle_negbin_data(X, y);
	MatrixXd X_pos = pos["X_pos"];
	VectorXi y_pos = pos["y_pos"];
	TruncatedNegBinCount fun(X_pos, y_pos);
	return -fun.hessian(params);
}

// [[Rcpp::export]]
List fast_hurdle_negbin_cpp(const Eigen::MatrixXd& X,
						   const Eigen::VectorXd& y,
						   const Eigen::MatrixXd& X_hurdle,
						   int maxit = 1000,
						   double tol = 1e-8,
						   Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
						   Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
						   std::string optimization_alg = "lbfgs") {
	const int n = X.rows();
	const int p = X.cols();
	const int p_hurdle = X_hurdle.cols();
	std::string alg = normalize_optimizer_algorithm(optimization_alg, "lbfgs", false);
	FixedParamSpec count_fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);

	VectorXd y_pos_ind = (y.array() > 0.0).cast<double>();
	Eigen::VectorXd hurdle_b = Eigen::VectorXd::Constant(p_hurdle, NA_REAL);
	bool hurdle_converged = false;
    
	if (y_pos_ind.minCoeff() < y_pos_ind.maxCoeff()) {
        ModelResult hurdle_res = fast_logistic_regression_internal(X_hurdle, y_pos_ind, Eigen::VectorXd(), R_NilValue, true, 100, 1e-8, R_NilValue, R_NilValue, alg);
		hurdle_b = hurdle_res.b;
		hurdle_converged = hurdle_res.converged;
	}

	std::vector<int> pos_rows;
	for (int i = 0; i < n; ++i) {
		if (y[i] > 0.0) pos_rows.push_back(i);
	}

	if (static_cast<int>(pos_rows.size()) <= p) {
		return List::create(
			Named("b") = NumericVector(p, NA_REAL),
			Named("theta_hat") = NA_REAL,
			Named("converged") = false,
			Named("hurdle_b") = hurdle_b,
			Named("hurdle_converged") = hurdle_converged
		);
	}

	MatrixXd X_pos(pos_rows.size(), p);
	VectorXi y_pos(pos_rows.size());
	for (size_t k = 0; k < pos_rows.size(); ++k) {
		const int i = pos_rows[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = static_cast<int>(y[i]);
	}

	VectorXd params = VectorXd::Zero(p + 1);
	params[0] = std::log(std::max(1.0, y_pos.cast<double>().mean()));

	TruncatedNegBinCount fun(X_pos, y_pos);
	double neg_ll = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(fun, params, count_fixed_spec, maxit, tol, alg, "lbfgs");
		params = fit.params;
		neg_ll = fit.value;
        converged = fit.converged;
	} catch (const std::exception& e) {
		Rcpp::warning(e.what());
	}

	VectorXd beta = params.head(p);
	double theta_hat = std::exp(params[p]);

	return List::create(
		Named("b") = beta,
		Named("theta_hat") = theta_hat,
		Named("converged") = converged,
		Named("hurdle_b") = hurdle_b,
		Named("hurdle_converged") = hurdle_converged
	);
}

// [[Rcpp::export]]
List fast_hurdle_negbin_with_var_cpp(const Eigen::MatrixXd& X,
									 const Eigen::VectorXd& y,
									 const Eigen::MatrixXd& X_hurdle,
									 int j = 2,
									 int maxit = 1000,
									 double tol = 1e-8,
									 Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
									 Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
									 std::string optimization_alg = "lbfgs") {
	List fit = fast_hurdle_negbin_cpp(X, y, X_hurdle, maxit, tol, fixed_idx, fixed_values, optimization_alg);
	SEXP b_sexp = fit["b"];
	NumericVector b_nv(b_sexp);
	const int p = b_nv.size();

	double hurdle_ssq_b_j = NA_REAL;
	double hurdle_ssq_b_2 = NA_REAL;
	if (fit.containsElementNamed("hurdle_b")) {
		VectorXd y_pos_ind = (y.array() > 0.0).cast<double>();
		if (y_pos_ind.minCoeff() < y_pos_ind.maxCoeff()) {
			ModelResult hurdle_res = fast_logistic_regression_internal(X_hurdle, y_pos_ind, Eigen::VectorXd(), R_NilValue, true, 100, 1e-8, R_NilValue, R_NilValue, optimization_alg);
			FixedParamSpec hurdle_spec = make_fixed_param_spec(X_hurdle.cols(), R_NilValue, R_NilValue);
			MatrixXd info_free = subset_matrix(hurdle_res.XtWX, hurdle_spec.free_idx, hurdle_spec.free_idx);
			MatrixXd vcov = expand_free_covariance(X_hurdle.cols(), hurdle_spec, info_free.inverse(), true);
			if (j > 0 && j <= X_hurdle.cols()) hurdle_ssq_b_j = vcov(j - 1, j - 1);
			if (X_hurdle.cols() >= 2) hurdle_ssq_b_2 = vcov(1, 1);
		}
	}

	double ssq_b_j = NA_REAL;
	double ssq_b_2 = NA_REAL;
	if (p > 0 && fit.containsElementNamed("theta_hat")) {
		SEXP theta_sexp = fit["theta_hat"];
		double theta_hat = as<double>(theta_sexp);
		if (R_finite(theta_hat) && p >= j) {
			std::vector<int> pos_rows;
			for (int i = 0; i < X.rows(); ++i) {
				if (y[i] > 0.0) pos_rows.push_back(i);
			}
			if (static_cast<int>(pos_rows.size()) > p) {
				MatrixXd X_pos(pos_rows.size(), p);
				VectorXi y_pos(pos_rows.size());
				for (size_t k = 0; k < pos_rows.size(); ++k) {
					const int i = pos_rows[k];
					X_pos.row(k) = X.row(i);
					y_pos[k] = static_cast<int>(y[i]);
				}

				VectorXd params(p + 1);
				for (int col = 0; col < p; ++col) params[col] = b_nv[col];
				params[p] = std::log(theta_hat);

				TruncatedNegBinCount fun(X_pos, y_pos);
				MatrixXd H = fun.hessian(params);
				if (H.allFinite()) {
					FixedParamSpec count_fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);
					MatrixXd H_free = subset_matrix(H, count_fixed_spec.free_idx, count_fixed_spec.free_idx);
					MatrixXd vcov = expand_free_covariance(p + 1, count_fixed_spec, H_free.inverse(), true);
					if (j > 0 && j <= p + 1) ssq_b_j = vcov(j - 1, j - 1);
					if (p >= 2) ssq_b_2 = vcov(1, 1);
				}
			}
		}
	}

	fit["ssq_b_j"] = ssq_b_j;
	fit["ssq_b_2"] = ssq_b_2;
	fit["hurdle_ssq_b_j"] = hurdle_ssq_b_j;
	fit["hurdle_ssq_b_2"] = hurdle_ssq_b_2;
	return fit;
}

// [[Rcpp::export]]
List fast_truncated_negbin_count_cpp(const Eigen::MatrixXd& X,
                                                                         const Eigen::VectorXd& y,
                                                                         Nullable<NumericVector> start_params = R_NilValue,
                                                                         bool estimate_only = false,
                                                                         int maxit = 1000,
                                                                         double tol = 1e-8,
                                                                         Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                                         Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                                         std::string optimization_alg = "lbfgs") {
        optimization_alg = normalize_optimizer_algorithm(optimization_alg, "lbfgs", false);

        List pos = build_positive_hurdle_negbin_data(X, y);
        MatrixXd X_pos = pos["X_pos"];
        VectorXi y_pos = pos["y_pos"];
        const int p = X_pos.cols();
        if (X_pos.rows() <= p){
                return List::create(
                        Named("b") = NumericVector(p, NA_REAL),
                        Named("params") = NumericVector(p + 1, NA_REAL),
                        Named("converged") = false,
                        Named("neg_ll") = NA_REAL
                );
        }

        VectorXd params(p + 1);
        if (start_params.isNotNull()) {
                params = as<VectorXd>(NumericVector(start_params));
        } else {
                params.setZero();
                params[0] = std::log(std::max(1.0, y_pos.cast<double>().mean()));
        }

        FixedParamSpec fixed_spec = make_fixed_param_spec(p + 1, fixed_idx, fixed_values);
        TruncatedNegBinCount fun(X_pos, y_pos);
        LikelihoodFitResult fit;
        try {
                fit = optimize_fixed_likelihood(fun, params, fixed_spec, maxit, tol, optimization_alg, "lbfgs");
        } catch (const std::exception& e) {
                Rcpp::warning(e.what());
                return List::create(
                        Named("b") = NumericVector(p, NA_REAL),
                        Named("params") = NumericVector(p + 1, NA_REAL),
                        Named("converged") = false,
                        Named("neg_ll") = NA_REAL
                );
        }

        params = fit.params;
        VectorXd beta = params.head(p);
        List out = List::create(
                Named("b") = beta,
                Named("params") = params,
                Named("converged") = fit.converged,
                Named("neg_ll") = fit.value
        );
        if (estimate_only || !fit.converged){
                return out;
        }

        VectorXd score = get_hurdle_negbin_count_score_cpp(X, y, params);
        MatrixXd H = fun.hessian(params);
        MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
        MatrixXd vcov = expand_free_covariance(p + 1, fixed_spec, H_free.inverse(), true);
        out["score"] = score;
        out["observed_information"] = H;
        out["information"] = H;
        out["hessian"] = -H;
        out["vcov"] = vcov;
        return out;
}
