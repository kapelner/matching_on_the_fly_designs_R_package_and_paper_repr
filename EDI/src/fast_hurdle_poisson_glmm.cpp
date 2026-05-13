// Zero-truncated Poisson GLMM for the count component of a hurdle Poisson GLMM.
//
// The hurdle likelihood factors into independent zi (logistic) and count
// (zero-truncated Poisson) components with separate random intercepts.
// This file handles the count component via Gauss-Hermite quadrature.
//
// Model (count component, positive observations only):
//   y_ij | u_i ~ TruncPoisson(exp(X_ij' beta + u_i)),  y_ij >= 1
//   u_i ~ N(0, sigma^2)
//
// Parameter vector: par = [beta_0, beta_1(treatment), ..., beta_{p-1}, log_sigma]
//   Total length: p + 1

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace Rcpp;

namespace {

struct GHRuleHP {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRuleHP gauss_hermite_rule_hp(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRuleHP rule;
	rule.nodes = es.eigenvalues();
	rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log()
		                    - 0.5 * std::log(M_PI);
	return rule;
}

inline double log_sum_exp_hp(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

inline double soft_barrier_hp(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return scale * d * d;
}

inline double soft_barrier_hp_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double soft_barrier_hp_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale;
}

// log(1 - exp(-lambda)) with numerical stability
inline double log_one_minus_exp_neg(double lambda) {
	if (lambda > 30.0) return std::log1p(-std::exp(-lambda));
	if (lambda < 1e-10) return std::log(lambda);
	return std::log1p(-std::exp(-lambda));
}

// Score d/d(eta) of log TruncPoisson(y; exp(eta)) = y - lambda / (1 - exp(-lambda))
inline double trunc_poisson_score(double y, double lambda) {
	if (lambda > 30.0) return y - lambda;
	if (lambda < 1e-8) return y - 1.0;
	return y - lambda / (1.0 - std::exp(-lambda));
}

inline double trunc_poisson_hessian(double lambda) {
	if (lambda > 30.0) return -lambda;
	if (lambda < 1e-8) return -lambda / 2.0;
	double exp_neg = std::exp(-lambda);
	double one_minus = 1.0 - exp_neg;
	return -lambda * (one_minus - lambda * exp_neg) / (one_minus * one_minus);
}

inline Eigen::ArrayXd clamp_eta_hp(const Eigen::ArrayXd& eta) {
	return eta.min(700.0);
}

inline Eigen::ArrayXd trunc_poisson_score_array(const Eigen::ArrayXd& y, const Eigen::ArrayXd& lambda) {
	const Eigen::Array<bool, Eigen::Dynamic, 1> large = (lambda > 30.0);
	const Eigen::Array<bool, Eigen::Dynamic, 1> small = (lambda < 1e-8);
	const Eigen::ArrayXd exp_neg = (-lambda).exp();
	const Eigen::ArrayXd mid = y - lambda / (1.0 - exp_neg);
	return large.select(y - lambda, small.select(y - 1.0, mid));
}

inline Eigen::ArrayXd trunc_poisson_hessian_array(const Eigen::ArrayXd& lambda) {
	const Eigen::Array<bool, Eigen::Dynamic, 1> large = (lambda > 30.0);
	const Eigen::Array<bool, Eigen::Dynamic, 1> small = (lambda < 1e-8);
	const Eigen::ArrayXd exp_neg = (-lambda).exp();
	const Eigen::ArrayXd one_minus = 1.0 - exp_neg;
	const Eigen::ArrayXd mid = -lambda * (one_minus - lambda * exp_neg) / one_minus.square();
	return large.select(-lambda, small.select(-lambda / 2.0, mid));
}

inline Eigen::ArrayXd log_one_minus_exp_neg_array(const Eigen::ArrayXd& lambda) {
	Eigen::ArrayXd out(lambda.size());
	for (Eigen::Index i = 0; i < lambda.size(); ++i) out[i] = log_one_minus_exp_neg(lambda[i]);
	return out;
}

struct HurdlePoissonGLMMData {
	Eigen::MatrixXd X_s;
	Eigen::VectorXd y_s;
	Eigen::VectorXd log_fact_y;
	std::vector<int> grp_start;
	std::vector<int> grp_size;
	int n, p, G;
	GHRuleHP gh;

	HurdlePoissonGLMMData(
		const Eigen::MatrixXd& X,
		const Eigen::VectorXd& y,
		const std::vector<int>& group_id,
		int n_gh
	) : n(X.rows()), p(X.cols()), gh(gauss_hermite_rule_hp(n_gh)) {

		std::vector<int> ord(n);
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(),
			[&](int a, int b){ return group_id[a] < group_id[b]; });

		X_s.resize(n, p);
		y_s.resize(n);
		log_fact_y.resize(n);
		for (int i = 0; i < n; ++i) {
			X_s.row(i) = X.row(ord[i]);
			y_s[i] = y[ord[i]];
			log_fact_y[i] = std::lgamma(y_s[i] + 1.0);
		}

		int prev = -1;
		for (int i = 0; i < n; ++i) {
			int g = group_id[ord[i]];
			if (g != prev) {
				grp_start.push_back(i);
				grp_size.push_back(1);
				prev = g;
			} else {
				grp_size.back()++;
			}
		}
		G = (int)grp_start.size();
	}
};

class HurdlePoissonGLMMObjective {
	const HurdlePoissonGLMMData& dat;

public:
	explicit HurdlePoissonGLMMObjective(const HurdlePoissonGLMMData& d) : dat(d) {}

	double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
		const double log_sigma = par[dat.p];
		const double sigma     = std::exp(log_sigma);
		const Eigen::VectorXd beta   = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		double total_nll = soft_barrier_hp(log_sigma);
		grad.setZero(dat.p + 1);
		const double center = 5.0, scale = 10.0;
		const double d = std::abs(log_sigma) - center;
		if (d > 0.0) grad[dat.p] += 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);

		for (int gi = 0; gi < dat.G; ++gi) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::MatrixXd Xg   = dat.X_s.middleRows(start, sz);
			const Eigen::VectorXd eta0 = Xg * beta;
			const Eigen::ArrayXd y_g = dat.y_s.segment(start, sz).array();
			const Eigen::ArrayXd log_fact_g = dat.log_fact_y.segment(start, sz).array();

			Eigen::VectorXd log_terms(n_nodes);
			std::vector<Eigen::VectorXd> lambda_nodes(n_nodes);  // lambda = exp(eta) per node

			for (int k = 0; k < n_nodes; ++k) {
				const Eigen::ArrayXd eta_k = eta0.array() + b_vals[k];
				lambda_nodes[k] = clamp_eta_hp(eta_k).exp().matrix();
				log_terms[k] = dat.gh.log_norm_weights[k] +
				               (y_g * eta_k - lambda_nodes[k].array() - log_fact_g -
				                log_one_minus_exp_neg_array(lambda_nodes[k].array())).sum();
			}

			const double ll_g = log_sum_exp_hp(log_terms);
			if (!std::isfinite(ll_g)) { grad.setZero(dat.p + 1); return 1e100; }
			total_nll -= ll_g;

			for (int k = 0; k < n_nodes; ++k) {
				const double post_k = std::exp(log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;

				Eigen::VectorXd res_k = trunc_poisson_score_array(y_g, lambda_nodes[k].array()).matrix();
				double res_sum = res_k.sum();

				grad.head(dat.p) -= post_k * (Xg.transpose() * res_k);
				grad[dat.p] -= post_k * res_sum * b_vals[k];
			}
		}
		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
		const int total = dat.p + 1;
		const double log_sigma = par[dat.p];
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total, total);
		H(dat.p, dat.p) = soft_barrier_hp_hessian(log_sigma);

		for (int gi = 0; gi < dat.G; gi++) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::MatrixXd Xg = dat.X_s.middleRows(start, sz);
			const Eigen::VectorXd eta0 = Xg * beta;
			const Eigen::ArrayXd y_g = dat.y_s.segment(start, sz).array();
			const Eigen::ArrayXd log_fact_g = dat.log_fact_y.segment(start, sz).array();

			Eigen::VectorXd log_terms(n_nodes);
			std::vector<Eigen::VectorXd> lambda_nodes(n_nodes);
			for (int k = 0; k < n_nodes; k++) {
				const Eigen::ArrayXd eta_k = eta0.array() + b_vals[k];
				lambda_nodes[k] = clamp_eta_hp(eta_k).exp().matrix();
				log_terms[k] = dat.gh.log_norm_weights[k] +
				               (y_g * eta_k - lambda_nodes[k].array() - log_fact_g -
				                log_one_minus_exp_neg_array(lambda_nodes[k].array())).sum();
			}
			const double ll_g = log_sum_exp_hp(log_terms);

			Eigen::MatrixXd E_Hik = Eigen::MatrixXd::Zero(total, total);
			Eigen::MatrixXd E_GiGiT = Eigen::MatrixXd::Zero(total, total);
			Eigen::VectorXd G_avg = Eigen::VectorXd::Zero(total);

			for (int k = 0; k < n_nodes; k++) {
				const double pk = std::exp(log_terms[k] - ll_g);
				if (pk < 1e-15) continue;

				Eigen::VectorXd G_ik = Eigen::VectorXd::Zero(total);
				Eigen::MatrixXd H_ik = Eigen::MatrixXd::Zero(total, total);
				Eigen::VectorXd res_k = trunc_poisson_score_array(y_g, lambda_nodes[k].array()).matrix();
				Eigen::VectorXd d2e_k = trunc_poisson_hessian_array(lambda_nodes[k].array()).matrix();
				const double sum_res = res_k.sum();
				const double sum_d2e = d2e_k.sum();

				G_ik.head(dat.p).noalias() = Xg.transpose() * res_k;
				H_ik.topLeftCorner(dat.p, dat.p).noalias() = weighted_crossprod(Xg, d2e_k);

				const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
				G_ik[dat.p] = sum_res * node_factor * sigma;

				H_ik(dat.p, dat.p) = (sum_d2e * node_factor * node_factor * sigma + sum_res * node_factor) * sigma;
				
				H_ik.block(0, dat.p, dat.p, 1).noalias() = (Xg.transpose() * d2e_k) * (node_factor * sigma);

				for (int r1 = 0; r1 < total; r1++) for (int c1 = 0; r1 > c1; c1++) H_ik(r1, c1) = H_ik(c1, r1);

				E_Hik += pk * H_ik;
				G_avg += pk * G_ik;
				E_GiGiT += pk * (G_ik * G_ik.transpose());
			}
			H -= (E_Hik + E_GiGiT - G_avg * G_avg.transpose());
		}
		return H;
	}
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_hurdle_poisson_glmm_score_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 7
) {
	std::vector<int> pos_idx;
	pos_idx.reserve(X.rows());
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = static_cast<int>(pos_idx.size());
	Eigen::MatrixXd X_pos(n_pos, X.cols());
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = y[i];
		gid_pos[k] = group_id[i];
	}
	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);
	HurdlePoissonGLMMObjective obj(dat);
	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	grad[X.cols()] -= soft_barrier_hp_grad(params[X.cols()]);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_hurdle_poisson_glmm_hessian_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 7
) {
	std::vector<int> pos_idx;
	pos_idx.reserve(X.rows());
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = static_cast<int>(pos_idx.size());
	Eigen::MatrixXd X_pos(n_pos, X.cols());
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = y[i];
		gid_pos[k] = group_id[i];
	}
	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);
	HurdlePoissonGLMMObjective obj(dat);
	Eigen::MatrixXd information = obj.hessian(params);
	information(X.cols(), X.cols()) -= soft_barrier_hp_hessian(params[X.cols()]);
	return -information;
}

// [[Rcpp::export]]
double get_hurdle_poisson_glmm_neg_loglik_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 7
) {
	std::vector<int> pos_idx;
	pos_idx.reserve(X.rows());
	for (int i = 0; i < X.rows(); ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = static_cast<int>(pos_idx.size());
	Eigen::MatrixXd X_pos(n_pos, X.cols());
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k] = y[i];
		gid_pos[k] = group_id[i];
	}
	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);
	HurdlePoissonGLMMObjective obj(dat);
	return likelihood_value(obj, params) - soft_barrier_hp(params[X.cols()]);
}

// [[Rcpp::export]]
List fast_hurdle_poisson_glmm_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	int j_T,
	Nullable<NumericVector> start_params = R_NilValue,
	bool estimate_only = false,
	int n_gh = 7,
	int maxit = 300,
	double eps_g = 1e-6,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
	const int n_all = X.rows();
	const int p     = X.cols();
	const int total = p + 1;

	// Filter to positive-count observations (truncated Poisson component)
	std::vector<int> pos_idx;
	pos_idx.reserve(n_all);
	for (int i = 0; i < n_all; ++i) {
		if (y[i] > 0.0) pos_idx.push_back(i);
	}
	const int n_pos = (int)pos_idx.size();

	if (n_pos <= p) {
		return List::create(
			Named("params")     = Eigen::VectorXd::Constant(total, NA_REAL),
			Named("b")          = Eigen::VectorXd::Constant(p, NA_REAL),
			Named("log_sigma")  = NA_REAL,
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	Eigen::MatrixXd X_pos(n_pos, p);
	Eigen::VectorXd y_pos(n_pos);
	std::vector<int> gid_pos(n_pos);
	for (int k = 0; k < n_pos; ++k) {
		const int i = pos_idx[k];
		X_pos.row(k) = X.row(i);
		y_pos[k]     = y[i];
		gid_pos[k]   = group_id[i];
	}

	HurdlePoissonGLMMData dat(X_pos, y_pos, gid_pos, n_gh);

	Eigen::VectorXd par(total);
	if (start_params.isNotNull()) {
		par = as<Eigen::VectorXd>(NumericVector(start_params));
	} else {
		par.setZero();
		Eigen::VectorXd log_y = y_pos.array().log().matrix();
		Eigen::VectorXd legacy_beta = Eigen::VectorXd::Zero(p);
		if (try_safe_ols_solve(X_pos, log_y, legacy_beta)) {
			par.head(p) = legacy_beta;
		}
		par[total - 1] = -3.0;
	}

	HurdlePoissonGLMMObjective obj(dat);
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	Eigen::MatrixXd info_start;
	Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	double neg_ll = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
		par       = fit.params;
		neg_ll    = fit.value;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("params")     = par,
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	const double pen        = soft_barrier_hp(par[total - 1]);
	const double true_neg_ll = neg_ll - pen;
	Eigen::VectorXd score(total);
	obj(par, score);
	score[total - 1] -= soft_barrier_hp_grad(par[total - 1]);
	score = -score;
	Eigen::MatrixXd information = obj.hessian(par);
	information(total - 1, total - 1) -= soft_barrier_hp_hessian(par[total - 1]);

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
	if (!estimate_only && converged) {
		Eigen::MatrixXd information_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::MatrixXd cov_free = covariance_from_information(information_free);
		vcov = expand_free_covariance(total, fixed_spec, cov_free, true);
		if (j_T < p) ssq_b_T = vcov(j_T, j_T);
	}

	return List::create(
		Named("params")     = par,
		Named("b")          = par.head(p),
		Named("log_sigma")  = par[total - 1],
		Named("ssq_b_T")    = ssq_b_T,
		Named("vcov")       = vcov,
		Named("score")      = score,
		Named("observed_information") = information,
		Named("information") = information,
		Named("information_type") = "observed",
		Named("hessian")    = -information,
		Named("converged")  = converged,
		Named("neg_loglik") = true_neg_ll,
		Named("neg_ll")     = true_neg_ll,
		Named("loglik")     = R_finite(true_neg_ll) ? -true_neg_ll : NA_REAL
	);
}
