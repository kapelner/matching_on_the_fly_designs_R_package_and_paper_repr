#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

namespace {

inline double log_sum_exp_wf(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

inline Eigen::ArrayXd clamp_weibull_wf(const Eigen::ArrayXd& w) {
	return w.min(700.0);
}

struct GHRuleWF {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRuleWF gauss_hermite_rule_wf(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRuleWF rule;
	rule.nodes = es.eigenvalues();
	Eigen::VectorXd weights = std::sqrt(M_PI) * es.eigenvectors().row(0).array().square().matrix();
	rule.log_norm_weights = weights.array().log() - 0.5 * std::log(M_PI);
	return rule;
}

// Weibull AFT GLMM with log-Normal random intercept per cluster.
// params: [beta(p), log_sigma_eps(1), log_sigma_u(1)]
// log_sigma_eps: log of Weibull scale (error SD on log scale)
// log_sigma_u: log of random-intercept SD
class WeibullFrailtyLikelihood {
private:
	Eigen::VectorXd m_y;
	Eigen::VectorXd m_dead;
	Eigen::MatrixXd m_X;
	Eigen::VectorXd m_log_y;
	const int m_n;
	const int m_p;
	const GHRuleWF m_gh;
	const double m_max_abs_log_sigma;

	std::vector<int> m_group_start;
	std::vector<int> m_group_end;
	Eigen::VectorXd m_group_dead_sum;
	Eigen::VectorXd m_group_dead_log_y_sum;
	int m_n_groups;
	int m_max_group_size;
	Eigen::VectorXd m_eta_all;
	Eigen::VectorXd m_u_vals;
	Eigen::MatrixXd m_log_terms_mat;
	Eigen::MatrixXd m_r_all_k_mat;
	Eigen::MatrixXd m_w_all_k_mat;
	Eigen::VectorXd m_ll_g_vec;
	Eigen::VectorXd m_weighted_r;

	void build_group_structure(const Eigen::VectorXi& group_id) {
		std::vector<int> order(m_n);
		std::iota(order.begin(), order.end(), 0);
		std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
			return group_id[a] < group_id[b];
		});

		Eigen::VectorXd y_s(m_n), dead_s(m_n), log_y_s(m_n);
		Eigen::MatrixXd X_s(m_n, m_p);
		for (int i = 0; i < m_n; ++i) {
			y_s[i]     = m_y[order[i]];
			dead_s[i]  = m_dead[order[i]];
			log_y_s[i] = m_log_y[order[i]];
			X_s.row(i) = m_X.row(order[i]);
		}
		m_y     = y_s;
		m_dead  = dead_s;
		m_log_y = log_y_s;
		m_X     = X_s;

		Eigen::VectorXi sorted_gid(m_n);
		for (int i = 0; i < m_n; ++i) sorted_gid[i] = group_id[order[i]];

			auto layout = build_contiguous_group_layout(m_n, [&](int idx) { return sorted_gid[idx]; });
			m_group_start = layout.warm_start_params;
			m_group_end.resize(layout.G);
			for (int gi = 0; gi < layout.G; ++gi) {
				m_group_end[gi] = layout.warm_start_params[gi] + layout.size[gi];
			}
			m_n_groups = layout.G;
			m_max_group_size = layout.max_size;
			m_group_dead_sum.resize(m_n_groups);
			m_group_dead_log_y_sum.resize(m_n_groups);
			for (int gi = 0; gi < m_n_groups; ++gi) {
				const int warm_start_params = m_group_start[gi];
				const int sz = m_group_end[gi] - warm_start_params;
				const auto dead_seg = m_dead.segment(warm_start_params, sz);
				m_group_dead_sum[gi] = dead_seg.sum();
				m_group_dead_log_y_sum[gi] = (dead_seg.array() * m_log_y.segment(warm_start_params, sz).array()).sum();
			}
		}

public:
	WeibullFrailtyLikelihood(
		const Eigen::VectorXd& y,
		const Eigen::VectorXd& dead,
		const Eigen::MatrixXd& X,
		const Eigen::VectorXi& group_id,
		int n_gh = 20,
		double max_abs_log_sigma = 8.0
	) : m_y(y), m_dead(dead), m_X(X),
		m_log_y(y.array().log().matrix()),
		m_n(y.size()), m_p(X.cols()),
		m_gh(gauss_hermite_rule_wf(n_gh)),
		m_max_abs_log_sigma(max_abs_log_sigma),
		m_n_groups(0),
		m_max_group_size(0),
		m_eta_all(m_n),
		m_u_vals(n_gh),
		m_log_terms_mat(1, n_gh),
		m_r_all_k_mat(m_n, n_gh),
		m_w_all_k_mat(m_n, n_gh),
		m_ll_g_vec(1),
		m_weighted_r(1)
	{
		build_group_structure(group_id);
		m_u_vals.resize(m_gh.nodes.size());
		m_log_terms_mat.resize(m_n_groups, m_gh.nodes.size());
		m_ll_g_vec.resize(m_n_groups);
		m_weighted_r.resize(std::max(1, m_max_group_size));
	}

	double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
		const Eigen::VectorXd beta     = par.head(m_p);
		const double log_sigma_eps     = par[m_p];
		const double log_sigma_u       = par[m_p + 1];

		grad.setZero(m_p + 2);

		if (!std::isfinite(log_sigma_eps) || std::abs(log_sigma_eps) > m_max_abs_log_sigma) return 1e100;
		if (!std::isfinite(log_sigma_u)   || std::abs(log_sigma_u)   > m_max_abs_log_sigma) return 1e100;

		const double sigma_eps = std::exp(log_sigma_eps);
		const double sigma_u   = std::exp(log_sigma_u);
		m_eta_all.noalias() = m_X * beta;
		const int n_nodes = (int)m_gh.nodes.size();

		for (int k = 0; k < n_nodes; ++k)
			m_u_vals[k] = std::sqrt(2.0) * sigma_u * m_gh.nodes[k];

		const Eigen::ArrayXd log_y_all = m_log_y.array();
		const Eigen::ArrayXd dead_all = m_dead.array();

			for (int k = 0; k < n_nodes; ++k) {
				const Eigen::ArrayXd wik = clamp_weibull_wf((log_y_all - m_eta_all.array() - m_u_vals[k]) / sigma_eps);
				const Eigen::ArrayXd ewik = wik.exp();
				m_w_all_k_mat.col(k) = wik.matrix();
				m_r_all_k_mat.col(k) = (ewik - dead_all).matrix();
				
				for (int gi = 0; gi < m_n_groups; ++gi) {
					const int warm_start_params = m_group_start[gi];
					const int sz = m_group_end[gi] - warm_start_params;
					const double sum_dead_wik = (dead_all.segment(warm_start_params, sz) * wik.segment(warm_start_params, sz)).sum();
					const double sum_ewik = ewik.segment(warm_start_params, sz).sum();
					m_log_terms_mat(gi, k) = m_gh.log_norm_weights[k]
						+ sum_dead_wik
						- m_group_dead_sum[gi] * log_sigma_eps
						- m_group_dead_log_y_sum[gi]
						- sum_ewik;
				}
			}

		double total_nll = 0.0;
		for (int gi = 0; gi < m_n_groups; ++gi) {
			m_ll_g_vec[gi] = log_sum_exp_wf(m_log_terms_mat.row(gi));
			if (!std::isfinite(m_ll_g_vec[gi])) { grad.setZero(m_p + 2); return 1e100; }
			total_nll -= m_ll_g_vec[gi];
		}

			Eigen::VectorXd grad_beta = Eigen::VectorXd::Zero(m_p);
			double grad_log_sigma_eps = 0.0;
			double grad_log_sigma_u = 0.0;
			const double inv_sigma_eps = 1.0 / sigma_eps;

			for (int gi = 0; gi < m_n_groups; ++gi) {
				const int warm_start_params = m_group_start[gi];
				const int sz = m_group_end[gi] - warm_start_params;
				m_weighted_r.head(sz).setZero();

				for (int k = 0; k < n_nodes; ++k) {
					const double pk = std::exp(m_log_terms_mat(gi, k) - m_ll_g_vec[gi]);
					if (pk < 1e-15) continue;
					const auto r_seg = m_r_all_k_mat.col(k).segment(warm_start_params, sz);
					const auto w_seg = m_w_all_k_mat.col(k).segment(warm_start_params, sz);
					const auto dead_seg = dead_all.segment(warm_start_params, sz);

					const double r_k_sum_gi = r_seg.sum();
					const double dll_dlse_gi = (w_seg.array() * (r_seg.array() + dead_seg) - dead_seg).sum();
					
					m_weighted_r.head(sz).noalias() += pk * r_seg;
					grad_log_sigma_u -= pk * (m_u_vals[k] * inv_sigma_eps) * r_k_sum_gi;
					grad_log_sigma_eps -= pk * dll_dlse_gi;
				}

				grad_beta.noalias() -= inv_sigma_eps * m_X.middleRows(warm_start_params, sz).transpose() * m_weighted_r.head(sz);
			}
		
			grad.head(m_p) = grad_beta;
		grad[m_p] = grad_log_sigma_eps;
		grad[m_p + 1] = grad_log_sigma_u;

		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
		const Eigen::VectorXd beta     = par.head(m_p);
		const double log_sigma_eps     = par[m_p];
		const double log_sigma_u       = par[m_p + 1];
		const int n_par = m_p + 2;
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n_par, n_par);

		if (!std::isfinite(log_sigma_eps) || std::abs(log_sigma_eps) > m_max_abs_log_sigma) {
			H.setConstant(NA_REAL); return H;
		}
		if (!std::isfinite(log_sigma_u) || std::abs(log_sigma_u) > m_max_abs_log_sigma) {
			H.setConstant(NA_REAL); return H;
		}

		const double sigma_eps = std::exp(log_sigma_eps);
		const double sigma_u = std::exp(log_sigma_u);
		const Eigen::VectorXd eta_all = m_X * beta;
		const int n_nodes = (int)m_gh.nodes.size();

		Eigen::VectorXd u_vals(n_nodes);
		for (int k = 0; k < n_nodes; ++k)
			u_vals[k] = std::sqrt(2.0) * sigma_u * m_gh.nodes[k];

		const Eigen::ArrayXd log_y_all = m_log_y.array();
		const Eigen::ArrayXd dead_all = m_dead.array();

		Eigen::MatrixXd log_terms_mat(m_n_groups, n_nodes);
		std::vector<Eigen::VectorXd> w_all_k_vec(n_nodes);
		std::vector<Eigen::VectorXd> ew_all_k_vec(n_nodes);

		for (int k = 0; k < n_nodes; ++k) {
			const Eigen::ArrayXd wik = clamp_weibull_wf((log_y_all - eta_all.array() - u_vals[k]) / sigma_eps);
			const Eigen::ArrayXd ewik = wik.exp();
			w_all_k_vec[k] = wik.matrix();
			ew_all_k_vec[k] = ewik.matrix();
			const Eigen::ArrayXd term_all_k = (dead_all > 0.5).select(wik - log_sigma_eps - log_y_all - ewik, -ewik);
			for (int gi = 0; gi < m_n_groups; ++gi) {
				log_terms_mat(gi, k) = m_gh.log_norm_weights[k] + 
				                       term_all_k.segment(m_group_start[gi], m_group_end[gi] - m_group_start[gi]).sum();
			}
		}

		Eigen::VectorXd ll_g_vec(m_n_groups);
		for (int gi = 0; gi < m_n_groups; ++gi) ll_g_vec[gi] = log_sum_exp_wf(log_terms_mat.row(gi));

		Eigen::MatrixXd E_Hik_sum = Eigen::MatrixXd::Zero(n_par, n_par);
		Eigen::MatrixXd E_GiGiT_sum = Eigen::MatrixXd::Zero(n_par, n_par);
		Eigen::MatrixXd G_avg_outer_sum = Eigen::MatrixXd::Zero(n_par, n_par);

		for (int k = 0; k < n_nodes; k++) {
			Eigen::VectorXd pk_vec(m_n_groups);
			for (int gi = 0; gi < m_n_groups; gi++) pk_vec[gi] = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
			
			Eigen::VectorXd pk_expanded(m_n);
			for (int gi = 0; gi < m_n_groups; gi++) 
				pk_expanded.segment(m_group_start[gi], m_group_end[gi] - m_group_start[gi]).setConstant(pk_vec[gi]);
			const double u = u_vals[k];
			const double v = u / sigma_eps;

			// Beta-Beta block
			Eigen::VectorXd w_Hik_beta = pk_expanded.cwiseProduct(ew_all_k_vec[k]) / (sigma_eps * sigma_eps);
			E_Hik_sum.topLeftCorner(m_p, m_p).noalias() -= weighted_crossprod(m_X, w_Hik_beta);

			// Sigma-Sigma and cross blocks
			for (int gi = 0; gi < m_n_groups; gi++) {
				const double pk = pk_vec[gi];
				if (pk < 1e-15) continue;
				const int warm_start_params = m_group_start[gi];
				const int sz = m_group_end[gi] - warm_start_params;
				
				const double sum_ew_k_gi = ew_all_k_vec[k].segment(warm_start_params, sz).sum();
				const double sum_w_ew_k_gi = (w_all_k_vec[k].array().segment(warm_start_params, sz) * ew_all_k_vec[k].array().segment(warm_start_params, sz)).sum();
				const double sum_w2_ew_k_gi = (w_all_k_vec[k].array().segment(warm_start_params, sz).square() * ew_all_k_vec[k].array().segment(warm_start_params, sz)).sum();
				const double sum_res_k_gi = (ew_all_k_vec[k].array().segment(warm_start_params, sz) - dead_all.segment(warm_start_params, sz)).sum();
				
				// Sigma_eps - Sigma_eps
				E_Hik_sum(m_p, m_p) += pk * (-sum_w2_ew_k_gi - 2.0 * sum_w_ew_k_gi + sum_res_k_gi);
				// Sigma_u - Sigma_u
				E_Hik_sum(m_p + 1, m_p + 1) += pk * (-(v * v) * sum_ew_k_gi + (v / sigma_u) * sum_res_k_gi * std::sqrt(2.0) * m_gh.nodes[k] * sigma_u);
				// Sigma_eps - Sigma_u
				E_Hik_sum(m_p, m_p + 1) += pk * (-v * (sum_w_ew_k_gi + sum_ew_k_gi));
				
				// Beta-Sigma blocks
				Eigen::VectorXd Xg_ew = m_X.middleRows(warm_start_params, sz).transpose() * ew_all_k_vec[k].segment(warm_start_params, sz);
				Eigen::VectorXd Xg_w_ew = m_X.middleRows(warm_start_params, sz).transpose() * (w_all_k_vec[k].array().segment(warm_start_params, sz) * ew_all_k_vec[k].array().segment(warm_start_params, sz)).matrix();
				
				E_Hik_sum.block(0, m_p, m_p, 1).noalias() += pk * (Xg_w_ew + Xg_ew) / sigma_eps;
				E_Hik_sum.block(0, m_p + 1, m_p, 1).noalias() += pk * (v / sigma_eps) * Xg_ew;
			}
		}

		for (int gi = 0; gi < m_n_groups; gi++) {
			Eigen::VectorXd G_avg_gi = Eigen::VectorXd::Zero(n_par);
			Eigen::MatrixXd E_GiGiT_gi = Eigen::MatrixXd::Zero(n_par, n_par);
			const int warm_start_params = m_group_start[gi];
			const int sz = m_group_end[gi] - warm_start_params;

			for (int k = 0; k < n_nodes; k++) {
				const double pk = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
				if (pk < 1e-15) continue;
				
				Eigen::VectorXd G_ik = Eigen::VectorXd::Zero(n_par);
				Eigen::VectorXd res_k_gi = (ew_all_k_vec[k].array().segment(warm_start_params, sz) - dead_all.segment(warm_start_params, sz)).matrix();
				G_ik.head(m_p).noalias() = m_X.middleRows(warm_start_params, sz).transpose() * res_k_gi / sigma_eps;
				G_ik[m_p] = (w_all_k_vec[k].array().segment(warm_start_params, sz) * ew_all_k_vec[k].array().segment(warm_start_params, sz) - dead_all.segment(warm_start_params, sz)).sum();
				G_ik[m_p + 1] = (u_vals[k] / sigma_eps) * res_k_gi.sum();

				G_avg_gi.noalias() += pk * G_ik;
				E_GiGiT_gi.noalias() += pk * (G_ik * G_ik.transpose());
			}
			E_GiGiT_sum.noalias() += E_GiGiT_gi;
			G_avg_outer_sum.noalias() += G_avg_gi * G_avg_gi.transpose();
		}

		E_Hik_sum.block(m_p, 0, 1, m_p) = E_Hik_sum.block(0, m_p, m_p, 1).transpose();
		E_Hik_sum.block(m_p + 1, 0, 1, m_p) = E_Hik_sum.block(0, m_p + 1, m_p, 1).transpose();
		E_Hik_sum(m_p + 1, m_p) = E_Hik_sum(m_p, m_p + 1);

		H -= (E_Hik_sum + E_GiGiT_sum - G_avg_outer_sum);
		return (H + H.transpose()) / 2.0;
	}
};

} // namespace

// [[Rcpp::export]]
double get_weibull_frailty_neg_loglik_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& dead,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);
	Eigen::VectorXd grad(params.size());
	return obj(params, grad);
}

// [[Rcpp::export]]
Eigen::VectorXd get_weibull_frailty_score_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& dead,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);
	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_weibull_frailty_hessian_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& dead,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);
	return -obj.hessian(params);
}

// [[Rcpp::export]]
List fast_weibull_frailty_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& dead,
	const Eigen::VectorXi& group_id,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_params = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
	bool estimate_only = false,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
	const int p     = X.cols();
	const int n_par = p + 2;

	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);

	Eigen::VectorXd par(n_par);
	if (warm_start_params.isNotNull()) {
		par = Rcpp::as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_params));
	} else if (warm_start_beta.isNotNull()) {
		VectorXd sb = as<VectorXd>(warm_start_beta);
		if (sb.size() == n_par) {
			par = sb;
		} else if (sb.size() == p) {
			par.head(p) = sb;
			par[p] = 0.0;
			par[p+1] = -3.0;
		}
	} else {
		Eigen::VectorXd log_y = y.array().log().matrix();
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
		par.head(p) = cod.solve(log_y);
		Eigen::VectorXd resid = log_y - X * par.head(p);
		const double std_resid = std::sqrt(resid.squaredNorm() / std::max(1, (int)y.size() - p));
		par[p]     = std::log(std_resid * 0.7797);
		par[p + 1] = 0.0;
	}

	FixedParamSpec fixed_spec = make_fixed_param_spec(n_par, fixed_idx, fixed_values);

	Eigen::MatrixXd info_start;
	Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	double neg_ll  = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
		par       = fit.params;
		neg_ll    = fit.value;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("b")           = par.head(p),
			Named("log_sigma_eps") = NA_REAL,
			Named("log_sigma_u") = NA_REAL,
			Named("ssq_b_T")     = NA_REAL,
			Named("converged")   = false,
			Named("neg_loglik")  = NA_REAL
		);
	}

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(n_par, n_par, NA_REAL);
	Eigen::VectorXd score = Eigen::VectorXd::Constant(n_par, NA_REAL);
	Eigen::MatrixXd information = Eigen::MatrixXd::Constant(n_par, n_par, NA_REAL);
	if (!estimate_only && converged) {
		Eigen::VectorXd grad(n_par);
		obj(par, grad);
		score = -grad;
		information = obj.hessian(par);
		Eigen::MatrixXd info_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(info_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(info_free.rows(), info_free.cols()));
			vcov = expand_free_covariance(n_par, fixed_spec, inv_free, true);
			if (vcov.allFinite()) ssq_b_T = vcov(0, 0);  // j_T = 0 usually
		}
	}

	return List::create(
		Named("params")        = par,
		Named("b")             = par.head(p),
		Named("log_sigma_eps") = par[p],
		Named("log_sigma_u")   = par[p + 1],
		Named("ssq_b_T")       = ssq_b_T,
		Named("vcov")          = vcov,
		Named("score")         = score,
		Named("observed_information") = information,
		Named("information")   = information,
		Named("information_type") = "observed",
		Named("hessian")       = -information,
		Named("converged")     = converged,
		Named("neg_loglik")    = neg_ll,
		Named("neg_ll")        = neg_ll,
		Named("loglik")        = R_finite(neg_ll) ? -neg_ll : NA_REAL
	);
}
