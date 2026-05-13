#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

inline double log1pexp_cpp(double x) {
	if (x > 0.0) return x + std::log1p(std::exp(-x));
	return std::log1p(std::exp(x));
}

inline double log_sum_exp_cpp(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

struct GHRule {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRule gauss_hermite_rule(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	Eigen::VectorXd nodes = es.eigenvalues();
	Eigen::VectorXd weights = std::sqrt(M_PI) * es.eigenvectors().row(0).array().square().matrix();
	GHRule rule;
	rule.nodes = nodes;
	rule.log_norm_weights = weights.array().log() - 0.5 * std::log(M_PI);
	return rule;
}

class ClogitPlusGLMMObjective {
private:
	const Eigen::MatrixXd X_disc;
	const Eigen::VectorXd y_disc;
	const Eigen::MatrixXd X_conc;
	const Eigen::VectorXd y_conc;
	const Eigen::VectorXi group_conc;
	const int q;
	const bool has_discordant;
	const bool has_concordant;
	const GHRule gh;
	const double max_abs_log_sigma;

public:
	ClogitPlusGLMMObjective(
		const Eigen::MatrixXd& X_disc_,
		const Eigen::VectorXd& y_disc_,
		const Eigen::MatrixXd& X_conc_,
		const Eigen::VectorXd& y_conc_,
		const Eigen::VectorXi& group_conc_,
		const bool has_discordant_,
		const bool has_concordant_,
		const int n_gh = 20,
		const double max_abs_log_sigma_ = 8.0
	) :
		X_disc(X_disc_), y_disc(y_disc_), X_conc(X_conc_), y_conc(y_conc_),
		group_conc(group_conc_), q(has_discordant_ ? (int)X_disc_.cols() : (int)X_conc_.cols()), has_discordant(has_discordant_),
		has_concordant(has_concordant_), gh(gauss_hermite_rule(n_gh)),
		max_abs_log_sigma(max_abs_log_sigma_) {}

	double neg_clogit(const Eigen::VectorXd& beta_no_intercept) const {
		if (!has_discordant) return 0.0;
		const Eigen::VectorXd eta = X_disc * beta_no_intercept;
		double nll = 0.0;
		for (int i = 0; i < eta.size(); ++i) {
			const double ll = y_disc[i] * eta[i] - log1pexp_cpp(eta[i]);
			if (!std::isfinite(ll)) return 1e100;
			nll -= ll;
		}
		return nll;
	}

	double neg_glmm(const Eigen::VectorXd& par_full) const {
		if (!has_concordant) return 0.0;
		const int p_full = par_full.size();
		const double log_sigma = par_full[p_full - 1];
		if (!std::isfinite(log_sigma) || std::abs(log_sigma) > max_abs_log_sigma) return 1e100;
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par_full.head(p_full - 1);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * gh.nodes;
		const int n_nodes = (int)b_vals.size();

		const Eigen::VectorXd eta_conc_all = X_conc * beta;
		const Eigen::ArrayXd y_conc_all = y_conc.array();

		std::vector<int> grp_start, grp_size;
		int i_idx = 0;
		while (i_idx < group_conc.size()) {
			const int g = group_conc[i_idx];
			int j_idx = i_idx + 1;
			while (j_idx < group_conc.size() && group_conc[j_idx] == g) ++j_idx;
			grp_start.push_back(i_idx);
			grp_size.push_back(j_idx - i_idx);
			i_idx = j_idx;
		}
		const int G_conc = (int)grp_start.size();

		Eigen::MatrixXd log_terms_mat(G_conc, n_nodes);
		for (int k = 0; k < n_nodes; ++k) {
			const Eigen::ArrayXd eta_k_all = eta_conc_all.array() + b_vals[k];
			const Eigen::ArrayXd term_k_all = y_conc_all * eta_k_all - log1pexp_array_safe(eta_k_all);
			for (int gi = 0; gi < G_conc; ++gi) {
				log_terms_mat(gi, k) = gh.log_norm_weights[k] + 
				                       term_k_all.segment(grp_start[gi], grp_size[gi]).sum();
			}
		}

		double total_ll = 0.0;
		for (int gi = 0; gi < G_conc; ++gi) {
			const double ll_g = log_sum_exp_cpp(log_terms_mat.row(gi));
			if (!std::isfinite(ll_g)) return 1e100;
			total_ll += ll_g;
		}
		return -total_ll;
	}

	double value(const Eigen::VectorXd& par) const {
		double out = 0.0;
		if (has_discordant) {
			const Eigen::VectorXd beta_no_intercept =
				has_concordant ? par.segment(1, q) : par.head(q);
			out += neg_clogit(beta_no_intercept);
		}
		if (has_concordant) out += neg_glmm(par);
		if (!std::isfinite(out)) return 1e100;
		return out;
	}

	double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
		const int p_full = par.size();
		grad.setZero(p_full);
		double total_nll = 0.0;

		// 1. Conditional Logistic component
		if (has_discordant) {
			const Eigen::VectorXd beta_no_intercept =
				has_concordant ? par.segment(1, q) : par.head(q);
			
			const Eigen::VectorXd eta_d = X_disc * beta_no_intercept;
			const Eigen::VectorXd mu_d = plogis_array_safe(eta_d.array()).matrix();
			total_nll += (log1pexp_array_safe(eta_d.array()) - y_disc.array() * eta_d.array()).sum();

			Eigen::VectorXd grad_clogit = X_disc.transpose() * (mu_d - y_disc);
			if (has_concordant) {
				grad.segment(1, q) += grad_clogit;
			} else {
				grad.head(q) += grad_clogit;
			}
		}

		// 2. GLMM component
		if (has_concordant) {
			const double log_sigma = par[p_full - 1];
			const double sigma = std::exp(log_sigma);
			const Eigen::VectorXd beta = par.head(p_full - 1);
			const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * gh.nodes;
			const int n_nodes = (int)b_vals.size();

			const int n_conc = X_conc.rows();
			const Eigen::VectorXd eta_conc_all = X_conc * beta;
			const Eigen::ArrayXd y_conc_all = y_conc.array();

			std::vector<int> grp_start, grp_size;
			int i_idx = 0;
			while (i_idx < group_conc.size()) {
				const int g = group_conc[i_idx];
				int j_idx = i_idx + 1;
				while (j_idx < group_conc.size() && group_conc[j_idx] == g) ++j_idx;
				grp_start.push_back(i_idx);
				grp_size.push_back(j_idx - i_idx);
				i_idx = j_idx;
			}
			const int G_conc = (int)grp_start.size();

			Eigen::MatrixXd log_terms_mat(G_conc, n_nodes);
			std::vector<Eigen::VectorXd> mu_conc_all_k_vec(n_nodes);

			for (int k = 0; k < n_nodes; ++k) {
				const Eigen::ArrayXd eta_k_all = eta_conc_all.array() + b_vals[k];
				mu_conc_all_k_vec[k] = plogis_array_safe(eta_k_all).matrix();
				const Eigen::ArrayXd term_k_all = y_conc_all * eta_k_all - log1pexp_array_safe(eta_k_all);
				for (int gi = 0; gi < G_conc; ++gi) {
					log_terms_mat(gi, k) = gh.log_norm_weights[k] + 
					                       term_k_all.segment(grp_start[gi], grp_size[gi]).sum();
				}
			}

			Eigen::VectorXd ll_g_vec(G_conc);
			for (int gi = 0; gi < G_conc; ++gi) {
				ll_g_vec[gi] = log_sum_exp_cpp(log_terms_mat.row(gi));
				total_nll -= ll_g_vec[gi];
			}

			Eigen::VectorXd grad_beta_conc = Eigen::VectorXd::Zero(p_full - 1);
			double grad_log_sigma_conc = 0.0;

			for (int k = 0; k < n_nodes; ++k) {
				Eigen::VectorXd post_k_expanded(n_conc);
				double dLL_dlog_sigma_k = 0.0;
				for (int gi = 0; gi < G_conc; ++gi) {
					const double pk = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
					if (pk < 1e-15) {
						post_k_expanded.segment(grp_start[gi], grp_size[gi]).setZero();
						continue;
					}
					post_k_expanded.segment(grp_start[gi], grp_size[gi]).setConstant(pk);
					const double res_sum_k_gi = (y_conc_all.array().segment(grp_start[gi], grp_size[gi]) - 
					                             mu_conc_all_k_vec[k].array().segment(grp_start[gi], grp_size[gi])).sum();

					dLL_dlog_sigma_k += pk * res_sum_k_gi * b_vals[k];
				}
				Eigen::VectorXd res_all_k = y_conc_all.matrix() - mu_conc_all_k_vec[k];
				grad_beta_conc.noalias() -= X_conc.transpose() * (post_k_expanded.cwiseProduct(res_all_k));
				grad_log_sigma_conc -= dLL_dlog_sigma_k;
			}
			grad.head(p_full - 1) += grad_beta_conc;
			grad[p_full - 1] += grad_log_sigma_conc;
		}

		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
		const int p_full = par.size();
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(p_full, p_full);

		// 1. Conditional Logistic component
		if (has_discordant) {
			const Eigen::VectorXd beta_no_intercept =
				has_concordant ? par.segment(1, q) : par.head(q);
			const Eigen::VectorXd eta_d = X_disc * beta_no_intercept;
			const Eigen::VectorXd mu_d = plogis_array_safe(eta_d.array()).matrix();
			Eigen::VectorXd w_d = (mu_d.array() * (1.0 - mu_d.array())).matrix();
			const Eigen::MatrixXd H_clogit = weighted_crossprod(X_disc, w_d);
			if (has_concordant) {
				H.block(1, 1, q, q).noalias() += H_clogit;
			} else {
				H.topLeftCorner(q, q).noalias() += H_clogit;
			}
		}

		// 2. GLMM component
		if (has_concordant) {
			const double log_sigma = par[p_full - 1];
			if (!std::isfinite(log_sigma) || std::abs(log_sigma) > max_abs_log_sigma) {
				H.setConstant(NA_REAL);
				return H;
			}

			const double sigma = std::exp(log_sigma);
			const Eigen::VectorXd beta = par.head(p_full - 1);
			const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * gh.nodes;
			const int n_nodes = (int)b_vals.size();
			const int n_beta = p_full - 1;
			const int n_conc = X_conc.rows();
			const Eigen::VectorXd eta_conc_all = X_conc * beta;
			const Eigen::ArrayXd y_conc_all = y_conc.array();

			std::vector<int> grp_start, grp_size;
			int i_idx = 0;
			while (i_idx < group_conc.size()) {
				const int g = group_conc[i_idx];
				int j_idx = i_idx + 1;
				while (j_idx < group_conc.size() && group_conc[j_idx] == g) ++j_idx;
				grp_start.push_back(i_idx);
				grp_size.push_back(j_idx - i_idx);
				i_idx = j_idx;
			}
			const int G_conc = (int)grp_start.size();

			Eigen::MatrixXd log_terms_mat(G_conc, n_nodes);
			std::vector<Eigen::VectorXd> mu_conc_all_k_vec(n_nodes);
			std::vector<Eigen::VectorXd> w_conc_all_k_vec(n_nodes);

			for (int k = 0; k < n_nodes; ++k) {
				const Eigen::ArrayXd eta_k_all = eta_conc_all.array() + b_vals[k];
				mu_conc_all_k_vec[k] = plogis_array_safe(eta_k_all).matrix();
				w_conc_all_k_vec[k] = (mu_conc_all_k_vec[k].array() * (1.0 - mu_conc_all_k_vec[k].array())).matrix();
				const Eigen::ArrayXd term_k_all = y_conc_all * eta_k_all - log1pexp_array_safe(eta_k_all);
				for (int gi = 0; gi < G_conc; ++gi) {
					log_terms_mat(gi, k) = gh.log_norm_weights[k] + 
					                       term_k_all.segment(grp_start[gi], grp_size[gi]).sum();
				}
			}

			Eigen::VectorXd ll_g_vec(G_conc);
			for (int gi = 0; gi < G_conc; ++gi) ll_g_vec[gi] = log_sum_exp_cpp(log_terms_mat.row(gi));

			Eigen::MatrixXd E_Hik_sum = Eigen::MatrixXd::Zero(p_full, p_full);
			Eigen::MatrixXd E_GiGiT_sum = Eigen::MatrixXd::Zero(p_full, p_full);
			Eigen::MatrixXd G_avg_outer_sum = Eigen::MatrixXd::Zero(p_full, p_full);

			for (int k = 0; k < n_nodes; k++) {
				Eigen::VectorXd pk_vec(G_conc);
				for (int gi = 0; gi < G_conc; gi++) pk_vec[gi] = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);

				Eigen::VectorXd pk_expanded(n_conc);
				for (int gi = 0; gi < G_conc; gi++) pk_expanded.segment(grp_start[gi], grp_size[gi]).setConstant(pk_vec[gi]);

				Eigen::VectorXd w_Hik_beta = pk_expanded.cwiseProduct(w_conc_all_k_vec[k]);
				E_Hik_sum.topLeftCorner(n_beta, n_beta).noalias() -= weighted_crossprod(X_conc, w_Hik_beta);

				for (int gi = 0; gi < G_conc; gi++) {
					const double pk = pk_vec[gi];
					if (pk < 1e-15) continue;
					const int start = grp_start[gi];
					const int sz = grp_size[gi];
					const double sum_w_k_gi = w_conc_all_k_vec[k].segment(start, sz).sum();
					const double sum_res_k_gi = (y_conc_all.array().segment(start, sz) - mu_conc_all_k_vec[k].array().segment(start, sz)).sum();
					const double b = b_vals[k];
					
					E_Hik_sum(p_full - 1, p_full - 1) += pk * (b * sum_res_k_gi - b * b * sum_w_k_gi);
					
					Eigen::VectorXd d2L_db_dsigma_k_gi = -b * (X_conc.middleRows(start, sz).transpose() * w_conc_all_k_vec[k].segment(start, sz));
					E_Hik_sum.block(0, p_full - 1, n_beta, 1).noalias() += pk * d2L_db_dsigma_k_gi;
				}
			}

			for (int gi = 0; gi < G_conc; gi++) {
				Eigen::VectorXd G_avg_gi = Eigen::VectorXd::Zero(p_full);
				Eigen::MatrixXd E_GiGiT_gi = Eigen::MatrixXd::Zero(p_full, p_full);
				const int start = grp_start[gi];
				const int sz = grp_size[gi];
				const Eigen::MatrixXd Xg = X_conc.middleRows(start, sz);
				const Eigen::VectorXd yg = y_conc_all.segment(start, sz);

				for (int k = 0; k < n_nodes; k++) {
					const double pk = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
					if (pk < 1e-15) continue;
					
					Eigen::VectorXd G_ik = Eigen::VectorXd::Zero(p_full);
					Eigen::VectorXd res_k_gi = (yg.array() - mu_conc_all_k_vec[k].array().segment(start, sz)).matrix();
					G_ik.head(n_beta).noalias() = Xg.transpose() * res_k_gi;
					const double b = b_vals[k];
					G_ik[p_full - 1] = res_k_gi.sum() * b;

					G_avg_gi.noalias() += pk * G_ik;
					E_GiGiT_gi.noalias() += pk * (G_ik * G_ik.transpose());
				}
				E_GiGiT_sum.noalias() += E_GiGiT_gi;
				G_avg_outer_sum.noalias() += G_avg_gi * G_avg_gi.transpose();
			}

			for (int r = 0; r < n_beta; r++) for (int c = 0; c < r; c++) E_Hik_sum(r, c) = E_Hik_sum(c, r);
			E_Hik_sum.block(p_full - 1, 0, 1, n_beta) = E_Hik_sum.block(0, p_full - 1, n_beta, 1).transpose();

			H.noalias() -= (E_Hik_sum + E_GiGiT_sum - G_avg_outer_sum);
		}

		return 0.5 * (H + H.transpose());
	}
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_clogit_plus_glmm_score_cpp(
	const Eigen::MatrixXd& X_disc,
	const Eigen::VectorXd& y_disc,
	const Eigen::MatrixXd& X_conc,
	const Eigen::VectorXd& y_conc,
	const Eigen::VectorXi& group_conc,
	const Eigen::VectorXd& params,
	bool has_discordant,
	bool has_concordant,
	double max_abs_log_sigma = 8.0
) {
	ClogitPlusGLMMObjective obj(
		X_disc, y_disc, X_conc, y_conc, group_conc,
		has_discordant, has_concordant, 20, max_abs_log_sigma
	);
	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_clogit_plus_glmm_hessian_cpp(
	const Eigen::MatrixXd& X_disc,
	const Eigen::VectorXd& y_disc,
	const Eigen::MatrixXd& X_conc,
	const Eigen::VectorXd& y_conc,
	const Eigen::VectorXi& group_conc,
	const Eigen::VectorXd& params,
	bool has_discordant,
	bool has_concordant,
	double max_abs_log_sigma = 8.0
) {
	ClogitPlusGLMMObjective obj(
		X_disc, y_disc, X_conc, y_conc, group_conc,
		has_discordant, has_concordant, 20, max_abs_log_sigma
	);
	return -obj.hessian(params);
}

// [[Rcpp::export]]
List fast_clogit_plus_glmm_cpp(
	const Eigen::MatrixXd& X_disc,
	const Eigen::VectorXd& y_disc,
	const Eigen::MatrixXd& X_conc,
	const Eigen::VectorXd& y_conc,
	const Eigen::VectorXi& group_conc,
	const Eigen::VectorXd& start,
	bool has_discordant,
	bool has_concordant,
	bool estimate_only = false,
	double max_abs_log_sigma = 8.0,
	int maxit = 200,
	double eps_g = 1e-5,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
	ClogitPlusGLMMObjective obj(
		X_disc, y_disc, X_conc, y_conc, group_conc,
		has_discordant, has_concordant, 20, max_abs_log_sigma
	);

	Eigen::VectorXd par = start;
	FixedParamSpec fixed_spec = make_fixed_param_spec(par.size(), fixed_idx, fixed_values);

	Eigen::MatrixXd info_start;
	Eigen::MatrixXd* info_start_ptr = nullptr;
	if (warm_start_fisher_info.isNotNull()) {
		info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
		info_start_ptr = &info_start;
	}

	double neg_ll = NA_REAL;
	int niter = maxit;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs", 0, info_start_ptr);
		par = fit.params;
		neg_ll = fit.value;
		niter = fit.niter;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("params") = par,
			Named("b") = par,
			Named("beta_T") = NA_REAL,
			Named("se_beta_T") = NA_REAL,
			Named("ssq_b_j") = NA_REAL,
			Named("converged") = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	const int j_beta_T = has_concordant ? 1 : 0; // 0-based
	double ssq_b_j = NA_REAL;
	Eigen::MatrixXd info = obj.hessian(par);
	Eigen::VectorXd score(par.size());
	obj(par, score);
	score = -score;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(par.size(), par.size(), NA_REAL);
	if (!estimate_only && converged) {
		Eigen::MatrixXd info_free = subset_matrix(info, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(info_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(info_free.rows(), info_free.cols()));
			vcov = expand_free_covariance(par.size(), fixed_spec, inv_free, true);
			if (vcov.allFinite()) ssq_b_j = vcov(j_beta_T, j_beta_T);
		}
	}
	double se_beta_T = (std::isfinite(ssq_b_j) && ssq_b_j > 0.0) ? std::sqrt(ssq_b_j) : NA_REAL;

	return List::create(
		Named("params") = par,
		Named("b") = par,
		Named("beta_T") = par[j_beta_T],
		Named("se_beta_T") = se_beta_T,
		Named("ssq_b_j") = ssq_b_j,
		Named("vcov") = vcov,
		Named("score") = score,
		Named("observed_information") = info,
		Named("information") = info,
		Named("information_type") = "observed",
		Named("hessian") = -info,
		Named("converged") = converged,
		Named("neg_loglik") = neg_ll,
		Named("neg_ll") = neg_ll,
		Named("loglik") = R_finite(neg_ll) ? -neg_ll : NA_REAL
	);
}
