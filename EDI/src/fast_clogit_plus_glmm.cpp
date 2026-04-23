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
		double total_ll = 0.0;
		int i = 0;
		while (i < group_conc.size()) {
			const int g = group_conc[i];
			int j = i + 1;
			while (j < group_conc.size() && group_conc[j] == g) ++j;
			const Eigen::VectorXd eta0 = X_conc.middleRows(i, j - i) * beta;
			Eigen::VectorXd log_terms(b_vals.size());
			for (int k = 0; k < b_vals.size(); ++k) {
				double ll = gh.log_norm_weights[k];
				for (int r = 0; r < eta0.size(); ++r) {
					const double eta = eta0[r] + b_vals[k];
					ll += y_conc[i + r] * eta - log1pexp_cpp(eta);
				}
				log_terms[k] = ll;
			}
			const double ll_g = log_sum_exp_cpp(log_terms);
			if (!std::isfinite(ll_g)) return 1e100;
			total_ll += ll_g;
			i = j;
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
			Eigen::VectorXd mu_d(eta_d.size());
			for (int i = 0; i < eta_d.size(); ++i) {
				double ed = eta_d[i];
				total_nll += log1pexp_cpp(ed) - y_disc[i] * ed;
				mu_d[i] = plogis_safe(ed);
			}

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

			int i = 0;
			while (i < group_conc.size()) {
				const int g = group_conc[i];
				int j = i + 1;
				while (j < group_conc.size() && group_conc[j] == g) ++j;
				const int sz = j - i;
				const Eigen::MatrixXd Xg = X_conc.middleRows(i, sz);
				const Eigen::VectorXd eta0 = Xg * beta;

				Eigen::VectorXd log_terms(n_nodes);
				std::vector<Eigen::VectorXd> mu_nodes(n_nodes);

				for (int k = 0; k < n_nodes; ++k) {
					double ll = gh.log_norm_weights[k];
					mu_nodes[k].resize(sz);
					for (int r = 0; r < sz; ++r) {
						const double eta = eta0[r] + b_vals[k];
						ll += y_conc[i + r] * eta - log1pexp_cpp(eta);
						mu_nodes[k][r] = plogis_safe(eta);
					}
					log_terms[k] = ll;
				}

				const double ll_g = log_sum_exp_cpp(log_terms);
				total_nll -= ll_g;

				for (int k = 0; k < n_nodes; ++k) {
					double post_k = std::exp(log_terms[k] - ll_g);
					if (post_k < 1e-15) continue;

					Eigen::VectorXd res_k(sz);
					for (int r = 0; r < sz; ++r) {
						res_k[r] = y_conc[i + r] - mu_nodes[k][r];
					}

					// dLL/dbeta
					grad.head(p_full - 1) -= post_k * (Xg.transpose() * res_k);

					// dLL/dsigma * dsigma/dlog_sigma
					double dll_dsigma = res_k.sum() * std::sqrt(2.0) * gh.nodes[k];
					grad[p_full - 1] -= post_k * dll_dsigma * sigma;
				}
				i = j;
			}
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
			Eigen::VectorXd w_d(eta_d.size());
			for (int i = 0; i < eta_d.size(); ++i) {
				const double mu = plogis_safe(eta_d[i]);
				w_d[i] = mu * (1.0 - mu);
			}
			const Eigen::MatrixXd H_clogit = X_disc.transpose() * w_d.asDiagonal() * X_disc;
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

			int i = 0;
			while (i < group_conc.size()) {
				const int g = group_conc[i];
				int j = i + 1;
				while (j < group_conc.size() && group_conc[j] == g) ++j;
				const int sz = j - i;
				const Eigen::MatrixXd Xg = X_conc.middleRows(i, sz);
				const Eigen::VectorXd yg = y_conc.segment(i, sz);
				const Eigen::VectorXd eta0 = Xg * beta;

				Eigen::VectorXd log_terms(n_nodes);
				std::vector<Eigen::VectorXd> node_score(n_nodes, Eigen::VectorXd::Zero(p_full));
				std::vector<Eigen::MatrixXd> node_hess(n_nodes, Eigen::MatrixXd::Zero(p_full, p_full));

				for (int k = 0; k < n_nodes; ++k) {
					double ll = gh.log_norm_weights[k];
					const double b = b_vals[k];
					Eigen::VectorXd mu(sz);
					Eigen::VectorXd res(sz);
					Eigen::VectorXd weight(sz);

					for (int r = 0; r < sz; ++r) {
						const double eta = eta0[r] + b;
						ll += yg[r] * eta - log1pexp_cpp(eta);
						mu[r] = plogis_safe(eta);
						res[r] = yg[r] - mu[r];
						weight[r] = mu[r] * (1.0 - mu[r]);
					}

					Eigen::VectorXd& a = node_score[k];
					Eigen::MatrixXd& B = node_hess[k];
					a.head(n_beta).noalias() = Xg.transpose() * res;
					a[p_full - 1] = b * res.sum();
					B.topLeftCorner(n_beta, n_beta).noalias() =
						-Xg.transpose() * weight.asDiagonal() * Xg;
					B.block(0, p_full - 1, n_beta, 1).noalias() =
						-b * (Xg.transpose() * weight);
					B.block(p_full - 1, 0, 1, n_beta) =
						B.block(0, p_full - 1, n_beta, 1).transpose();
					B(p_full - 1, p_full - 1) =
						b * res.sum() - b * b * weight.sum();
					log_terms[k] = ll;
				}

				const double ll_g = log_sum_exp_cpp(log_terms);
				if (!std::isfinite(ll_g)) {
					H.setConstant(NA_REAL);
					return H;
				}

				Eigen::VectorXd mean_score = Eigen::VectorXd::Zero(p_full);
				Eigen::MatrixXd mean_hess = Eigen::MatrixXd::Zero(p_full, p_full);
				Eigen::MatrixXd mean_outer = Eigen::MatrixXd::Zero(p_full, p_full);
				for (int k = 0; k < n_nodes; ++k) {
					const double post_k = std::exp(log_terms[k] - ll_g);
					if (post_k < 1e-15) continue;
					mean_score.noalias() += post_k * node_score[k];
					mean_hess.noalias() += post_k * node_hess[k];
					mean_outer.noalias() += post_k * (node_score[k] * node_score[k].transpose());
				}

				H.noalias() -= mean_hess + mean_outer - mean_score * mean_score.transpose();
				i = j;
			}
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
	std::string optimization_alg = "lbfgs"
) {
	ClogitPlusGLMMObjective obj(
		X_disc, y_disc, X_conc, y_conc, group_conc,
		has_discordant, has_concordant, 20, max_abs_log_sigma
	);

	Eigen::VectorXd par = start;
	FixedParamSpec fixed_spec = make_fixed_param_spec(par.size(), fixed_idx, fixed_values);

	double neg_ll = NA_REAL;
	int niter = maxit;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs");
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
