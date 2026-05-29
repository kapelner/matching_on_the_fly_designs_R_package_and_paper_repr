// Poisson GLMM for KK designs via Gauss-Hermite quadrature.
//
// Model:  log E[Y_ij | u_i] = X_ij' beta + u_i
//   u_i ~ N(0, sigma^2)    (random intercept per matched pair / singleton)
//   y_ij in {0, 1, 2, ...} count responses
//   X includes an intercept column
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

struct GHRuleP {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRuleP gauss_hermite_rule_poisson(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRuleP rule;
	rule.nodes = es.eigenvalues();
	rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log()
		                    - 0.5 * std::log(M_PI);
	return rule;
}

inline double log_sum_exp_p(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

inline double soft_barrier(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return scale * d * d;
}

inline double soft_barrier_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double soft_barrier_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale;
}

inline Eigen::ArrayXd clamp_eta_pg(const Eigen::ArrayXd& eta) {
	return eta.min(700.0);
}

struct PoissonGLMMData {
	Eigen::MatrixXd X_s;
	Eigen::VectorXd y_s;
	Eigen::VectorXd log_fact_y;  // precomputed log(y!) = lgamma(y+1)
	std::vector<int> grp_start;
	std::vector<int> grp_size;
	int n, p, G;
	GHRuleP gh;

	PoissonGLMMData(
		const Eigen::Ref<const Eigen::MatrixXd>& X,
		const Eigen::Ref<const Eigen::VectorXd>& y,
		const std::vector<int>& group_id,
		int n_gh
	) : n(X.rows()), p(X.cols()), gh(gauss_hermite_rule_poisson(n_gh)) {

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

class PoissonGLMMObjective {
	const PoissonGLMMData& dat;

public:
	explicit PoissonGLMMObjective(const PoissonGLMMData& d) : dat(d) {}

	double operator()(const Eigen::Ref<const Eigen::VectorXd>& par, Eigen::VectorXd& grad) {
		const double log_sigma = par[dat.p];
		const double sigma     = std::exp(log_sigma);
		const Eigen::VectorXd beta   = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		const Eigen::VectorXd eta_all = dat.X_s * beta;
		const Eigen::ArrayXd y_all = dat.y_s.array();
		
		Eigen::MatrixXd log_terms_mat(dat.G, n_nodes);
		std::vector<Eigen::VectorXd> mu_all_k_vec(n_nodes);

		for (int k = 0; k < n_nodes; ++k) {
			const Eigen::ArrayXd eta_all_k = eta_all.array() + b_vals[k];
			mu_all_k_vec[k] = clamp_eta_pg(eta_all_k).exp().matrix();
			const Eigen::ArrayXd term_all_k = y_all * eta_all_k - mu_all_k_vec[k].array() - dat.log_fact_y.array();
			
			for (int gi = 0; gi < dat.G; ++gi) {
				log_terms_mat(gi, k) = dat.gh.log_norm_weights[k] + 
				                       term_all_k.segment(dat.grp_start[gi], dat.grp_size[gi]).sum();
			}
		}

		Eigen::VectorXd ll_g_vec(dat.G);
		double total_nll = soft_barrier(log_sigma);
		for (int gi = 0; gi < dat.G; ++gi) {
			ll_g_vec[gi] = log_sum_exp_p(log_terms_mat.row(gi));
			if (!std::isfinite(ll_g_vec[gi])) { grad.setZero(dat.p + 1); return 1e100; }
			total_nll -= ll_g_vec[gi];
		}

		grad.setZero(dat.p + 1);
		const double center = 5.0, scale = 10.0;
		const double d_pen = std::abs(log_sigma) - center;
		if (d_pen > 0.0) grad[dat.p] += 2.0 * scale * d_pen * (log_sigma > 0 ? 1.0 : -1.0);

		Eigen::VectorXd grad_beta = Eigen::VectorXd::Zero(dat.p);
		double grad_log_sigma = 0.0;

		for (int k = 0; k < n_nodes; ++k) {
			Eigen::VectorXd post_k_expanded(dat.n);
			double dLL_dlog_sigma_k = 0.0;
			
			for (int gi = 0; gi < dat.G; ++gi) {
				const double pk = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
				if (pk < 1e-15) {
					post_k_expanded.segment(dat.grp_start[gi], dat.grp_size[gi]).setZero();
					continue;
				}
				post_k_expanded.segment(dat.grp_start[gi], dat.grp_size[gi]).setConstant(pk);
				
				const double res_sum_k_gi = (y_all.array().segment(dat.grp_start[gi], dat.grp_size[gi]) - 
				                             mu_all_k_vec[k].array().segment(dat.grp_start[gi], dat.grp_size[gi])).sum();
				dLL_dlog_sigma_k += pk * res_sum_k_gi * b_vals[k];
			}
			
			Eigen::VectorXd res_all_k = y_all.matrix() - mu_all_k_vec[k];
			grad_beta.noalias() -= dat.X_s.transpose() * (post_k_expanded.cwiseProduct(res_all_k));
			grad_log_sigma -= dLL_dlog_sigma_k;
		}
		
		grad.head(dat.p) = grad_beta;
		grad[dat.p] += grad_log_sigma;

		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::Ref<const Eigen::VectorXd>& par) {
		const int total = dat.p + 1;
		const double log_sigma = par[dat.p];
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		const Eigen::VectorXd eta_all = dat.X_s * beta;
		const Eigen::ArrayXd y_all = dat.y_s.array();

		Eigen::MatrixXd log_terms_mat(dat.G, n_nodes);
		std::vector<Eigen::VectorXd> mu_all_k_vec(n_nodes);
		for (int k = 0; k < n_nodes; ++k) {
			const Eigen::ArrayXd eta_all_k = eta_all.array() + b_vals[k];
			mu_all_k_vec[k] = clamp_eta_pg(eta_all_k).exp().matrix();
			const Eigen::ArrayXd term_all_k = y_all * eta_all_k - mu_all_k_vec[k].array() - dat.log_fact_y.array();
			for (int gi = 0; gi < dat.G; ++gi) {
				log_terms_mat(gi, k) = dat.gh.log_norm_weights[k] + 
				                       term_all_k.segment(dat.grp_start[gi], dat.grp_size[gi]).sum();
			}
		}

		Eigen::VectorXd ll_g_vec(dat.G);
		for (int gi = 0; gi < dat.G; ++gi) ll_g_vec[gi] = log_sum_exp_p(log_terms_mat.row(gi));

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total, total);
		H(dat.p, dat.p) = soft_barrier_hessian(log_sigma);

		Eigen::MatrixXd E_Hik_sum = Eigen::MatrixXd::Zero(total, total);
		Eigen::MatrixXd E_GiGiT_sum = Eigen::MatrixXd::Zero(total, total);
		Eigen::MatrixXd G_avg_outer_sum = Eigen::MatrixXd::Zero(total, total);

		for (int k = 0; k < n_nodes; k++) {
			Eigen::VectorXd pk_vec(dat.G);
			for (int gi = 0; gi < dat.G; gi++) pk_vec[gi] = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);

			Eigen::VectorXd pk_expanded(dat.n);
			for (int gi = 0; gi < dat.G; gi++) pk_expanded.segment(dat.grp_start[gi], dat.grp_size[gi]).setConstant(pk_vec[gi]);

			// 1. E[H_ik]
			// Beta-Beta block: sum_gi pk_gi * (-Xg^T * diag(mu_k_gi) * Xg) = -X_s^T * diag(pk_expanded * mu_all_k) * X_s
			Eigen::VectorXd w_Hik_beta = pk_expanded.cwiseProduct(mu_all_k_vec[k]);
			E_Hik_sum.topLeftCorner(dat.p, dat.p).noalias() -= weighted_crossprod(dat.X_s, w_Hik_beta);
			// Sigma-Sigma block: sum_gi pk_gi * H_ik(p,p)
			const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
			for (int gi = 0; gi < dat.G; ++gi) {
				const double pk = pk_vec[gi];
				if (pk < 1e-15) continue;
				const int start = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				const double sum_mu_k_gi = mu_all_k_vec[k].segment(start, sz).sum();
				const double sum_res_k_gi = (y_all.array().segment(start, sz) - mu_all_k_vec[k].array().segment(start, sz)).sum();
				E_Hik_sum(dat.p, dat.p) += pk * ((-sum_mu_k_gi * node_factor * node_factor * sigma + sum_res_k_gi * node_factor) * sigma);
				
				// Beta-Sigma block
				Eigen::VectorXd d2L_db_dlogsigma_k_gi = -(dat.X_s.middleRows(start, sz).transpose() * mu_all_k_vec[k].segment(start, sz)) * (node_factor * sigma);
				E_Hik_sum.block(0, dat.p, dat.p, 1).noalias() += pk * d2L_db_dlogsigma_k_gi;
			}
		}

		for (int gi = 0; gi < dat.G; ++gi) {
			Eigen::VectorXd G_avg_gi = Eigen::VectorXd::Zero(total);
			Eigen::MatrixXd E_GiGiT_gi = Eigen::MatrixXd::Zero(total, total);
			for (int k = 0; k < n_nodes; k++) {
				const double pk = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
				if (pk < 1e-15) continue;
				
				const int start = dat.grp_start[gi];
				const int sz = dat.grp_size[gi];
				Eigen::VectorXd G_ik = Eigen::VectorXd::Zero(total);
				Eigen::VectorXd res_k_gi = (y_all.array().segment(start, sz) - mu_all_k_vec[k].array().segment(start, sz)).matrix();
				G_ik.head(dat.p).noalias() = dat.X_s.middleRows(start, sz).transpose() * res_k_gi;
				const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
				G_ik[dat.p] = res_k_gi.sum() * node_factor * sigma;

				G_avg_gi.noalias() += pk * G_ik;
				E_GiGiT_gi.noalias() += pk * (G_ik * G_ik.transpose());
			}
			E_GiGiT_sum.noalias() += E_GiGiT_gi;
			G_avg_outer_sum.noalias() += G_avg_gi * G_avg_gi.transpose();
		}

		for (int r = 0; r < dat.p; r++) for (int c = 0; c < r; c++) E_Hik_sum(r, c) = E_Hik_sum(c, r);
		E_Hik_sum.block(dat.p, 0, 1, dat.p) = E_Hik_sum.block(0, dat.p, dat.p, 1).transpose();

		H -= (E_Hik_sum + E_GiGiT_sum - G_avg_outer_sum);
		return H;
	}
};

} // namespace

// [[Rcpp::export]]
SEXP fast_poisson_glmm_cpp(
	const NumericMatrix& X_r,
	const NumericVector& y_r,
	const IntegerVector& group_id_r,
	int j_T,
	Nullable<NumericVector> warm_start_params = R_NilValue,
	bool smart_cold_start = true,
	bool estimate_only = false,
	int n_gh = 20,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue
) {
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.rows(), X_r.cols());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());

	const int n = X.rows();
	const int p = X.cols();
	const int total = p + 1;

	std::vector<int> gid_v(n);
	for (int i = 0; i < n; ++i) gid_v[i] = group_id[i];

	PoissonGLMMData dat(X, y, gid_v, n_gh);

	// Initialize
	Eigen::VectorXd par(total);
	if (warm_start_params.isNotNull()) {
		NumericVector sp(warm_start_params);
		if (sp.size() == total) {
			for (int i = 0; i < total; ++i) par[i] = sp[i];
		}
	} else if (smart_cold_start) {
		// Init: beta via OLS on log(y+0.5), log_sigma = -3
		Eigen::VectorXd log_y_safe = (y.array() + 0.5).log().matrix();
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
		par.head(p) = cod.solve(log_y_safe);
		par[total - 1] = -3.0;
	} else {
		par.head(p).setZero();
		par[total - 1] = -3.0;
	}

	PoissonGLMMObjective obj(dat);
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
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

		const double pen = soft_barrier(par[total - 1]);
		const double true_neg_ll = neg_ll - pen;
		if (estimate_only) {
			return List::create(
				Named("b")          = par.head(p),
				Named("log_sigma")  = par[total - 1],
				Named("ssq_b_T")    = NA_REAL,
				Named("vcov")       = R_NilValue,
				Named("converged")  = converged,
				Named("neg_loglik") = true_neg_ll,
				Named("fisher_information") = R_NilValue
			);
		}

		Eigen::MatrixXd information = obj.hessian(par);
		information(total - 1, total - 1) -= soft_barrier_hessian(par[total - 1]);

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
		if (converged) {
		Eigen::MatrixXd H_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(H_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(H_free.rows(), H_free.cols()));
			vcov = expand_free_covariance(total, fixed_spec, inv_free, true);
			if (vcov.allFinite() && j_T < p) ssq_b_T = vcov(j_T, j_T);
		}
	}

	return List::create(
		Named("b")          = par.head(p),
		Named("log_sigma")  = par[total - 1],
		Named("ssq_b_T")    = ssq_b_T,
		Named("vcov")       = vcov,
		Named("converged")  = converged,
		Named("neg_loglik") = true_neg_ll,
		Named("fisher_information") = information
	);
}

// ── R-exported: score (gradient of log_lik) at arbitrary par ─────────────────
// [[Rcpp::export]]
SEXP get_poisson_glmm_score_cpp(
	const NumericMatrix& X_r,
	const NumericVector& y_r,
	const IntegerVector& group_id_r,
	const NumericVector& par_r,
	int n_gh = 20
) {
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.rows(), X_r.cols());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());
	Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

	std::vector<int> gid_v(group_id.size());
	for (int i = 0; i < group_id.size(); ++i) gid_v[i] = group_id[i];
	PoissonGLMMData dat(X, y, gid_v, n_gh);
	PoissonGLMMObjective obj(dat);
	Eigen::VectorXd grad;
	obj(par, grad);
	return wrap(-grad);
}

// ── R-exported: observed information (Hessian of neg_ll) at par ─────────────
// [[Rcpp::export]]
SEXP get_poisson_glmm_hessian_cpp(
	const NumericMatrix& X_r,
	const NumericVector& y_r,
	const IntegerVector& group_id_r,
	const NumericVector& par_r,
	int n_gh = 20
) {
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.rows(), X_r.cols());
	Eigen::Map<const Eigen::VectorXd> y(y_r.begin(), y_r.size());
	Eigen::Map<const Eigen::VectorXi> group_id(group_id_r.begin(), group_id_r.size());
	Eigen::Map<const Eigen::VectorXd> par(par_r.begin(), par_r.size());

	std::vector<int> gid_v(group_id.size());
	for (int i = 0; i < group_id.size(); ++i) gid_v[i] = group_id[i];
	PoissonGLMMData dat(X, y, gid_v, n_gh);
	PoissonGLMMObjective obj(dat);
	return wrap(-obj.hessian(par));
}
