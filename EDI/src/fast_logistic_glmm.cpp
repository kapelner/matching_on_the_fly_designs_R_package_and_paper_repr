// Logistic GLMM for KK designs via Gauss-Hermite quadrature.
//
// Model:  logit P(Y_ij = 1 | u_i) = X_ij' beta + u_i
//   u_i ~ N(0, sigma^2)    (random intercept per matched pair / singleton)
//   y_ij in [0,1]          (binary 0/1 or continuous proportion)
//   X includes an intercept column
//
// Parameter vector: par = [beta_0, beta_1(treatment), ..., beta_{p-1}, log_sigma]
//   Total length: p + 1
//
// log_sigma is parameterized with a soft-barrier penalty rather than a hard cut.
// This avoids the infinite-gradient issue at the boundary.

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace Rcpp;

namespace {

struct GHRule {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRule gauss_hermite_rule_log(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRule rule;
	rule.nodes = es.eigenvalues();
	rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log()
		                    - 0.5 * std::log(M_PI);
	return rule;
}

inline double log_sum_exp_v(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

inline double log1pexp_s(double x) {
	if (x > 0.0) return x + std::log1p(std::exp(-x));
	return std::log1p(std::exp(x));
}

// Soft barrier: adds a smooth penalty when |log_sigma| is large,
// avoiding the hard-wall numerical catastrophe.
// penalty(t) = 0 when |t| <= center, smoothly increases beyond.
inline double log_sigma_penalty(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return scale * d * d;
}

inline double log_sigma_penalty_grad(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
}

inline double log_sigma_penalty_hessian(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return 2.0 * scale;
}

struct LogisticGLMMData {
	Eigen::MatrixXd X_s;        // n x p (includes intercept)
	std::vector<double> y_s;    // responses in [0,1], length n
	std::vector<int> grp_start; // start index of each group
	std::vector<int> grp_size;  // size of each group
	int n, p, G;
	GHRule gh;

	LogisticGLMMData(
		const Eigen::MatrixXd& X,
		const std::vector<double>& y,
		const std::vector<int>& group_id,
		int n_gh
	) : n(X.rows()), p(X.cols()), gh(gauss_hermite_rule_log(n_gh)) {

		// Sort by group_id for contiguous group access
		std::vector<int> ord(n);
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(),
			[&](int a, int b){ return group_id[a] < group_id[b]; });

		X_s.resize(n, p);
		y_s.resize(n);
		for (int i = 0; i < n; ++i) {
			X_s.row(i) = X.row(ord[i]);
			y_s[i] = y[ord[i]];
		}

		// Build group structure
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

class LogisticGLMMObjective {
	const LogisticGLMMData& dat;

public:
	explicit LogisticGLMMObjective(const LogisticGLMMData& d) : dat(d) {}

	double value(const Eigen::VectorXd& par) const {
		const double log_sigma = par[dat.p];
		// Soft barrier for |log_sigma| > 5 (instead of hard cut)
		const double pen = log_sigma_penalty(log_sigma);
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		double total_ll = 0.0;
		for (int gi = 0; gi < dat.G; ++gi) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::VectorXd eta0 = dat.X_s.middleRows(start, sz) * beta;

			Eigen::VectorXd log_terms(n_nodes);
			for (int k = 0; k < n_nodes; ++k) {
				double ll = dat.gh.log_norm_weights[k];
				for (int r = 0; r < sz; ++r) {
					const double eta = eta0[r] + b_vals[k];
					ll += dat.y_s[start + r] * eta - log1pexp_s(eta);
				}
				log_terms[k] = ll;
			}
			const double ll_g = log_sum_exp_v(log_terms);
			if (!std::isfinite(ll_g)) return 1e50;
			total_ll += ll_g;
		}
		return -total_ll + pen;
	}

	double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
		const double log_sigma = par[dat.p];
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		double total_nll = log_sigma_penalty(log_sigma);
		grad.setZero();
		// Gradient of penalty
		const double center = 5.0, scale = 10.0;
		const double d = std::abs(log_sigma) - center;
		if (d > 0.0) {
			grad[dat.p] += 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
		}

		for (int gi = 0; gi < dat.G; ++gi) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::MatrixXd Xg = dat.X_s.middleRows(start, sz);
			const Eigen::VectorXd eta0 = Xg * beta;
			const Eigen::Map<const Eigen::VectorXd> y_g(dat.y_s.data() + start, sz);

			Eigen::VectorXd log_terms(n_nodes);
			std::vector<Eigen::VectorXd> mu_nodes(n_nodes);

			for (int k = 0; k < n_nodes; ++k) {
				const Eigen::ArrayXd eta_k = eta0.array() + b_vals[k];
				mu_nodes[k] = plogis_array_safe(eta_k).matrix();
				log_terms[k] = dat.gh.log_norm_weights[k] +
				               (y_g.array() * eta_k - log1pexp_array_safe(eta_k)).sum();
			}

			const double ll_g = log_sum_exp_v(log_terms);
			total_nll -= ll_g;

			// Posterior weights for gradient: exp(log_terms[k] - ll_g)
			for (int k = 0; k < n_nodes; ++k) {
				double post_k = std::exp(log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;

				Eigen::VectorXd res_k = y_g - mu_nodes[k];

				// dLL/dbeta
				grad.head(dat.p) -= post_k * (Xg.transpose() * res_k);

				// dLL/dsigma * dsigma/dlog_sigma
				// dLL/dsigma = sum_j (y_ij - mu_ijk) * sqrt(2) * nodes[k]
				double dll_dsigma = res_k.sum() * std::sqrt(2.0) * dat.gh.nodes[k];
				grad[dat.p] -= post_k * dll_dsigma * sigma;
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
		H(dat.p, dat.p) = log_sigma_penalty_hessian(log_sigma);

		for (int gi = 0; gi < dat.G; gi++) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::MatrixXd Xg = dat.X_s.middleRows(start, sz);
			const Eigen::VectorXd eta0 = Xg * beta;
			const Eigen::Map<const Eigen::VectorXd> y_g(dat.y_s.data() + start, sz);

			Eigen::VectorXd log_terms(n_nodes);
			std::vector<Eigen::VectorXd> mu_nodes(n_nodes);
			for (int k = 0; k < n_nodes; k++) {
				const Eigen::ArrayXd eta_k = eta0.array() + b_vals[k];
				mu_nodes[k] = plogis_array_safe(eta_k).matrix();
				log_terms[k] = dat.gh.log_norm_weights[k] +
				               (y_g.array() * eta_k - log1pexp_array_safe(eta_k)).sum();
			}
			const double ll_g = log_sum_exp_v(log_terms);

			Eigen::VectorXd post_weights(n_nodes);
			for (int k = 0; k < n_nodes; k++) post_weights[k] = std::exp(log_terms[k] - ll_g);

			Eigen::MatrixXd E_Hik = Eigen::MatrixXd::Zero(total, total);
			Eigen::MatrixXd E_GiGiT = Eigen::MatrixXd::Zero(total, total);
			Eigen::VectorXd G_avg = Eigen::VectorXd::Zero(total);

			for (int k = 0; k < n_nodes; k++) {
				const double pk = post_weights[k];
				if (pk < 1e-15) continue;

				Eigen::VectorXd G_ik = Eigen::VectorXd::Zero(total);
				Eigen::MatrixXd H_ik = Eigen::MatrixXd::Zero(total, total);

				Eigen::VectorXd res_k = y_g - mu_nodes[k];
				Eigen::VectorXd w_k = (mu_nodes[k].array() * (1.0 - mu_nodes[k].array())).matrix();
				const double sum_res = res_k.sum();
				const double sum_w = w_k.sum();

				G_ik.head(dat.p) = Xg.transpose() * res_k;
				const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
				G_ik[dat.p] = sum_res * node_factor * sigma;

				H_ik.topLeftCorner(dat.p, dat.p).noalias() = -weighted_crossprod(Xg, w_k);
				H_ik(dat.p, dat.p) = (-sum_w * node_factor * node_factor * sigma + sum_res * node_factor) * sigma;
				
				Eigen::VectorXd d2L_db_dsigma = -Xg.transpose() * w_k * node_factor * sigma;
				H_ik.block(0, dat.p, dat.p, 1) = d2L_db_dsigma;
				H_ik.block(dat.p, 0, 1, dat.p) = d2L_db_dsigma.transpose();

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
Eigen::VectorXd get_logistic_glmm_score_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 20
) {
	std::vector<double> y_v(y.size());
	std::vector<int> gid_v(group_id.size());
	for (int i = 0; i < y.size(); ++i) { y_v[i] = y[i]; gid_v[i] = group_id[i]; }
	
	LogisticGLMMData dat(X, y_v, gid_v, n_gh);
	LogisticGLMMObjective obj(dat);
	
	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	grad[X.cols()] -= log_sigma_penalty_grad(params[X.cols()]);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_logistic_glmm_hessian_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 20
) {
	std::vector<double> y_v(y.size());
	std::vector<int> gid_v(group_id.size());
	for (int i = 0; i < y.size(); ++i) { y_v[i] = y[i]; gid_v[i] = group_id[i]; }
	
	LogisticGLMMData dat(X, y_v, gid_v, n_gh);
	LogisticGLMMObjective obj(dat);
	
	Eigen::MatrixXd information = obj.hessian(params);
	information(X.cols(), X.cols()) -= log_sigma_penalty_hessian(params[X.cols()]);
	return -information;
}

// [[Rcpp::export]]
double get_logistic_glmm_neg_loglik_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int n_gh = 20
) {
	std::vector<double> y_v(y.size());
	std::vector<int> gid_v(group_id.size());
	for (int i = 0; i < y.size(); ++i) { y_v[i] = y[i]; gid_v[i] = group_id[i]; }

	LogisticGLMMData dat(X, y_v, gid_v, n_gh);
	LogisticGLMMObjective obj(dat);
	return likelihood_value(obj, params) - log_sigma_penalty(params[X.cols()]);
}

// [[Rcpp::export]]
List fast_logistic_glmm_cpp(
	const Eigen::MatrixXd& X,       // n x p, includes intercept; treatment at col j_T (0-based)
	const Eigen::VectorXd& y,       // responses in [0,1], length n
	const Eigen::VectorXi& group_id,// group IDs, sorted internally
	int j_T,                        // 0-based treatment column index in X
	bool estimate_only = false,
	int n_gh = 20,
	int maxit = 300,
	double eps_g = 1e-6,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
	const int n = X.rows();
	const int p = X.cols();
	const int total = p + 1; // betas + log_sigma

	std::vector<double> y_v(n);
	std::vector<int> gid_v(n);
	for (int i = 0; i < n; ++i) { y_v[i] = y[i]; gid_v[i] = group_id[i]; }

	LogisticGLMMData dat(X, y_v, gid_v, n_gh);

	// Initialize: betas = 0, log_sigma = -3 (sigma ≈ 0.05, small but away from boundary)
	Eigen::VectorXd par(total);
	par.head(p).setZero();
	par[total - 1] = -3.0;

	LogisticGLMMObjective obj(dat);
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

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
			Named("params")     = par,
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	// Remove soft-barrier penalty from neg_ll so it reflects the true neg log-likelihood
	const double pen = log_sigma_penalty(par[total - 1]);
	const double true_neg_ll = neg_ll - pen;
	Eigen::VectorXd score(total);
	obj(par, score);
	score[total - 1] -= log_sigma_penalty_grad(par[total - 1]);
	score = -score;
	Eigen::MatrixXd information = obj.hessian(par);
	information(total - 1, total - 1) -= log_sigma_penalty_hessian(par[total - 1]);

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
	if (!estimate_only && converged) {
		Eigen::MatrixXd information_free = subset_matrix(information, fixed_spec.free_idx, fixed_spec.free_idx);
		Eigen::MatrixXd cov_free = covariance_from_information(information_free);
		vcov = expand_free_covariance(total, fixed_spec, cov_free, true);
		if (j_T < p) {
			ssq_b_T = vcov(j_T, j_T);
		}
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
