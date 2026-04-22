// Ordinal Cumulative-Logit GLMM for KK designs via Gauss-Hermite quadrature.
//
// Model:  logit P(Y_ij <= k | u_i) = alpha_k - X_ij' beta - u_i
//   u_i ~ N(0, sigma^2)    (random intercept per matched pair / singleton)
//   y_ij in {1, ..., K}    (1-indexed integer ordinal outcome)
//   X has NO intercept column (cutpoints alpha_1..alpha_{K-1} serve as intercepts)
//
// Parameter vector: par = [alpha_1, log_d2, ..., log_d_{K-1}, beta_0..beta_{p-1}, log_sigma]
//   where alpha_k = alpha_1 + sum_{j=2}^{k} exp(log_dj)   (ensures monotonicity)
//   Total length: (K-1) + p + 1 = K + p - 1 + 1 = K + p

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

GHRule gauss_hermite_rule_ord(int n) {
	Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n);
	for (int i = 0; i < n - 1; ++i) {
		const double v = std::sqrt((i + 1.0) / 2.0);
		J(i, i + 1) = v;
		J(i + 1, i) = v;
	}
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(J);
	GHRule rule;
	rule.nodes = es.eigenvalues();
	rule.log_norm_weights = (std::sqrt(M_PI) * es.eigenvectors().row(0).array().square()).log() - 0.5 * std::log(M_PI);
	return rule;
}

inline double log_sum_exp_v(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

// Recover the k-th cutpoint (1-indexed k in 1..K-1) from par[0..K-2]
inline double get_alpha(const Eigen::VectorXd& par, int k) {
	// k in 1..(K-1); par[0]=alpha_1; par[j]=log_diff_{j+1} for j>=1
	double a = par[0];
	for (int j = 1; j < k; ++j) a += std::exp(par[j]);
	return a;
}

// log P(Y=y_k | eta) under cumulative logit, y_k in 1..K, K levels
inline double cumlogit_log_prob(
	int y_k, int K, const double* alpha, double eta
) {
	// alpha is a precomputed array alpha[0..K-2], indexed 0-based (alpha[0]=alpha_1...)
	double p_upper = (y_k >= K) ? 1.0 : plogis_safe(alpha[y_k - 1] - eta);
	double p_lower = (y_k <= 1)  ? 0.0 : plogis_safe(alpha[y_k - 2] - eta);
	return std::log(std::max(1e-15, p_upper - p_lower));
}

struct OrdinalGLMMData {
	// Sorted so observations within each group are contiguous
	Eigen::MatrixXd X_s;       // n x p (no intercept)
	std::vector<int> y_s;      // 1-indexed ordinal outcomes, length n
	std::vector<int> grp_start;// start index of each group in sorted arrays
	std::vector<int> grp_size; // number of obs in each group
	int n, p, G, K;
	GHRule gh;
	double max_abs_log_sigma;

	OrdinalGLMMData(
		const Eigen::MatrixXd& X,
		const std::vector<int>& y_int_r, // 1-indexed
		const std::vector<int>& group_id_r, // 1-indexed
		int K_,
		int n_gh,
		double max_abs_log_sigma_
	) : n(X.rows()), p(X.cols()), K(K_),
	    gh(gauss_hermite_rule_ord(n_gh)),
	    max_abs_log_sigma(max_abs_log_sigma_) {

		// Sort by group_id
		std::vector<int> ord(n);
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(),
			[&](int a, int b){ return group_id_r[a] < group_id_r[b]; });

		X_s.resize(n, p);
		y_s.resize(n);
		for (int i = 0; i < n; ++i) {
			X_s.row(i) = X.row(ord[i]);
			y_s[i] = y_int_r[ord[i]];
		}

		// Build group structure
		int prev = -1;
		for (int i = 0; i < n; ++i) {
			int g = group_id_r[ord[i]];
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

class OrdinalGLMMObjective {
	const OrdinalGLMMData& dat;

public:
	explicit OrdinalGLMMObjective(const OrdinalGLMMData& d) : dat(d) {}

	double value(const Eigen::VectorXd& par) const {
		const int n_alpha = dat.K - 1;
		const double log_sigma = par[n_alpha + dat.p];
		if (!std::isfinite(log_sigma) || std::abs(log_sigma) > dat.max_abs_log_sigma)
			return 1e100;
		const double sigma = std::exp(log_sigma);

		// Precompute cutpoints
		std::vector<double> alpha(n_alpha);
		alpha[0] = par[0];
		for (int k = 1; k < n_alpha; ++k) alpha[k] = alpha[k-1] + std::exp(par[k]);

		const Eigen::VectorXd beta = par.segment(n_alpha, dat.p);
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
					ll += cumlogit_log_prob(dat.y_s[start + r], dat.K, alpha.data(), eta);
				}
				log_terms[k] = ll;
			}
			const double ll_g = log_sum_exp_v(log_terms);
			if (!std::isfinite(ll_g)) return 1e100;
			total_ll += ll_g;
		}
		return -total_ll;
	}

	double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
		const int n_alpha = dat.K - 1;
		const double log_sigma = par[n_alpha + dat.p];
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.segment(n_alpha, dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		// Precompute cutpoints and their gradients w.r.t parameterization
		std::vector<double> alpha(n_alpha);
		alpha[0] = par[0];
		for (int k = 1; k < n_alpha; ++k) alpha[k] = alpha[k-1] + std::exp(par[k]);

		double total_nll = 0.0;
		grad.setZero();

		for (int gi = 0; gi < dat.G; ++gi) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::MatrixXd Xg = dat.X_s.middleRows(start, sz);
			const Eigen::VectorXd eta0 = Xg * beta;

			Eigen::VectorXd log_terms(n_nodes);
			// For each node, store d(log P_ij)/d(alpha_k) and d(log P_ij)/d(eta_ij)
			// But simpler: just dL_i/d(theta) = sum_k post_ik * sum_j d(log P_ijk)/d(theta)
			
			std::vector<std::vector<Eigen::VectorXd>> dlogp_dalpha_nodes(n_nodes, std::vector<Eigen::VectorXd>(sz, Eigen::VectorXd::Zero(n_alpha)));
			std::vector<Eigen::VectorXd> dlogp_deta_nodes(n_nodes, Eigen::VectorXd::Zero(sz));

			for (int k = 0; k < n_nodes; ++k) {
				double ll = dat.gh.log_norm_weights[k];
				for (int r = 0; r < sz; ++r) {
					int y_ir = dat.y_s[start + r];
					double eta_ijk = eta0[r] + b_vals[k];
					
					double F_up = (y_ir >= dat.K) ? 1.0 : plogis_safe(alpha[y_ir - 1] - eta_ijk);
					double F_lo = (y_ir <= 1)      ? 0.0 : plogis_safe(alpha[y_ir - 2] - eta_ijk);
					double prob = std::max(1e-15, F_up - F_lo);
					ll += std::log(prob);

					// d(log P)/d(cutpoint)
					double f_up = (y_ir >= dat.K) ? 0.0 : F_up * (1.0 - F_up);
					double f_lo = (y_ir <= 1)      ? 0.0 : F_lo * (1.0 - F_lo);
					
					if (y_ir < dat.K) dlogp_dalpha_nodes[k][r][y_ir - 1] += f_up / prob;
					if (y_ir > 1)      dlogp_dalpha_nodes[k][r][y_ir - 2] -= f_lo / prob;

					// d(log P)/d(eta) = - [f_up - f_lo] / prob
					dlogp_deta_nodes[k][r] = -(f_up - f_lo) / prob;
				}
				log_terms[k] = ll;
			}

			const double ll_g = log_sum_exp_v(log_terms);
			total_nll -= ll_g;

			for (int k = 0; k < n_nodes; ++k) {
				double post_k = std::exp(log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;

				Eigen::VectorXd dLi_dalpha = Eigen::VectorXd::Zero(n_alpha);
				double dLi_deta_sum = 0.0;
				Eigen::VectorXd dLi_dbeta = Eigen::VectorXd::Zero(dat.p);

				for (int r = 0; r < sz; ++r) {
					dLi_dalpha += dlogp_dalpha_nodes[k][r];
					double de = dlogp_deta_nodes[k][r];
					dLi_deta_sum += de;
					dLi_dbeta += de * Xg.row(r).transpose();
				}

				// dLL/dalpha_k (preliminary, still need chain rule for log-diff parameterization)
				for (int j = 0; j < n_alpha; ++j) {
					grad[j] -= post_k * dLi_dalpha[j];
				}

				// dLL/dbeta
				grad.segment(n_alpha, dat.p) -= post_k * dLi_dbeta;

				// dLL/dlog_sigma
				double dll_dsigma = dLi_deta_sum * std::sqrt(2.0) * dat.gh.nodes[k];
				grad[n_alpha + dat.p] -= post_k * dll_dsigma * sigma;
			}
		}

		// Apply chain rule for cutpoint parameterization: 
		// par[0]=alpha_1, par[j]=log(alpha_{j+1}-alpha_j)
		// dLL/d par[j] = sum_{k=j+1}^{n_alpha} (dLL/d alpha_k) * (d alpha_k / d par[j])
		// For j=0: d alpha_k / d par[0] = 1 for all k
		// For j>0: d alpha_k / d par[j] = exp(par[j]) if k >= j+1, else 0
		
		Eigen::VectorXd dLL_dalpha = grad.head(n_alpha);
		grad[0] = dLL_dalpha.sum();
		for (int j = 1; j < n_alpha; ++j) {
			double sum_higher = 0.0;
			for (int k = j; k < n_alpha; ++k) sum_higher += dLL_dalpha[k];
			grad[j] = sum_higher * std::exp(par[j]);
		}

		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
		return numerical_hessian(*this, par);
	}
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_ordinal_glmm_score_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXi& y_int,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int K,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	std::vector<int> y_v(y_int.size()), gid_v(group_id.size());
	for (int i = 0; i < y_int.size(); ++i) { y_v[i] = y_int[i]; gid_v[i] = group_id[i]; }

	OrdinalGLMMData dat(X, y_v, gid_v, K, n_gh, max_abs_log_sigma);
	OrdinalGLMMObjective obj(dat);

	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_ordinal_glmm_hessian_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXi& y_int,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int K,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	std::vector<int> y_v(y_int.size()), gid_v(group_id.size());
	for (int i = 0; i < y_int.size(); ++i) { y_v[i] = y_int[i]; gid_v[i] = group_id[i]; }

	OrdinalGLMMData dat(X, y_v, gid_v, K, n_gh, max_abs_log_sigma);
	OrdinalGLMMObjective obj(dat);

	return -obj.hessian(params);
}

// [[Rcpp::export]]
List fast_ordinal_glmm_cpp(
	const Eigen::MatrixXd& X,     // n x p, NO intercept column; treatment at column j_T (0-based)
	const Eigen::VectorXi& y_int, // 1-indexed ordinal outcomes, length n
	const Eigen::VectorXi& group_id, // group IDs, length n (sorted internally)
	int K,                        // number of ordinal levels
	int j_T,                      // 0-based treatment column index in X
	bool estimate_only = false,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::NumericVector> start = R_NilValue,  // optional warm start
	std::string optimization_alg = "lbfgs"
) {
	const int n = X.rows();
	const int p = X.cols();
	const int n_alpha = K - 1;
	const int total = n_alpha + p + 1; // cutpoint params + betas + log_sigma

	// Convert Eigen/R vectors to std::vector for OrdinalGLMMData
	std::vector<int> y_v(n), gid_v(n);
	for (int i = 0; i < n; ++i) { y_v[i] = y_int[i]; gid_v[i] = group_id[i]; }

	OrdinalGLMMData dat(X, y_v, gid_v, K, n_gh, max_abs_log_sigma);

	// Initialize parameters
	Eigen::VectorXd par(total);
	if (start.isNotNull()) {
		Rcpp::NumericVector sv(start);
		if (sv.size() == total) {
			for (int i = 0; i < total; ++i) par[i] = sv[i];
		} else {
			// Fallback: zero-initialize
			par.head(n_alpha).setZero();
			par.segment(n_alpha, p).setZero();
			par[total - 1] = -3.0;
		}
	} else {
		// Cutpoints: alpha_1 = 0, log_diffs = 0 (evenly spaced by 1), log_sigma = -3
		par.head(n_alpha).setZero();
		par.segment(n_alpha, p).setZero();
		par[total - 1] = -3.0;
	}

	OrdinalGLMMObjective obj(dat);

	double neg_ll = NA_REAL;
	int niter = maxit;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_likelihood(obj, par, maxit, eps_g, optimization_alg, "lbfgs");
		par = fit.params;
		neg_ll = fit.value;
		niter = fit.niter;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("b")          = par.segment(n_alpha, p),
			Named("alpha")      = par.head(n_alpha),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	// Treatment coefficient: 0-based index j_T in X → index (n_alpha + j_T) in par
	const int j_T_full = n_alpha + j_T;

	double ssq_b_T = NA_REAL;
	if (!estimate_only && converged) {
		Eigen::MatrixXd H = obj.hessian(par);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv = ldlt.solve(Eigen::MatrixXd::Identity(total, total));
			if (inv.allFinite() && j_T_full < total) {
				ssq_b_T = inv(j_T_full, j_T_full);
			}
		}
	}

	return List::create(
		Named("b")          = par.segment(n_alpha, p),
		Named("alpha")      = par.head(n_alpha),
		Named("log_sigma")  = par[total - 1],
		Named("ssq_b_T")    = ssq_b_T,
		Named("converged")  = converged,
		Named("neg_loglik") = neg_ll
	);
}
