// Ordinal Cumulative-Link GLMM for KK designs via Gauss-Hermite quadrature.
// Supports link functions: logit, probit, cauchit, cloglog.
//
// Model:  F^{-1}(P(Y_ij <= k | u_i)) = alpha_k - X_ij' beta - u_i
//   u_i ~ N(0, sigma^2)    (random intercept per matched pair / singleton)
//   y_ij in {1, ..., K}    (1-indexed integer ordinal outcome)
//   X has NO intercept column
//
// Parameter vector: par = [alpha_1, log_d2, ..., log_d_{K-1}, beta_0..beta_{p-1}, log_sigma]
//   alpha_k = alpha_1 + sum_{j=2}^{k} exp(log_dj)  (ensures monotonicity)
//   Total length: (K-1) + p + 1

#include "_helper_functions.h"
#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>

using namespace Rcpp;

namespace {

// ---- Link functions ----

enum class LinkType { LOGIT, PROBIT, CAUCHIT, CLOGLOG };

inline LinkType parse_link(const std::string& s) {
	if (s == "logit")   return LinkType::LOGIT;
	if (s == "probit")  return LinkType::PROBIT;
	if (s == "cauchit") return LinkType::CAUCHIT;
	if (s == "cloglog") return LinkType::CLOGLOG;
	Rcpp::stop("fast_ordinal_clmm_cpp: unknown link '%s'", s.c_str());
}

inline double link_cdf(double x, LinkType lt) {
	switch (lt) {
		case LinkType::LOGIT:   return plogis_safe(x);
		case LinkType::PROBIT:  return R::pnorm(x, 0.0, 1.0, 1, 0);
		case LinkType::CAUCHIT: return 0.5 + std::atan(x) / M_PI;
		case LinkType::CLOGLOG:
			if (x >  100.0) return 1.0;
			if (x < -700.0) return 0.0;
			return 1.0 - std::exp(-std::exp(x));
	}
	return 0.0;
}

inline double link_pdf(double x, LinkType lt) {
	switch (lt) {
		case LinkType::LOGIT: { double p = plogis_safe(x); return p * (1.0 - p); }
		case LinkType::PROBIT:  return R::dnorm(x, 0.0, 1.0, 0);
		case LinkType::CAUCHIT: return 1.0 / (M_PI * (1.0 + x * x));
		case LinkType::CLOGLOG: {
			if (x > 700.0 || x < -700.0) return 0.0;
			const double ex = std::exp(x);
			return ex * std::exp(-ex);
		}
	}
	return 0.0;
}

// ---- Gauss-Hermite quadrature ----

struct GHRule {
	Eigen::VectorXd nodes;
	Eigen::VectorXd log_norm_weights;
};

GHRule gauss_hermite_rule_clmm(int n) {
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

inline double log_sum_exp_clmm(const Eigen::VectorXd& x) {
	const double m = x.maxCoeff();
	if (!std::isfinite(m)) return m;
	return m + std::log((x.array() - m).exp().sum());
}

// ---- Data ----

struct OrdinalCLMMData {
	Eigen::MatrixXd X_s;
	std::vector<int> y_s;
	std::vector<int> grp_start, grp_size;
	int n, p, G, K;
	GHRule gh;
	double max_abs_log_sigma;
	LinkType link;

	OrdinalCLMMData(
		const Eigen::MatrixXd& X,
		const std::vector<int>& y_r,
		const std::vector<int>& gid_r,
		int K_, int n_gh,
		double max_abs_log_sigma_,
		LinkType link_
	) : n(X.rows()), p(X.cols()), K(K_),
	    gh(gauss_hermite_rule_clmm(n_gh)),
	    max_abs_log_sigma(max_abs_log_sigma_),
	    link(link_) {

		std::vector<int> ord(n);
		std::iota(ord.begin(), ord.end(), 0);
		std::stable_sort(ord.begin(), ord.end(),
			[&](int a, int b){ return gid_r[a] < gid_r[b]; });

		X_s.resize(n, p);
		y_s.resize(n);
		for (int i = 0; i < n; ++i) {
			X_s.row(i) = X.row(ord[i]);
			y_s[i] = y_r[ord[i]];
		}

		int prev = -1;
		for (int i = 0; i < n; ++i) {
			int g = gid_r[ord[i]];
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

// Soft barrier on log_sigma: penalizes |log_sigma| > center smoothly.
// Used in both value() and operator() so that L-BFGS never escapes the barrier.
inline double clmm_sigma_penalty(double log_sigma, double center = 5.0, double scale = 10.0) {
	const double d = std::abs(log_sigma) - center;
	if (d <= 0.0) return 0.0;
	return scale * d * d;
}

// ---- Objective ----

class OrdinalCLMMObjective {
	const OrdinalCLMMData& dat;

public:
	explicit OrdinalCLMMObjective(const OrdinalCLMMData& d) : dat(d) {}

	double value(const Eigen::VectorXd& par) const {
		const int na = dat.K - 1;
		const double log_sigma = par[na + dat.p];
		if (!std::isfinite(log_sigma) || std::abs(log_sigma) > dat.max_abs_log_sigma)
			return 1e100;
		const double pen = clmm_sigma_penalty(log_sigma);
		const double sigma = std::exp(log_sigma);

		std::vector<double> alpha(na);
		alpha[0] = par[0];
		for (int k = 1; k < na; ++k) alpha[k] = alpha[k-1] + std::exp(par[k]);

		const Eigen::VectorXd beta   = par.segment(na, dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int nn = (int)b_vals.size();

		double total_ll = 0.0;
		for (int gi = 0; gi < dat.G; ++gi) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::VectorXd eta0 = dat.X_s.middleRows(start, sz) * beta;

			Eigen::VectorXd log_terms(nn);
			for (int k = 0; k < nn; ++k) {
				double ll = dat.gh.log_norm_weights[k];
				for (int r = 0; r < sz; ++r) {
					const double eta = eta0[r] + b_vals[k];
					const int y_ir   = dat.y_s[start + r];
					double F_up = (y_ir >= dat.K) ? 1.0 : link_cdf(alpha[y_ir - 1] - eta, dat.link);
					double F_lo = (y_ir <= 1)     ? 0.0 : link_cdf(alpha[y_ir - 2] - eta, dat.link);
					ll += std::log(std::max(1e-15, F_up - F_lo));
				}
				log_terms[k] = ll;
			}
			const double ll_g = log_sum_exp_clmm(log_terms);
			if (!std::isfinite(ll_g)) return 1e100;
			total_ll += ll_g;
		}
		return -total_ll + pen;
	}

	// Analytical gradient. grad is accumulated as negative log-likelihood gradient.
	double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
		const int na = dat.K - 1;
		const double log_sigma = par[na + dat.p];
		const double sigma     = std::exp(log_sigma);

		const double pen = clmm_sigma_penalty(log_sigma);

		std::vector<double> alpha(na);
		alpha[0] = par[0];
		for (int k = 1; k < na; ++k) alpha[k] = alpha[k-1] + std::exp(par[k]);

		const Eigen::VectorXd beta   = par.segment(na, dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int nn = (int)b_vals.size();

		double total_nll = pen;
		grad.setZero();
		// Gradient of soft barrier on log_sigma
		{
			const double center = 5.0, scale = 10.0;
			const double d = std::abs(log_sigma) - center;
			if (d > 0.0)
				grad[na + dat.p] += 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);
		}

		for (int gi = 0; gi < dat.G; ++gi) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::MatrixXd Xg   = dat.X_s.middleRows(start, sz);
			const Eigen::VectorXd eta0 = Xg * beta;

			Eigen::VectorXd log_terms(nn);
			// per-node per-obs gradients w.r.t. alpha and eta
			std::vector<std::vector<Eigen::VectorXd>> dlogp_da(nn,
				std::vector<Eigen::VectorXd>(sz, Eigen::VectorXd::Zero(na)));
			std::vector<Eigen::VectorXd> dlogp_deta(nn, Eigen::VectorXd::Zero(sz));

			for (int k = 0; k < nn; ++k) {
				double ll = dat.gh.log_norm_weights[k];
				for (int r = 0; r < sz; ++r) {
					const int y_ir     = dat.y_s[start + r];
					const double eta   = eta0[r] + b_vals[k];
					double F_up = (y_ir >= dat.K) ? 1.0 : link_cdf(alpha[y_ir - 1] - eta, dat.link);
					double F_lo = (y_ir <= 1)     ? 0.0 : link_cdf(alpha[y_ir - 2] - eta, dat.link);
					double prob = std::max(1e-15, F_up - F_lo);
					ll += std::log(prob);

					double f_up = (y_ir >= dat.K) ? 0.0 : link_pdf(alpha[y_ir - 1] - eta, dat.link);
					double f_lo = (y_ir <= 1)     ? 0.0 : link_pdf(alpha[y_ir - 2] - eta, dat.link);

					if (y_ir < dat.K) dlogp_da[k][r][y_ir - 1] +=  f_up / prob;
					if (y_ir > 1)     dlogp_da[k][r][y_ir - 2] -= f_lo / prob;

					// d(log P)/d(eta) = -(f_up - f_lo) / prob
					dlogp_deta[k][r] = -(f_up - f_lo) / prob;
				}
				log_terms[k] = ll;
			}

			const double ll_g = log_sum_exp_clmm(log_terms);
			total_nll -= ll_g;

			for (int k = 0; k < nn; ++k) {
				double post_k = std::exp(log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;

				Eigen::VectorXd dLi_da   = Eigen::VectorXd::Zero(na);
				Eigen::VectorXd dLi_dbeta = Eigen::VectorXd::Zero(dat.p);
				double dLi_deta_sum = 0.0;

				for (int r = 0; r < sz; ++r) {
					dLi_da   += dlogp_da[k][r];
					const double de = dlogp_deta[k][r];
					dLi_deta_sum   += de;
					dLi_dbeta      += de * Xg.row(r).transpose();
				}

				for (int j = 0; j < na; ++j)
					grad[j] -= post_k * dLi_da[j];

				grad.segment(na, dat.p) -= post_k * dLi_dbeta;

				// dLL/d(log_sigma) = dLL/d(sigma) * sigma
				double dll_dsigma = dLi_deta_sum * std::sqrt(2.0) * dat.gh.nodes[k];
				grad[na + dat.p] -= post_k * dll_dsigma * sigma;
			}
		}

		// Chain rule: log-diff parameterization
		// par[0]=alpha_1; par[j]=log(alpha_{j+1}-alpha_j) for j>=1
		// d/d par[0] = sum_k dLL/d alpha_k
		// d/d par[j] = exp(par[j]) * sum_{k>=j} dLL/d alpha_k   (j >= 1)
		const Eigen::VectorXd dLL_da = grad.head(na);
		grad[0] = dLL_da.sum();
		for (int j = 1; j < na; ++j) {
			double s = 0.0;
			for (int k = j; k < na; ++k) s += dLL_da[k];
			grad[j] = s * std::exp(par[j]);
		}

		return total_nll;
	}

	// Hessian via finite differences of analytical gradient (more accurate than
	// double finite-differencing the objective value).
	Eigen::MatrixXd hessian(const Eigen::VectorXd& par) const {
		const int n = par.size();
		Eigen::MatrixXd H(n, n);
		OrdinalCLMMObjective* self = const_cast<OrdinalCLMMObjective*>(this);
		for (int i = 0; i < n; ++i) {
			double h = 1e-4 * std::max(1.0, std::abs(par[i]));
			Eigen::VectorXd pp = par, pm = par;
			pp[i] += h; pm[i] -= h;
			Eigen::VectorXd gp(n), gm(n);
			(*self)(pp, gp);
			(*self)(pm, gm);
			H.col(i) = (gp - gm) / (2.0 * h);
		}
		return (H + H.transpose()) / 2.0;
	}
};

} // namespace

// [[Rcpp::export]]
List fast_ordinal_clmm_cpp(
	const Eigen::MatrixXd& X,         // n x p, NO intercept; treatment at col j_T (0-based)
	const Eigen::VectorXi& y_int,     // 1-indexed ordinal outcomes
	const Eigen::VectorXi& group_id,  // group IDs (sorted internally)
	int K,                            // number of ordinal levels
	int j_T,                          // 0-based treatment column in X
	std::string link = "logit",       // link: "logit"|"probit"|"cauchit"|"cloglog"
	bool estimate_only = false,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::NumericVector> start = R_NilValue,
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
	const int n  = X.rows();
	const int p  = X.cols();
	const int na = K - 1;
	const int total = na + p + 1;

	const LinkType lt = parse_link(link);

	std::vector<int> y_v(n), gid_v(n);
	for (int i = 0; i < n; ++i) { y_v[i] = y_int[i]; gid_v[i] = group_id[i]; }

	OrdinalCLMMData dat(X, y_v, gid_v, K, n_gh, max_abs_log_sigma, lt);

	Eigen::VectorXd par(total);
	if (start.isNotNull()) {
		Rcpp::NumericVector sv(start);
		if (sv.size() == total) {
			for (int i = 0; i < total; ++i) par[i] = sv[i];
		} else {
			par.head(na).setZero();
			par.segment(na, p).setZero();
			par[total - 1] = -3.0;
		}
	} else {
		par.head(na).setZero();
		par.segment(na, p).setZero();
		par[total - 1] = -3.0;
	}

	OrdinalCLMMObjective obj(dat);
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	double neg_ll = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs");
		par       = fit.params;
		neg_ll    = fit.value;
		converged = std::isfinite(neg_ll) && fit.converged;
	} catch (...) {
		return List::create(
			Named("b")          = par.segment(na, p),
			Named("alpha")      = par.head(na),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	const int j_T_full = na + j_T;

	double ssq_b_T = NA_REAL;
	if (!estimate_only && converged) {
		FixedParameterFunctor<OrdinalCLMMObjective> fixed_obj(obj, fixed_spec, par);
		Eigen::VectorXd params_free = subset_vector(par, fixed_spec.free_idx);
		Eigen::MatrixXd H_free = fixed_obj.hessian(params_free);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(H_free);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(H_free.rows(), H_free.cols()));
			if (inv_free.allFinite()) {
				Eigen::MatrixXd inv = expand_free_covariance(total, fixed_spec, inv_free, true);
				if (j_T_full < total) {
					ssq_b_T = inv(j_T_full, j_T_full);
				}
			}
		}
	}

	return List::create(
		Named("b")          = par.segment(na, p),
		Named("alpha")      = par.head(na),
		Named("log_sigma")  = par[total - 1],
		Named("ssq_b_T")    = ssq_b_T,
		Named("converged")  = converged,
		Named("neg_loglik") = neg_ll
	);
}
