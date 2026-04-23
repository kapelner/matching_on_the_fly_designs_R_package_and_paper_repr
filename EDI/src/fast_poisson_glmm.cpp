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

struct PoissonGLMMData {
	Eigen::MatrixXd X_s;
	Eigen::VectorXd y_s;
	Eigen::VectorXd log_fact_y;  // precomputed log(y!) = lgamma(y+1)
	std::vector<int> grp_start;
	std::vector<int> grp_size;
	int n, p, G;
	GHRuleP gh;

	PoissonGLMMData(
		const Eigen::MatrixXd& X,
		const Eigen::VectorXd& y,
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

	double operator()(const Eigen::VectorXd& par, Eigen::VectorXd& grad) {
		const double log_sigma = par[dat.p];
		const double sigma     = std::exp(log_sigma);
		const Eigen::VectorXd beta   = par.head(dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		double total_nll = soft_barrier(log_sigma);
		grad.setZero(dat.p + 1);
		const double center = 5.0, scale = 10.0;
		const double d = std::abs(log_sigma) - center;
		if (d > 0.0) grad[dat.p] += 2.0 * scale * d * (log_sigma > 0 ? 1.0 : -1.0);

		for (int gi = 0; gi < dat.G; ++gi) {
			const int start = dat.grp_start[gi];
			const int sz    = dat.grp_size[gi];
			const Eigen::MatrixXd Xg   = dat.X_s.middleRows(start, sz);
			const Eigen::VectorXd eta0 = Xg * beta;

			Eigen::VectorXd log_terms(n_nodes);
			// mu_nodes[k][r] = exp(eta0[r] + b_vals[k])
			std::vector<Eigen::VectorXd> mu_nodes(n_nodes);

			for (int k = 0; k < n_nodes; ++k) {
				double ll = dat.gh.log_norm_weights[k];
				mu_nodes[k].resize(sz);
				for (int r = 0; r < sz; ++r) {
					const double eta_rk = eta0[r] + b_vals[k];
					// Cap exp to avoid overflow; Poisson log-lik: y*eta - exp(eta) - log(y!)
					const double cap_eta = std::min(eta_rk, 700.0);
					const double mu_rk   = std::exp(cap_eta);
					mu_nodes[k][r] = mu_rk;
					ll += dat.y_s[start + r] * eta_rk - mu_rk - dat.log_fact_y[start + r];
				}
				log_terms[k] = ll;
			}

			const double ll_g = log_sum_exp_p(log_terms);
			if (!std::isfinite(ll_g)) { grad.setZero(dat.p + 1); return 1e100; }
			total_nll -= ll_g;

			for (int k = 0; k < n_nodes; ++k) {
				const double post_k = std::exp(log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;

				// residual = y - mu at this GH node
				Eigen::VectorXd res_k(sz);
				double res_sum = 0.0;
				for (int r = 0; r < sz; ++r) {
					res_k[r] = dat.y_s[start + r] - mu_nodes[k][r];
					res_sum += res_k[r];
				}

				// d(neg_ll)/d(beta)
				grad.head(dat.p) -= post_k * (Xg.transpose() * res_k);

				// d(neg_ll)/d(log_sigma): chain rule u_k = sqrt(2)*sigma*t_k
				// d ll_k / d log_sigma = sum_r (y_r - mu_rk) * u_k
				grad[dat.p] -= post_k * res_sum * b_vals[k];
			}
		}
		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
		return numerical_hessian(*this, par);
	}
};

} // namespace

// [[Rcpp::export]]
List fast_poisson_glmm_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Eigen::VectorXi& group_id,
	int j_T,
	bool estimate_only = false,
	int n_gh = 20,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
	std::string optimization_alg = "lbfgs"
) {
	const int n = X.rows();
	const int p = X.cols();
	const int total = p + 1;

	std::vector<int> gid_v(n);
	for (int i = 0; i < n; ++i) gid_v[i] = group_id[i];

	PoissonGLMMData dat(X, y, gid_v, n_gh);

	// Init: beta via OLS on log(y+0.5), log_sigma = -3
	Eigen::VectorXd par(total);
	{
		Eigen::VectorXd log_y_safe = (y.array() + 0.5).log().matrix();
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
		par.head(p) = cod.solve(log_y_safe);
	}
	par[total - 1] = -3.0;

	PoissonGLMMObjective obj(dat);
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
			Named("b")          = par.head(p),
			Named("log_sigma")  = par[total - 1],
			Named("ssq_b_T")    = NA_REAL,
			Named("converged")  = false,
			Named("neg_loglik") = NA_REAL
		);
	}

	const double pen = soft_barrier(par[total - 1]);
	const double true_neg_ll = neg_ll - pen;

	double ssq_b_T = NA_REAL;
	Eigen::MatrixXd vcov = Eigen::MatrixXd::Constant(total, total, NA_REAL);
	if (!estimate_only && converged) {
		Eigen::MatrixXd H = obj.hessian(par);
		Eigen::MatrixXd H_free = subset_matrix(H, fixed_spec.free_idx, fixed_spec.free_idx);
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
		Named("neg_loglik") = true_neg_ll
	);
}
