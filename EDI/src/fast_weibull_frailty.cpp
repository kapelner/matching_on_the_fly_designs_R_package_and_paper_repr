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
	int m_n_groups;

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

		m_group_start.clear();
		m_group_end.clear();
		int i = 0;
		while (i < m_n) {
			const int g = sorted_gid[i];
			m_group_start.push_back(i);
			int j = i + 1;
			while (j < m_n && sorted_gid[j] == g) ++j;
			m_group_end.push_back(j);
			i = j;
		}
		m_n_groups = (int)m_group_start.size();
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
		m_n_groups(0)
	{
		build_group_structure(group_id);
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
		const Eigen::VectorXd eta = m_X * beta;
		const int n_nodes = (int)m_gh.nodes.size();

		Eigen::VectorXd u_vals(n_nodes);
		for (int k = 0; k < n_nodes; ++k)
			u_vals[k] = std::sqrt(2.0) * sigma_u * m_gh.nodes[k];

		double total_nll = 0.0;

		for (int g = 0; g < m_n_groups; ++g) {
			const int i_start = m_group_start[g];
			const int sz      = m_group_end[g] - i_start;
			const Eigen::MatrixXd Xg = m_X.middleRows(i_start, sz);

			Eigen::VectorXd log_terms(n_nodes);
			Eigen::MatrixXd w_mat(sz, n_nodes);
			Eigen::MatrixXd ew_mat(sz, n_nodes);

			for (int k = 0; k < n_nodes; ++k) {
				double ll = m_gh.log_norm_weights[k];
				for (int r = 0; r < sz; ++r) {
					const int i = i_start + r;
					double wik = (m_log_y[i] - eta[i] - u_vals[k]) / sigma_eps;
					if (wik > 700.0) wik = 700.0;
					const double ewik = std::exp(wik);
					w_mat(r, k)  = wik;
					ew_mat(r, k) = ewik;
					if (m_dead[i] > 0.5) {
						ll += wik - log_sigma_eps - m_log_y[i] - ewik;
					} else {
						ll -= ewik;
					}
				}
				log_terms[k] = ll;
			}

			const double ll_g = log_sum_exp_wf(log_terms);
			if (!std::isfinite(ll_g)) { grad.setZero(m_p + 2); return 1e100; }
			total_nll -= ll_g;

			for (int k = 0; k < n_nodes; ++k) {
				const double post_k = std::exp(log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;

				Eigen::VectorXd r_k(sz);
				double r_k_sum  = 0.0;
				double dll_dlse = 0.0;  // sum_i [w_ik*(exp(w_ik)-dead_i) - dead_i]

				for (int r = 0; r < sz; ++r) {
					const int i  = i_start + r;
					r_k[r]    = ew_mat(r, k) - m_dead[i];
					r_k_sum  += r_k[r];
					dll_dlse += w_mat(r, k) * r_k[r] - m_dead[i];
				}

				// d(neg_ll)/d(beta)
				grad.head(m_p) -= post_k * (Xg.transpose() * r_k) / sigma_eps;
				// d(neg_ll)/d(log_sigma_eps)
				grad[m_p]     -= post_k * dll_dlse;
				// d(neg_ll)/d(log_sigma_u)
				grad[m_p + 1] -= post_k * (u_vals[k] / sigma_eps) * r_k_sum;
			}
		}

		return total_nll;
	}

	Eigen::MatrixXd hessian(const Eigen::VectorXd& par) {
		const Eigen::VectorXd beta     = par.head(m_p);
		const double log_sigma_eps     = par[m_p];
		const double log_sigma_u       = par[m_p + 1];
		const int n_par = m_p + 2;
		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n_par, n_par);

		if (!std::isfinite(log_sigma_eps) || std::abs(log_sigma_eps) > m_max_abs_log_sigma) {
			H.setConstant(NA_REAL);
			return H;
		}
		if (!std::isfinite(log_sigma_u) || std::abs(log_sigma_u) > m_max_abs_log_sigma) {
			H.setConstant(NA_REAL);
			return H;
		}

		const double sigma_eps = std::exp(log_sigma_eps);
		const double sigma_u = std::exp(log_sigma_u);
		const Eigen::VectorXd eta = m_X * beta;
		const int n_nodes = (int)m_gh.nodes.size();

		Eigen::VectorXd u_vals(n_nodes);
		for (int k = 0; k < n_nodes; ++k)
			u_vals[k] = std::sqrt(2.0) * sigma_u * m_gh.nodes[k];

		for (int g = 0; g < m_n_groups; ++g) {
			const int i_start = m_group_start[g];
			const int sz = m_group_end[g] - i_start;
			const Eigen::MatrixXd Xg = m_X.middleRows(i_start, sz);

			Eigen::VectorXd log_terms(n_nodes);
			std::vector<Eigen::VectorXd> node_score(n_nodes, Eigen::VectorXd::Zero(n_par));
			std::vector<Eigen::MatrixXd> node_hess(n_nodes, Eigen::MatrixXd::Zero(n_par, n_par));

			for (int k = 0; k < n_nodes; ++k) {
				double ll = m_gh.log_norm_weights[k];
				const double u = u_vals[k];
				const double v = u / sigma_eps;
				Eigen::VectorXd& a = node_score[k];
				Eigen::MatrixXd& B = node_hess[k];

				for (int r = 0; r < sz; ++r) {
					const int i = i_start + r;
					double wik = (m_log_y[i] - eta[i] - u) / sigma_eps;
					if (wik > 700.0) wik = 700.0;
					const double ewik = std::exp(wik);
					const double dead_i = m_dead[i];
					const double resid = ewik - dead_i;

					if (dead_i > 0.5) {
						ll += wik - log_sigma_eps - m_log_y[i] - ewik;
					} else {
						ll -= ewik;
					}

					Eigen::VectorXd w_grad = Eigen::VectorXd::Zero(n_par);
					w_grad.head(m_p).noalias() = -Xg.row(r).transpose() / sigma_eps;
					w_grad[m_p] = -wik;
					w_grad[m_p + 1] = -v;

					Eigen::MatrixXd w_hess = Eigen::MatrixXd::Zero(n_par, n_par);
					w_hess.block(0, m_p, m_p, 1).noalias() = Xg.row(r).transpose() / sigma_eps;
					w_hess.block(m_p, 0, 1, m_p) = w_hess.block(0, m_p, m_p, 1).transpose();
					w_hess(m_p, m_p) = wik;
					w_hess(m_p, m_p + 1) = v;
					w_hess(m_p + 1, m_p) = v;
					w_hess(m_p + 1, m_p + 1) = -v;

					a.noalias() -= resid * w_grad;
					a[m_p] -= dead_i;
					B.noalias() -= ewik * (w_grad * w_grad.transpose()) + resid * w_hess;
				}
				log_terms[k] = ll;
			}

			const double ll_g = log_sum_exp_wf(log_terms);
			if (!std::isfinite(ll_g)) {
				H.setConstant(NA_REAL);
				return H;
			}

			Eigen::VectorXd mean_score = Eigen::VectorXd::Zero(n_par);
			Eigen::MatrixXd mean_hess = Eigen::MatrixXd::Zero(n_par, n_par);
			Eigen::MatrixXd mean_outer = Eigen::MatrixXd::Zero(n_par, n_par);
			for (int k = 0; k < n_nodes; ++k) {
				const double post_k = std::exp(log_terms[k] - ll_g);
				if (post_k < 1e-15) continue;
				mean_score.noalias() += post_k * node_score[k];
				mean_hess.noalias() += post_k * node_hess[k];
				mean_outer.noalias() += post_k * (node_score[k] * node_score[k].transpose());
			}

			H.noalias() -= mean_hess + mean_outer - mean_score * mean_score.transpose();
		}

		return (H + H.transpose()) / 2.0;
	}
};

} // namespace

// [[Rcpp::export]]
List fast_weibull_frailty_cpp(
	const Eigen::VectorXd& y,
	const Eigen::VectorXd& dead,
	const Eigen::MatrixXd& X,
	const Eigen::VectorXi& group_id,
	Rcpp::Nullable<Rcpp::NumericVector> start = R_NilValue,
	bool estimate_only = false,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0,
	int maxit = 300,
	double eps_g = 1e-6,
	std::string optimization_alg = "lbfgs"
) {
	const int p     = X.cols();
	const int n_par = p + 2;

	WeibullFrailtyLikelihood obj(y, dead, X, group_id, n_gh, max_abs_log_sigma);

	Eigen::VectorXd par(n_par);
	if (start.isNotNull()) {
		par = Rcpp::as<Eigen::VectorXd>(Rcpp::NumericVector(start));
	} else {
		Eigen::VectorXd log_y = y.array().log().matrix();
		Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X);
		par.head(p) = cod.solve(log_y);
		Eigen::VectorXd resid = log_y - X * par.head(p);
		const double std_resid = std::sqrt(resid.squaredNorm() / std::max(1, (int)y.size() - p));
		par[p]     = std::log(std_resid * 0.7797);
		par[p + 1] = 0.0;
	}

	FixedParamSpec fixed_spec = make_fixed_param_spec(n_par);

	double neg_ll  = NA_REAL;
	bool converged = false;
	try {
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs");
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
	if (!estimate_only && converged) {
		Eigen::MatrixXd info = obj.hessian(par);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(info);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv = ldlt.solve(Eigen::MatrixXd::Identity(info.rows(), info.cols()));
			if (inv.allFinite()) ssq_b_T = inv(0, 0);  // j_T = 0
		}
	}

	return List::create(
		Named("b")             = par.head(p),
		Named("log_sigma_eps") = par[p],
		Named("log_sigma_u")   = par[p + 1],
		Named("ssq_b_T")       = ssq_b_T,
		Named("converged")     = converged,
		Named("neg_loglik")    = neg_ll
	);
}
