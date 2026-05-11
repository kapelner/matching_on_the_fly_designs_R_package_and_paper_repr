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
		const std::vector<int>& y_r, // 1-indexed
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
			y_s[i] = y_r[ord[i]];
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

		const Eigen::VectorXd eta_all = dat.X_s * beta;
		Eigen::Map<const Eigen::VectorXi> y_all(dat.y_s.data(), dat.n);

		double total_nll = 0.0;
		grad.setZero();

		Eigen::MatrixXd log_terms_mat(dat.G, n_nodes);
		std::vector<Eigen::VectorXd> dlogp_deta_all_nodes(n_nodes);
		std::vector<std::vector<Eigen::VectorXd>> dlogp_dalpha_all_nodes(n_nodes, std::vector<Eigen::VectorXd>(n_alpha));

		for (int k = 0; k < n_nodes; k++) {
			const double u = b_vals[k];
			const Eigen::VectorXd eta_k_all = eta_all.array() + u;
			
			// Matrix of size n x n_alpha storing F_k = plogis(alpha_j - eta_k)
			Eigen::MatrixXd F_mat(dat.n, n_alpha);
			for (int j = 0; j < n_alpha; j++) {
				F_mat.col(j) = plogis_array_safe(alpha[j] - eta_k_all.array()).matrix();
			}

			Eigen::VectorXd prob_all(dat.n);
			Eigen::VectorXd deta_all(dat.n);
			std::vector<Eigen::VectorXd> dalpha_all(n_alpha, Eigen::VectorXd::Zero(dat.n));

			for (int i = 0; i < dat.n; i++) {
				int yi = y_all[i];
				double F_up = (yi >= dat.K) ? 1.0 : F_mat(i, yi - 1);
				double F_lo = (yi <= 1)      ? 0.0 : F_mat(i, yi - 2);
				double prob = std::max(1e-15, F_up - F_lo);
				prob_all[i] = prob;

				double f_up = (yi >= dat.K) ? 0.0 : F_up * (1.0 - F_up);
				double f_lo = (yi <= 1)      ? 0.0 : F_lo * (1.0 - F_lo);
				deta_all[i] = -(f_up - f_lo) / prob;
				if (yi < dat.K) dalpha_all[yi - 1][i] = f_up / prob;
				if (yi > 1)      dalpha_all[yi - 2][i] = -f_lo / prob;
			}

			dlogp_deta_all_nodes[k] = deta_all;
			for(int j=0; j<n_alpha; j++) dlogp_dalpha_all_nodes[k][j] = dalpha_all[j];

			Eigen::ArrayXd log_prob_all = prob_all.array().log();
			for (int gi = 0; gi < dat.G; gi++) {
				log_terms_mat(gi, k) = dat.gh.log_norm_weights[k] + 
				                       log_prob_all.segment(dat.grp_start[gi], dat.grp_size[gi]).sum();
			}
		}

		Eigen::VectorXd ll_g_vec(dat.G);
		for (int gi = 0; gi < dat.G; gi++) {
			ll_g_vec[gi] = log_sum_exp_v(log_terms_mat.row(gi));
			total_nll -= ll_g_vec[gi];
		}

		for (int k = 0; k < n_nodes; k++) {
			Eigen::VectorXd pk_vec(dat.G);
			for (int gi = 0; gi < dat.G; gi++) pk_vec[gi] = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);

			Eigen::VectorXd pk_expanded(dat.n);
			for (int gi = 0; gi < dat.G; gi++) 
				pk_expanded.segment(dat.grp_start[gi], dat.grp_size[gi]).setConstant(pk_vec[gi]);

			// dLL/dalpha
			for (int j = 0; j < n_alpha; j++) {
				grad[j] -= pk_expanded.dot(dlogp_dalpha_all_nodes[k][j]);
			}

			// dLL/dbeta
			Eigen::VectorXd de_pk = pk_expanded.cwiseProduct(dlogp_deta_all_nodes[k]);
			grad.segment(n_alpha, dat.p).noalias() -= dat.X_s.transpose() * de_pk;

			// dLL/dlog_sigma
			double dll_dsigma_k = de_pk.sum() * std::sqrt(2.0) * dat.gh.nodes[k];
			grad[n_alpha + dat.p] -= dll_dsigma_k * sigma;
		}

		// Apply chain rule for cutpoint parameterization
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
		const int n_alpha = dat.K - 1;
		const int total = n_alpha + dat.p + 1;
		const double log_sigma = par[n_alpha + dat.p];
		const double sigma = std::exp(log_sigma);
		const Eigen::VectorXd beta = par.segment(n_alpha, dat.p);
		const Eigen::VectorXd b_vals = std::sqrt(2.0) * sigma * dat.gh.nodes;
		const int n_nodes = (int)b_vals.size();

		std::vector<double> alpha(n_alpha);
		alpha[0] = par[0];
		for (int k = 1; k < n_alpha; ++k) alpha[k] = alpha[k - 1] + std::exp(par[k]);

		const Eigen::VectorXd eta_all = dat.X_s * beta;
		Eigen::Map<const Eigen::VectorXi> y_all(dat.y_s.data(), dat.n);

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(total, total);

		Eigen::MatrixXd log_terms_mat(dat.G, n_nodes);
		std::vector<Eigen::VectorXd> prob_all_nodes(n_nodes);
		std::vector<Eigen::VectorXd> de_all_nodes(n_nodes);
		std::vector<Eigen::VectorXd> d2e_all_nodes(n_nodes);
		std::vector<std::vector<Eigen::VectorXd>> dL_da_all_nodes(n_nodes, std::vector<Eigen::VectorXd>(n_alpha));
		std::vector<std::vector<Eigen::VectorXd>> d2L_daa_all_nodes(n_nodes, std::vector<Eigen::VectorXd>(n_alpha));
		std::vector<std::vector<Eigen::VectorXd>> d2L_dab_all_nodes(n_nodes, std::vector<Eigen::VectorXd>(n_alpha));
		std::vector<Eigen::VectorXd> d2L_da_prev_all_nodes(n_nodes, Eigen::VectorXd::Zero(dat.n)); // for H(j, j-1)

		for (int k = 0; k < n_nodes; k++) {
			const double u = b_vals[k];
			const Eigen::VectorXd eta_k_all = eta_all.array() + u;
			
			Eigen::MatrixXd F_mat(dat.n, n_alpha);
			for (int j = 0; j < n_alpha; j++) F_mat.col(j) = plogis_array_safe(alpha[j] - eta_k_all.array()).matrix();

			Eigen::VectorXd prob_all(dat.n);
			Eigen::VectorXd de_all(dat.n), d2e_all(dat.n);
			for(int j=0; j<n_alpha; j++) {
				dL_da_all_nodes[k][j].setZero(dat.n);
				d2L_daa_all_nodes[k][j].setZero(dat.n);
				d2L_dab_all_nodes[k][j].setZero(dat.n);
			}

			for (int i = 0; i < dat.n; i++) {
				int yi = y_all[i];
				double Fup = (yi >= dat.K) ? 1.0 : F_mat(i, yi - 1);
				double Flo = (yi <= 1)      ? 0.0 : F_mat(i, yi - 2);
				double prob = std::max(1e-15, Fup - Flo);
				prob_all[i] = prob;

				double fup = (yi >= dat.K) ? 0.0 : Fup * (1.0 - Fup);
				double flo = (yi <= 1)      ? 0.0 : Flo * (1.0 - Flo);
				double hup = (yi >= dat.K) ? 0.0 : fup * (1.0 - 2.0 * Fup);
				double hlo = (yi <= 1)      ? 0.0 : flo * (1.0 - 2.0 * Flo);

				de_all[i] = -(fup - flo) / prob;
				d2e_all[i] = -(hup - hlo) / prob + (fup-flo)*(fup-flo)/(prob*prob);

				if (yi < dat.K) {
					dL_da_all_nodes[k][yi-1][i] = fup/prob;
					d2L_daa_all_nodes[k][yi-1][i] = hup/prob - (fup*fup)/(prob*prob);
					d2L_dab_all_nodes[k][yi-1][i] = -hup/prob + (fup*(fup-flo))/(prob*prob);
				}
				if (yi > 1) {
					dL_da_all_nodes[k][yi-2][i] = -flo/prob;
					d2L_daa_all_nodes[k][yi-2][i] = -hlo/prob - (flo*flo)/(prob*prob);
					d2L_dab_all_nodes[k][yi-2][i] = hlo/prob - (flo*(fup-flo))/(prob*prob);
				}
				if (yi > 1 && yi < dat.K) {
					d2L_da_prev_all_nodes[k][i] = (fup * flo) / (prob * prob);
				}
			}
			prob_all_nodes[k] = prob_all;
			de_all_nodes[k] = de_all;
			d2e_all_nodes[k] = d2e_all;

			Eigen::ArrayXd log_prob_all = prob_all.array().log();
			for (int gi = 0; gi < dat.G; gi++) {
				log_terms_mat(gi, k) = dat.gh.log_norm_weights[k] + 
				                       log_prob_all.segment(dat.grp_start[gi], dat.grp_size[gi]).sum();
			}
		}

		Eigen::VectorXd ll_g_vec(dat.G);
		for (int gi = 0; gi < dat.G; gi++) ll_g_vec[gi] = log_sum_exp_v(log_terms_mat.row(gi));

		Eigen::MatrixXd E_Hik_sum = Eigen::MatrixXd::Zero(total, total);
		Eigen::MatrixXd E_GiGiT_sum = Eigen::MatrixXd::Zero(total, total);
		Eigen::MatrixXd G_avg_outer_sum = Eigen::MatrixXd::Zero(total, total);

		for (int k = 0; k < n_nodes; k++) {
			Eigen::VectorXd pk_vec(dat.G);
			for (int gi = 0; gi < dat.G; gi++) pk_vec[gi] = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
			Eigen::VectorXd pk_exp(dat.n);
			for (int gi = 0; gi < dat.G; gi++) pk_exp.segment(dat.grp_start[gi], dat.grp_size[gi]).setConstant(pk_vec[gi]);

			Eigen::MatrixXd H_ik_k = Eigen::MatrixXd::Zero(total, total);
			for(int j=0; j<n_alpha; j++) {
				H_ik_k(j, j) = pk_exp.dot(d2L_daa_all_nodes[k][j]);
				H_ik_k.block(j, n_alpha, 1, dat.p).noalias() = (pk_exp.cwiseProduct(d2L_dab_all_nodes[k][j])).transpose() * dat.X_s;
			}
			for(int j=1; j<n_alpha; j++) {
				// This is H(j, j-1)
				H_ik_k(j, j-1) = pk_exp.dot(d2L_da_prev_all_nodes[k]); // only non-zero when yi = j+1
				H_ik_k(j-1, j) = H_ik_k(j, j-1);
			}

			H_ik_k.block(n_alpha, n_alpha, dat.p, dat.p).noalias() = weighted_crossprod(dat.X_s, pk_exp.cwiseProduct(d2e_all_nodes[k]));
			
			const double node_factor = std::sqrt(2.0) * dat.gh.nodes[k];
			H_ik_k.block(n_alpha, total-1, dat.p, 1).noalias() = (dat.X_s.transpose() * pk_exp.cwiseProduct(d2e_all_nodes[k])) * (node_factor * sigma);
			for(int j=0; j<n_alpha; j++) H_ik_k(j, total-1) = pk_exp.dot(d2L_dab_all_nodes[k][j]) * (node_factor * sigma);
			
			H_ik_k(total-1, total-1) = pk_exp.dot(d2e_all_nodes[k]) * (node_factor * node_factor * sigma * sigma);
			// plus G_ik[total-1] part (sigma chain rule)
			double g_sigma_k = (pk_exp.cwiseProduct(de_all_nodes[k])).sum() * node_factor * sigma;
			H_ik_k(total-1, total-1) += g_sigma_k;

			E_Hik_sum += H_ik_k;
		}

		// Outer product part (G_avg)
		for (int gi = 0; gi < dat.G; gi++) {
			Eigen::VectorXd G_avg_gi = Eigen::VectorXd::Zero(total);
			Eigen::MatrixXd E_GiGiT_gi = Eigen::MatrixXd::Zero(total, total);
			const int start = dat.grp_start[gi];
			const int sz = dat.grp_size[gi];

			for (int k = 0; k < n_nodes; k++) {
				const double pk = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
				if (pk < 1e-15) continue;
				
				Eigen::VectorXd G_ik = Eigen::VectorXd::Zero(total);
				for(int j=0; j<n_alpha; j++) G_ik[j] = dL_da_all_nodes[k][j].segment(start, sz).sum();
				G_ik.segment(n_alpha, dat.p).noalias() = dat.X_s.middleRows(start, sz).transpose() * de_all_nodes[k].segment(start, sz);
				G_ik[total-1] = de_all_nodes[k].segment(start, sz).sum() * (std::sqrt(2.0) * dat.gh.nodes[k] * sigma);

				G_avg_gi.noalias() += pk * G_ik;
				E_GiGiT_gi.noalias() += pk * (G_ik * G_ik.transpose());
			}
			E_GiGiT_sum += E_GiGiT_gi;
			G_avg_outer_sum += G_avg_gi * G_avg_gi.transpose();
		}

		// Symmetrize E_Hik_sum and apply cutpoint Jacobian J
		for(int r1=0; r1<total; r1++) for(int c1=0; c1<r1; c1++) E_Hik_sum(r1, c1) = E_Hik_sum(c1, r1);
		
		Eigen::MatrixXd J = Eigen::MatrixXd::Zero(total, total);
		J(0, 0) = 1.0;
		for (int r1 = 1; r1 < n_alpha; r1++) { J(r1, 0) = 1.0; for (int c1 = 1; c1 <= r1; c1++) J(r1, c1) = std::exp(par[c1]); }
		for (int i = n_alpha; i < total; i++) J(i, i) = 1.0;

		Eigen::MatrixXd E_Hik_final = J.transpose() * E_Hik_sum * J;
		// Add diagonal term for cutpoint exp-parameterization
		for (int j = 1; j < n_alpha; j++) {
			double sum_dL_da = 0.0;
			for (int k = 0; k < n_nodes; k++) {
				Eigen::VectorXd pk_vec(dat.G);
				for (int gi = 0; gi < dat.G; gi++) pk_vec[gi] = std::exp(log_terms_mat(gi, k) - ll_g_vec[gi]);
				for (int r1 = j; r1 < n_alpha; r1++) {
					for(int gi=0; gi<dat.G; gi++) sum_dL_da += pk_vec[gi] * dL_da_all_nodes[k][r1].segment(dat.grp_start[gi], dat.grp_size[gi]).sum();
				}
			}
			E_Hik_final(j, j) += sum_dL_da * std::exp(par[j]);
		}

		H = -(E_Hik_final + J.transpose() * (E_GiGiT_sum - G_avg_outer_sum) * J);
		return (H + H.transpose()) / 2.0;
	}
};

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_ordinal_glmm_score_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXi& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int K,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	std::vector<int> y_v(y.size()), gid_v(group_id.size());
	for (int i = 0; i < y.size(); ++i) { y_v[i] = y[i]; gid_v[i] = group_id[i]; }

	OrdinalGLMMData dat(X, y_v, gid_v, K, n_gh, max_abs_log_sigma);
	OrdinalGLMMObjective obj(dat);

	Eigen::VectorXd grad(params.size());
	obj(params, grad);
	return -grad;
}

// [[Rcpp::export]]
Eigen::MatrixXd get_ordinal_glmm_hessian_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXi& y,
	const Eigen::VectorXi& group_id,
	const Eigen::VectorXd& params,
	int K,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0
) {
	std::vector<int> y_v(y.size()), gid_v(group_id.size());
	for (int i = 0; i < y.size(); ++i) { y_v[i] = y[i]; gid_v[i] = group_id[i]; }

	OrdinalGLMMData dat(X, y_v, gid_v, K, n_gh, max_abs_log_sigma);
	OrdinalGLMMObjective obj(dat);

	return -obj.hessian(params);
}

// [[Rcpp::export]]
List fast_ordinal_glmm_cpp(
	const Eigen::MatrixXd& X,     // n x p, NO intercept column; treatment at column j_T (0-based)
	const Eigen::VectorXi& y,     // 1-indexed ordinal outcomes, length n
	const Eigen::VectorXi& group_id, // group IDs, length n (sorted internally)
	int K,                        // number of ordinal levels
	int j_T,                      // 0-based treatment column index in X
	bool estimate_only = false,
	int n_gh = 20,
	double max_abs_log_sigma = 8.0,
	int maxit = 300,
	double eps_g = 1e-6,
	Rcpp::Nullable<Rcpp::NumericVector> start = R_NilValue,  // optional warm start
	std::string optimization_alg = "lbfgs",
	Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
	Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
	const int n = X.rows();
	const int p = X.cols();
	const int n_alpha = K - 1;
	const int total = n_alpha + p + 1; // cutpoint params + betas + log_sigma
	FixedParamSpec fixed_spec = make_fixed_param_spec(total, fixed_idx, fixed_values);

	// Convert Eigen/R vectors to std::vector for OrdinalGLMMData
	std::vector<int> y_v(n), gid_v(n);
	for (int i = 0; i < n; ++i) { y_v[i] = y[i]; gid_v[i] = group_id[i]; }

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
		LikelihoodFitResult fit = optimize_fixed_likelihood(obj, par, fixed_spec, maxit, eps_g, optimization_alg, "lbfgs");
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
		FixedParameterFunctor<OrdinalGLMMObjective> fixed_obj(obj, fixed_spec, par);
		Eigen::VectorXd params_free = subset_vector(par, fixed_spec.free_idx);
		Eigen::MatrixXd H = fixed_obj.hessian(params_free);
		Eigen::LDLT<Eigen::MatrixXd> ldlt(H);
		if (ldlt.info() == Eigen::Success) {
			Eigen::MatrixXd inv_free = ldlt.solve(Eigen::MatrixXd::Identity(H.rows(), H.cols()));
			if (inv_free.allFinite()) {
				Eigen::MatrixXd inv = expand_free_covariance(total, fixed_spec, inv_free, true);
				if (j_T_full < total) {
					ssq_b_T = inv(j_T_full, j_T_full);
				}
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
