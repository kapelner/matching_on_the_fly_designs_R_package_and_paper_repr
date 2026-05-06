#include "_helper_functions.h"
#include <RcppEigen.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_kk_compound_distr_parallel_cpp(
	const Eigen::VectorXd& y,
	const Eigen::MatrixXi& w_mat,
	const Eigen::MatrixXi& m_mat,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = y.size();
	std::vector<double> results_vec(nsim, NA_REAL);
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
	{
		std::vector<double> diffs;
		std::vector<int> treated_idx;
		std::vector<int> control_idx;

		#pragma omp for schedule(static)
		for (int b = 0; b < nsim; ++b) {
			const int* w_col = w_mat.data() + static_cast<size_t>(b) * n;
			const int* m_col = m_mat.data() + static_cast<size_t>(b) * n;

			int m = 0;
			for (int i = 0; i < n; ++i) {
				if (m_col[i] > m) m = m_col[i];
			}

			double d_bar = NA_REAL;
			double ssqD_bar = NA_REAL;

			if (m > 0) {
				treated_idx.assign(m, -1);
				control_idx.assign(m, -1);
				for (int i = 0; i < n; ++i) {
					const int match_id = m_col[i];
					if (match_id <= 0) continue;
					const int slot = match_id - 1;
					if (w_col[i] == 1) treated_idx[slot] = i;
					else control_idx[slot] = i;
				}
				diffs.assign(m, 0.0);
				for (int match_id = 1; match_id <= m; ++match_id) {
					const int slot = match_id - 1;
					diffs[slot] = y[treated_idx[slot]] - y[control_idx[slot]];
				}
				double sum_d = 0.0;
				for (double diff : diffs) sum_d += diff;
				d_bar = sum_d / m;
				if (m > 1) {
					double ss = 0.0;
					for (double diff : diffs) {
						const double centered = diff - d_bar;
						ss += centered * centered;
					}
					ssqD_bar = ss / (m - 1) / m;
				}
			}

			int nRT = 0, nRC = 0;
			double sum_T = 0, sum_C = 0;
			for (int i = 0; i < n; ++i) {
				if (m_col[i] == 0) {
					if (w_col[i] == 1) { nRT++; sum_T += y[i]; }
					else { nRC++; sum_C += y[i]; }
				}
			}

			double r_bar = NA_REAL;
			double ssqR = NA_REAL;
			if (nRT > 0 && nRC > 0) {
				double mean_T = sum_T / nRT;
				double mean_C = sum_C / nRC;
				r_bar = mean_T - mean_C;
				if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
					double sq_diff_T = 0, sq_diff_C = 0;
					for (int i = 0; i < n; ++i) {
						if (m_col[i] == 0) {
							if (w_col[i] == 1) sq_diff_T += (y[i] - mean_T) * (y[i] - mean_T);
							else sq_diff_C += (y[i] - mean_C) * (y[i] - mean_C);
						}
					}
					double var_T = sq_diff_T / (nRT - 1);
					double var_C = sq_diff_C / (nRC - 1);
					int nR = nRT + nRC;
					ssqR = (var_T * (nRT - 1) + var_C * (nRC - 1)) / (nR - 2) * (1.0 / nRT + 1.0 / nRC);
				}
			}

			double beta_hat_T = NA_REAL;
			if (nRT <= 1 || nRC <= 1) {
				beta_hat_T = d_bar;
			} else if (m == 0) {
				beta_hat_T = r_bar;
			} else {
				if (!std::isfinite(ssqD_bar) || ssqD_bar <= 0) beta_hat_T = r_bar;
				else if (!std::isfinite(ssqR) || ssqR <= 0) beta_hat_T = d_bar;
				else {
					double w_star = ssqR / (ssqR + ssqD_bar);
					beta_hat_T = w_star * d_bar + (1.0 - w_star) * r_bar;
				}
			}
			res_ptr[b] = beta_hat_T;
		}
	}
	return wrap(results_vec);
}

// [[Rcpp::export]]
NumericVector compute_kk_compound_bootstrap_parallel_cpp(
	const Eigen::MatrixXd& y_mat,
	const Eigen::MatrixXi& w_mat,
	const Eigen::MatrixXi& m_mat,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = w_mat.rows();
	std::vector<double> results_vec(nsim, NA_REAL);
	double* res_ptr = results_vec.data();
	const bool use_parallel = should_parallelize_replicates(nsim, n, num_cores);

#ifdef _OPENMP
	if (use_parallel) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel if(use_parallel)
	{
		std::vector<double> diffs;
		std::vector<int> treated_idx;
		std::vector<int> control_idx;

		#pragma omp for schedule(static)
		for (int b = 0; b < nsim; ++b) {
			const double* y_col = y_mat.data() + static_cast<size_t>(b) * n;
			const int* w_col = w_mat.data() + static_cast<size_t>(b) * n;
			const int* m_col = m_mat.data() + static_cast<size_t>(b) * n;

			int m = 0;
			for (int i = 0; i < n; ++i) {
				if (m_col[i] > m) m = m_col[i];
			}

			double d_bar = NA_REAL;
			double ssqD_bar = NA_REAL;

			if (m > 0) {
				treated_idx.assign(m, -1);
				control_idx.assign(m, -1);
				for (int i = 0; i < n; ++i) {
					const int match_id = m_col[i];
					if (match_id <= 0) continue;
					const int slot = match_id - 1;
					if (w_col[i] == 1) treated_idx[slot] = i;
					else control_idx[slot] = i;
				}
				diffs.assign(m, 0.0);
				for (int match_id = 1; match_id <= m; ++match_id) {
					const int slot = match_id - 1;
					diffs[slot] = y_col[treated_idx[slot]] - y_col[control_idx[slot]];
				}
				double sum_d = 0.0;
				for (double diff : diffs) sum_d += diff;
				d_bar = sum_d / m;
				if (m > 1) {
					double ss = 0.0;
					for (double diff : diffs) {
						const double centered = diff - d_bar;
						ss += centered * centered;
					}
					ssqD_bar = ss / (m - 1) / m;
				}
			}

			int nRT = 0, nRC = 0;
			double sum_T = 0, sum_C = 0;
			for (int i = 0; i < n; ++i) {
				if (m_col[i] == 0) {
					if (w_col[i] == 1) { nRT++; sum_T += y_col[i]; }
					else { nRC++; sum_C += y_col[i]; }
				}
			}

			double r_bar = NA_REAL;
			double ssqR = NA_REAL;
			if (nRT > 0 && nRC > 0) {
				double mean_T = sum_T / nRT;
				double mean_C = sum_C / nRC;
				r_bar = mean_T - mean_C;
				if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
					double sq_diff_T = 0, sq_diff_C = 0;
					for (int i = 0; i < n; ++i) {
						if (m_col[i] == 0) {
							if (w_col[i] == 1) sq_diff_T += (y_col[i] - mean_T) * (y_col[i] - mean_T);
							else sq_diff_C += (y_col[i] - mean_C) * (y_col[i] - mean_C);
						}
					}
					double var_T = sq_diff_T / (nRT - 1);
					double var_C = sq_diff_C / (nRC - 1);
					int nR = nRT + nRC;
					ssqR = (var_T * (nRT - 1) + var_C * (nRC - 1)) / (nR - 2) * (1.0 / nRT + 1.0 / nRC);
				}
			}

			double beta_hat_T = NA_REAL;
			if (nRT <= 1 || nRC <= 1) {
				beta_hat_T = d_bar;
			} else if (m == 0) {
				beta_hat_T = r_bar;
			} else {
				if (!std::isfinite(ssqD_bar) || ssqD_bar <= 0) beta_hat_T = r_bar;
				else if (!std::isfinite(ssqR) || ssqR <= 0) beta_hat_T = d_bar;
				else {
					double w_star = ssqR / (ssqR + ssqD_bar);
					beta_hat_T = w_star * d_bar + (1.0 - w_star) * r_bar;
				}
			}
			res_ptr[b] = beta_hat_T;
		}
	}
	return wrap(results_vec);
}
