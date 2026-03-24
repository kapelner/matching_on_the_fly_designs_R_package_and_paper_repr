#include <RcppEigen.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_ols_distr_parallel_cpp(
	const NumericVector& y,
	const Eigen::MatrixXd& X_covars,
	const IntegerMatrix& w_mat,
	double delta,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = y.size();
	int p_covars = X_covars.cols();
	int p_full = p_covars + 2; // Intercept + w + covars
	
	std::vector<double> results_vec(nsim);
	
	const double* y_ptr = y.begin();
	const int* w_ptr = w_mat.begin();
	double* res_ptr = results_vec.data();

#ifdef _OPENMP
	omp_set_num_threads(num_cores);
#endif

	// MEMOIZATION
	double sum_1 = (double)n;
	Eigen::VectorXd Xt_1 = X_covars.colwise().sum();
	Eigen::MatrixXd XtX_c = X_covars.transpose() * X_covars;
	
#pragma omp parallel for schedule(static)
	for (int b = 0; b < nsim; ++b) {
		const int* w_col = w_ptr + (size_t)b * n;
		
		Eigen::VectorXd w_d(n);
		Eigen::VectorXd y_sim(n);
		double sum_w = 0;
		double sum_y = 0;
		
		for (int i = 0; i < n; ++i) {
			double w_val = (double)w_col[i];
			w_d[i] = w_val;
			sum_w += w_val;
			double y_val = y_ptr[i] + (w_col[i] == 1 ? delta : 0.0);
			y_sim[i] = y_val;
			sum_y += y_val;
		}

		Eigen::VectorXd Xt_w = X_covars.transpose() * w_d;

		Eigen::MatrixXd XtX(p_full, p_full);
		XtX(0, 0) = sum_1;
		XtX(0, 1) = sum_w;
		XtX.row(0).tail(p_covars) = Xt_1.transpose();
		
		XtX(1, 0) = sum_w;
		XtX(1, 1) = sum_w;
		XtX.row(1).tail(p_covars) = Xt_w.transpose();
		
		XtX.col(0).tail(p_covars) = Xt_1;
		XtX.col(1).tail(p_covars) = Xt_w;
		XtX.bottomRightCorner(p_covars, p_covars) = XtX_c;

		Eigen::VectorXd Xty(p_full);
		Xty[0] = sum_y;
		Xty[1] = w_d.dot(y_sim);
		Xty.tail(p_covars) = X_covars.transpose() * y_sim;

		Eigen::VectorXd beta = XtX.ldlt().solve(Xty);
		res_ptr[b] = beta[1];
	}

	return wrap(results_vec);
}

// Bootstrap OLS: for each column of indices_mat (0-based row indices, -1 = NA bootstrap),
// resample y/w/X_covars and return the OLS treatment coefficient.
// [[Rcpp::export]]
NumericVector compute_ols_bootstrap_parallel_cpp(
	const NumericVector& y,
	const Eigen::MatrixXd& X_covars,
	const IntegerVector& w,
	const IntegerMatrix& indices_mat,
	int num_cores) {

	int B = indices_mat.cols();
	int n_boot = indices_mat.rows(); // bootstrap sample size (= n for simple bootstrap)
	int p_covars = X_covars.cols();
	int p_full = p_covars + 2; // intercept + w + covars

	std::vector<double> results_vec(B, NA_REAL);

	const double* y_ptr = y.begin();
	const int* w_ptr = w.begin();
	const int* idx_ptr = indices_mat.begin();

#ifdef _OPENMP
	omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(static)
	for (int b = 0; b < B; ++b) {
		const int* idx_col = idx_ptr + (size_t)b * n_boot;
		if (idx_col[0] < 0) { results_vec[b] = NA_REAL; continue; }

		Eigen::VectorXd y_b(n_boot);
		Eigen::VectorXd w_b(n_boot);
		Eigen::MatrixXd X_b(n_boot, p_covars);

		double sum_1 = (double)n_boot;
		double sum_w = 0.0, sum_y = 0.0;

		for (int i = 0; i < n_boot; ++i) {
			int idx = idx_col[i];
			y_b[i] = y_ptr[idx];
			double wv = (double)w_ptr[idx];
			w_b[i] = wv;
			sum_w += wv;
			sum_y += y_b[i];
			X_b.row(i) = X_covars.row(idx);
		}

		Eigen::VectorXd Xt_1 = X_b.colwise().sum();
		Eigen::VectorXd Xt_w = X_b.transpose() * w_b;
		Eigen::MatrixXd XtX_c = X_b.transpose() * X_b;

		Eigen::MatrixXd XtX(p_full, p_full);
		XtX(0, 0) = sum_1;
		XtX(0, 1) = sum_w;
		XtX.row(0).tail(p_covars) = Xt_1.transpose();

		XtX(1, 0) = sum_w;
		XtX(1, 1) = sum_w;
		XtX.row(1).tail(p_covars) = Xt_w.transpose();

		XtX.col(0).tail(p_covars) = Xt_1;
		XtX.col(1).tail(p_covars) = Xt_w;
		XtX.bottomRightCorner(p_covars, p_covars) = XtX_c;

		Eigen::VectorXd Xty(p_full);
		Xty[0] = sum_y;
		Xty[1] = w_b.dot(y_b);
		Xty.tail(p_covars) = X_b.transpose() * y_b;

		Eigen::VectorXd beta = XtX.ldlt().solve(Xty);
		results_vec[b] = beta[1];
	}

	return wrap(results_vec);
}
