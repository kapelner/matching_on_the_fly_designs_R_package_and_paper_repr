#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_ols_distr_parallel_cpp(
	const Eigen::VectorXd& y,
	const Eigen::MatrixXd& X_covars,
	const Eigen::MatrixXi& w_mat,
	double delta,
	int num_cores) {

	int nsim = w_mat.cols();
	int n = y.size();
	int p_covars = X_covars.cols();
	int p_full = p_covars + 2; // Intercept + w + covars
	NumericVector results(nsim);

#ifdef _OPENMP
	if (num_cores > 1) {
		omp_set_num_threads(num_cores);
	}
#endif

	// MEMOIZATION Strategy:
	// Let X = [1 | w | X_c]. X'X is:
	// [ n      | sum(w) | 1'X_c   ]
	// [ sum(w) | sum(w) | w'X_c   ]
	// [ X_c'1  | X_c'w  | X_c'X_c ]
	//
	// Most blocks are constant.
	double sum_1 = (double)n;
	Eigen::VectorXd Xt_1 = X_covars.colwise().sum(); // X_c'1
	Eigen::MatrixXd XtX_c = X_covars.transpose() * X_covars; // X_c'X_c
	
#pragma omp parallel for schedule(dynamic)
	for (int b = 0; b < nsim; ++b) {
		Eigen::VectorXi w = w_mat.col(b);
		Eigen::VectorXd w_d = w.cast<double>();
		
		// Shift y if delta != 0
		Eigen::VectorXd y_sim = y;
		if (delta != 0) {
			for (int i = 0; i < n; ++i) {
				if (w[i] == 1) y_sim[i] += delta;
			}
		}

		double sum_w = w_d.sum();
		Eigen::VectorXd Xt_w = X_covars.transpose() * w_d; // X_c'w

		Eigen::MatrixXd XtX(p_full, p_full);
		XtX(0, 0) = sum_1;
		XtX(0, 1) = sum_w;
		XtX.row(0).tail(p_covars) = Xt_1.transpose();
		
		XtX(1, 0) = sum_w;
		XtX(1, 1) = sum_w; // w'w = sum(w) since w is binary
		XtX.row(1).tail(p_covars) = Xt_w.transpose();
		
		XtX.col(0).tail(p_covars) = Xt_1;
		XtX.col(1).tail(p_covars) = Xt_w;
		XtX.bottomRightCorner(p_covars, p_covars) = XtX_c;

		Eigen::VectorXd Xty(p_full);
		Xty[0] = y_sim.sum();
		Xty[1] = w_d.dot(y_sim);
		Xty.tail(p_covars) = X_covars.transpose() * y_sim;

		Eigen::VectorXd beta = XtX.ldlt().solve(Xty);
		results[b] = beta[1];
	}

	return results;
}

// [[Rcpp::export]]
NumericVector compute_ols_bootstrap_parallel_cpp(
	const Eigen::VectorXd& y,
	const Eigen::MatrixXd& X_covars,
	const Eigen::VectorXi& w,
	const Eigen::MatrixXi& indices_mat, // B columns, each is n rows of indices (0-based)
	int num_cores) {

	int B = indices_mat.cols();
	int n = y.size();
	int p_covars = X_covars.cols();
	int p_full = p_covars + 2; // Intercept + w + covars
	NumericVector results(B);

#ifdef _OPENMP
	if (num_cores > 1) {
		omp_set_num_threads(num_cores);
	}
#endif

#pragma omp parallel for schedule(dynamic)
	for (int b = 0; b < B; ++b) {
		Eigen::VectorXi idx = indices_mat.col(b);

		// Check if there are any NAs (which we encoded as -1 if failed max attempts)
		if (idx[0] == -1) {
			results[b] = NA_REAL;
			continue;
		}

		// Construct the resampled design matrix and y vector
		Eigen::VectorXd y_b(n);
		Eigen::MatrixXd X_full_b(n, p_full);

		for (int i = 0; i < n; ++i) {
			int row_idx = idx[i];
			y_b[i] = y[row_idx];
			X_full_b(i, 0) = 1.0; // Intercept
			X_full_b(i, 1) = static_cast<double>(w[row_idx]);
			if (p_covars > 0) {
				X_full_b.row(i).tail(p_covars) = X_covars.row(row_idx);
			}
		}

		// OLS via LDLT
		Eigen::MatrixXd XtX = X_full_b.transpose() * X_full_b;
		Eigen::VectorXd Xty = X_full_b.transpose() * y_b;

		Eigen::VectorXd beta = XtX.ldlt().solve(Xty);
		results[b] = beta[1]; // Treatment effect is the 2nd coefficient
	}

	return results;
}
