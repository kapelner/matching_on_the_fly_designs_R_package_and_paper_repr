#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

double clamp_eta_for_exp(double eta) {
	return std::min(eta, 700.0);
}

ModelResult fast_poisson_internal(const Eigen::MatrixXd& X,
							 const Eigen::VectorXd& y,
							 int maxit = 100,
							 double tol = 1e-8) {
	const int n = X.rows();
	const int p = X.cols();
    ModelResult res;

	res.b = VectorXd::Zero(p);
	VectorXd mu = (y.array() + 0.1).matrix();
	VectorXd eta = mu.array().log().matrix();

	VectorXd XtWz(p);

	for (int iter = 0; iter < maxit; ++iter) {
		for (int i = 0; i < n; ++i) {
			eta[i] = clamp_eta_for_exp(eta[i]);
		}
		mu = eta.array().exp().matrix();
		mu = mu.array().max(1e-10);

		VectorXd z = eta + (y - mu).cwiseQuotient(mu);
		res.XtWX.noalias() = X.transpose() * mu.asDiagonal() * X;
		XtWz.noalias() = X.transpose() * (mu.cwiseProduct(z));

		VectorXd beta_new = res.XtWX.ldlt().solve(XtWz);
		if (!beta_new.allFinite()) {
			break;
		}

		if ((beta_new - res.b).norm() < tol) {
			res.b = beta_new;
			res.converged = true;
			break;
		}

		res.b = beta_new;
		eta.noalias() = X * res.b;
	}

	eta.noalias() = X * res.b;
	for (int i = 0; i < n; ++i) {
		eta[i] = clamp_eta_for_exp(eta[i]);
	}
	res.mu = eta.array().exp().matrix();
	res.mu = res.mu.array().max(1e-10);
	res.XtWX.noalias() = X.transpose() * res.mu.asDiagonal() * X;

	return res;
}

}

// [[Rcpp::export]]
List fast_poisson_regression_cpp(const Eigen::MatrixXd& X,
									 const Eigen::VectorXd& y,
									 int maxit = 100,
									 double tol = 1e-8) {
	ModelResult res = fast_poisson_internal(X, y, maxit, tol);
    return List::create(
		Named("b") = res.b,
		Named("mu") = res.mu,
		Named("XtWX") = res.XtWX,
		Named("converged") = res.converged
	);
}

// [[Rcpp::export]]
List fast_poisson_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
											  const Eigen::VectorXd& y,
											  int j = 2,
											  int maxit = 100,
											  double tol = 1e-8) {
	ModelResult res = fast_poisson_internal(Xmm, y, maxit, tol);
	res.ssq_b_j = compute_diagonal_inverse_entry(res.XtWX, j);
	res.ssq_b_2 = (Xmm.cols() >= 2) ? compute_diagonal_inverse_entry(res.XtWX, 2) : NA_REAL;

	return List::create(
		Named("b") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("mu") = res.mu,
		Named("converged") = res.converged
	);
}

// [[Rcpp::export]]
List fast_quasipoisson_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
												   const Eigen::VectorXd& y,
												   int j = 2,
												   int maxit = 100,
												   double tol = 1e-8) {
	ModelResult res = fast_poisson_internal(Xmm, y, maxit, tol);

	const int df_resid = Xmm.rows() - Xmm.cols();
	if (df_resid > 0) {
		ArrayXd pearson_terms = ((y - res.mu).array().square()) / res.mu.array();
		res.dispersion = pearson_terms.sum() / static_cast<double>(df_resid);
		if (std::isfinite(res.dispersion) && res.dispersion > 0) {
			res.ssq_b_j = res.dispersion * compute_diagonal_inverse_entry(res.XtWX, j);
			if (Xmm.cols() >= 2) {
				res.ssq_b_2 = res.dispersion * compute_diagonal_inverse_entry(res.XtWX, 2);
			}
		}
	}

	return List::create(
		Named("b") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("dispersion") = res.dispersion,
		Named("mu") = res.mu,
		Named("converged") = res.converged
	);
}

//' Parallel Poisson Randomization Distribution
//'
//' @param y Numeric vector of response values (pre-null-shifted for treated).
//' @param X_covars Matrix of covariates (without intercept or treatment).
//' @param w_mat Integer matrix of permuted treatment assignments (n x nsim).
//' @param delta Null treatment effect shift.
//' @param log_transform If TRUE, apply multiplicative delta shift (exp scale); otherwise additive.
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of length nsim with treatment coefficients.
// [[Rcpp::export]]
NumericVector compute_poisson_distr_parallel_cpp(
	const Eigen::VectorXd& y,
	const Eigen::MatrixXd& X_covars,
	const Rcpp::IntegerMatrix& w_mat,
	double delta,
	bool log_transform,
	int num_cores
) {
	int nsim = w_mat.cols();
	int n = y.size();
	int p_covars = X_covars.cols();
	int p_full = p_covars + 2; // intercept + treatment + covars

	std::vector<double> results(nsim, NA_REAL);
	const int* w_ptr = w_mat.begin();

#ifdef _OPENMP
	omp_set_num_threads(num_cores);
#endif

	const double exp_delta = std::exp(delta);

#pragma omp parallel for schedule(static)
	for (int b = 0; b < nsim; ++b) {
		const int* w_col = w_ptr + (size_t)b * n;

		Eigen::MatrixXd X_full(n, p_full);
		Eigen::VectorXd y_shifted(n);

		for (int i = 0; i < n; ++i) {
			X_full(i, 0) = 1.0;
			X_full(i, 1) = (double)w_col[i];
			for (int k = 0; k < p_covars; ++k) {
				X_full(i, 2 + k) = X_covars(i, k);
			}
			bool treated = (w_col[i] == 1);
			if (log_transform) {
				y_shifted[i] = treated ? y[i] * exp_delta : y[i];
			} else {
				y_shifted[i] = treated ? y[i] + delta : y[i];
			}
		}

		ModelResult res = fast_poisson_internal(X_full, y_shifted);
		if (res.converged && p_full >= 2 && std::isfinite(res.b[1])) {
			results[b] = res.b[1];
		}
	}

	return wrap(results);
}
