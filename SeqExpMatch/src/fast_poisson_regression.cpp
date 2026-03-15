#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

namespace {

double clamp_eta_for_exp(double eta) {
	return std::min(eta, 700.0);
}

List fast_poisson_core_cpp(const Eigen::MatrixXd& X,
							 const Eigen::VectorXd& y,
							 int maxit = 100,
							 double tol = 1e-8) {
	const int n = X.rows();
	const int p = X.cols();

	VectorXd beta = VectorXd::Zero(p);
	VectorXd mu = (y.array() + 0.1).matrix();
	VectorXd eta = mu.array().log().matrix();

	MatrixXd XtWX(p, p);
	VectorXd XtWz(p);
	bool converged = false;

	for (int iter = 0; iter < maxit; ++iter) {
		for (int i = 0; i < n; ++i) {
			eta[i] = clamp_eta_for_exp(eta[i]);
		}
		mu = eta.array().exp().matrix();
		mu = mu.array().max(1e-10);

		VectorXd z = eta + (y - mu).cwiseQuotient(mu);
		XtWX.noalias() = X.transpose() * mu.asDiagonal() * X;
		XtWz.noalias() = X.transpose() * (mu.cwiseProduct(z));

		VectorXd beta_new = XtWX.ldlt().solve(XtWz);
		if (!beta_new.allFinite()) {
			break;
		}

		if ((beta_new - beta).norm() < tol) {
			beta = beta_new;
			converged = true;
			break;
		}

		beta = beta_new;
		eta.noalias() = X * beta;
	}

	eta.noalias() = X * beta;
	for (int i = 0; i < n; ++i) {
		eta[i] = clamp_eta_for_exp(eta[i]);
	}
	mu = eta.array().exp().matrix();
	mu = mu.array().max(1e-10);
	XtWX.noalias() = X.transpose() * mu.asDiagonal() * X;

	return List::create(
		Named("b") = beta,
		Named("mu") = mu,
		Named("XtWX") = XtWX,
		Named("converged") = converged
	);
}

}

//' Fast Poisson Regression using Rcpp and IRLS
//'
//' @param X Design matrix.
//' @param y Response vector.
//' @param maxit Maximum IRLS iterations.
//' @param tol Convergence tolerance for coefficient updates.
//' @return A list with coefficients, fitted means, Fisher information, and convergence flag.
// [[Rcpp::export]]
List fast_poisson_regression_cpp(const Eigen::MatrixXd& X,
									 const Eigen::VectorXd& y,
									 int maxit = 100,
									 double tol = 1e-8) {
	return fast_poisson_core_cpp(X, y, maxit, tol);
}

//' Fast Poisson Regression with variance using Rcpp and IRLS
//'
//' @param Xmm Design matrix.
//' @param y Response vector.
//' @param j Index of the coefficient for which to compute the variance.
//' @param maxit Maximum IRLS iterations.
//' @param tol Convergence tolerance for coefficient updates.
//' @return A list with coefficients and selected coefficient variances.
// [[Rcpp::export]]
List fast_poisson_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
											  const Eigen::VectorXd& y,
											  int j = 2,
											  int maxit = 100,
											  double tol = 1e-8) {
	List fit = fast_poisson_core_cpp(Xmm, y, maxit, tol);
	MatrixXd XtWX = fit["XtWX"];
	VectorXd b = fit["b"];

	double ssq_b_j = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, j);
	double ssq_b_2 = (Xmm.cols() >= 2) ? eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2) : NA_REAL;

	return List::create(
		Named("b") = b,
		Named("ssq_b_j") = ssq_b_j,
		Named("ssq_b_2") = ssq_b_2,
		Named("mu") = fit["mu"],
		Named("converged") = fit["converged"]
	);
}

//' Fast Quasi-Poisson Regression with variance using Rcpp and IRLS
//'
//' @param Xmm Design matrix.
//' @param y Response vector.
//' @param j Index of the coefficient for which to compute the variance.
//' @param maxit Maximum IRLS iterations.
//' @param tol Convergence tolerance for coefficient updates.
//' @return A list with coefficients, dispersion estimate, and selected coefficient variances.
// [[Rcpp::export]]
List fast_quasipoisson_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
												   const Eigen::VectorXd& y,
												   int j = 2,
												   int maxit = 100,
												   double tol = 1e-8) {
	List fit = fast_poisson_core_cpp(Xmm, y, maxit, tol);
	MatrixXd XtWX = fit["XtWX"];
	VectorXd b = fit["b"];
	VectorXd mu = fit["mu"];

	const int df_resid = Xmm.rows() - Xmm.cols();
	double dispersion = NA_REAL;
	double ssq_b_j = NA_REAL;
	double ssq_b_2 = NA_REAL;

	if (df_resid > 0) {
		ArrayXd pearson_terms = ((y - mu).array().square()) / mu.array();
		dispersion = pearson_terms.sum() / static_cast<double>(df_resid);
		if (std::isfinite(dispersion) && dispersion > 0) {
			ssq_b_j = dispersion * eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, j);
			if (Xmm.cols() >= 2) {
				ssq_b_2 = dispersion * eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2);
			}
		}
	}

	return List::create(
		Named("b") = b,
		Named("ssq_b_j") = ssq_b_j,
		Named("ssq_b_2") = ssq_b_2,
		Named("dispersion") = dispersion,
		Named("mu") = mu,
		Named("converged") = fit["converged"]
	);
}
