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

class PoissonNegLogLik {
private:
    const MatrixXd& m_X;
    const VectorXd& m_y;
    const VectorXd& m_weights;
    const bool m_use_weights;
    const int m_n;

public:
    PoissonNegLogLik(const MatrixXd& X,
                     const VectorXd& y,
                     const VectorXd& weights) :
        m_X(X), m_y(y), m_weights(weights), m_use_weights(weights.size() == X.rows()), m_n(X.rows()) {}

    double operator()(const VectorXd& beta, VectorXd& grad) {
        VectorXd eta = m_X * beta;
        for (int i = 0; i < m_n; ++i) eta[i] = clamp_eta_for_exp(eta[i]);
        VectorXd mu = eta.array().exp().matrix();
        VectorXd wt = m_use_weights ? m_weights : VectorXd::Ones(m_n);
        double neg_ll = (wt.array() * (mu.array() - m_y.array() * eta.array())).sum();
        grad = m_X.transpose() * wt.cwiseProduct(mu - m_y);
        return neg_ll;
    }

    MatrixXd hessian(const VectorXd& beta) {
        VectorXd eta = m_X * beta;
        for (int i = 0; i < m_n; ++i) eta[i] = clamp_eta_for_exp(eta[i]);
        VectorXd w = eta.array().exp().matrix();
        if (m_use_weights) w = w.cwiseProduct(m_weights);
        return m_X.transpose() * w.asDiagonal() * m_X;
    }
};

ModelResult fast_poisson_internal(const Eigen::MatrixXd& X,
							 const Eigen::VectorXd& y,
                             const Eigen::VectorXd& weights = Eigen::VectorXd(),
							 int maxit = 100,
							 double tol = 1e-8,
                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                             std::string optimization_alg = "irls") {
	const int n = X.rows();
	const int p = X.cols();
    bool use_weights = (weights.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    std::string alg = normalize_optimizer_algorithm(optimization_alg, "irls", true);

    if (alg != "irls") {
        VectorXd beta = VectorXd::Zero(p);
        beta = apply_fixed_values(beta, fixed_spec);
        PoissonNegLogLik fun(X, y, weights);
        LikelihoodFitResult fit = (alg == "lbfgs") ?
            optimize_fixed_likelihood_lbfgs(fun, beta, fixed_spec, maxit, tol) :
            optimize_fixed_likelihood_newton(fun, beta, fixed_spec, maxit, tol);

        ModelResult res;
        res.b = fit.params;
        VectorXd eta = X * res.b;
        for (int i = 0; i < n; ++i) eta[i] = clamp_eta_for_exp(eta[i]);
        res.mu = eta.array().exp().matrix();
        res.XtWX = fun.hessian(res.b);
        res.converged = fit.converged;
        return res;
    }

    const int p_free = fixed_spec.free_idx.size();
    MatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
    VectorXd beta_free = VectorXd::Zero(p_free);
    VectorXd eta_fixed = VectorXd::Zero(n);
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }
    ModelResult res;

	res.b = VectorXd::Zero(p);
	VectorXd mu = (y.array() + 0.1).matrix();
	VectorXd eta = eta_fixed + X_free * beta_free;

	VectorXd XtWz(p_free);
    VectorXd w(n);

	for (int iter = 0; iter < maxit; ++iter) {
        eta = eta_fixed + X_free * beta_free;
		for (int i = 0; i < n; ++i) {
			eta[i] = clamp_eta_for_exp(eta[i]);
		}
		mu = eta.array().exp().matrix();
		mu = mu.array().max(1e-10);
        
        if (use_weights) {
            w = mu.cwiseProduct(weights);
        } else {
            w = mu;
		}

		VectorXd z = eta + (y - mu).cwiseQuotient(mu);
        VectorXd z_adj = z - eta_fixed;
        MatrixXd XtWX_free = X_free.transpose() * w.asDiagonal() * X_free;
		XtWz.noalias() = X_free.transpose() * (w.cwiseProduct(z_adj));

		VectorXd beta_new = XtWX_free.ldlt().solve(XtWz);
		if (!beta_new.allFinite()) {
			break;
		}

		if ((beta_new - beta_free).norm() < tol) {
			beta_free = beta_new;
			res.converged = true;
			break;
		}

		beta_free = beta_new;
	}

    for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = beta_free[j];
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];

	eta.noalias() = X * res.b;
	for (int i = 0; i < n; ++i) {
		eta[i] = clamp_eta_for_exp(eta[i]);
	}
	res.mu = eta.array().exp().matrix();
	res.mu = res.mu.array().max(1e-10);
    
    if (use_weights) {
        w = res.mu.cwiseProduct(weights);
    } else {
        w = res.mu;
    }
    MatrixXd info_free = X_free.transpose() * w.asDiagonal() * X_free;
	res.XtWX = expand_free_covariance(p, fixed_spec, info_free, false);

	return res;
}

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_poisson_regression_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return X.transpose() * (y - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_poisson_regression_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return - (X.transpose() * mu.asDiagonal() * X);
}

// [[Rcpp::export]]
Eigen::VectorXd get_poisson_regression_weighted_score_cpp(const Eigen::MatrixXd& X,
                                                          const Eigen::VectorXd& y,
                                                          const Eigen::VectorXd& weights,
                                                          const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return X.transpose() * weights.cwiseProduct(y - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_poisson_regression_weighted_hessian_cpp(const Eigen::MatrixXd& X,
                                                            const Eigen::VectorXd& weights,
                                                            const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    Eigen::VectorXd w = weights.cwiseProduct(mu);
    return - (X.transpose() * w.asDiagonal() * X);
}

// [[Rcpp::export]]
List fast_poisson_regression_cpp(const Eigen::MatrixXd& X,
									 const Eigen::VectorXd& y,
									 int maxit = 100,
									 double tol = 1e-8,
                                     Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "irls") {
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), maxit, tol, fixed_idx, fixed_values, optimization_alg);
    return List::create(
		Named("b") = res.b,
		Named("mu") = res.mu,
		Named("XtWX") = res.XtWX,
		Named("converged") = res.converged
	);
}

// [[Rcpp::export]]
List fast_poisson_regression_weighted_cpp(const Eigen::MatrixXd& X,
                                          const Eigen::VectorXd& y,
                                          const Eigen::VectorXd& weights,
                                          int maxit = 100,
                                          double tol = 1e-8,
                                          Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                          Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                          std::string optimization_alg = "irls") {
	ModelResult res = fast_poisson_internal(X, y, weights, maxit, tol, fixed_idx, fixed_values, optimization_alg);
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
											  double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls") {
	ModelResult res = fast_poisson_internal(Xmm, y, Eigen::VectorXd(), maxit, tol, fixed_idx, fixed_values, optimization_alg);
    FixedParamSpec fixed_spec = make_fixed_param_spec(Xmm.cols(), fixed_idx, fixed_values);
    MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    MatrixXd cov_free = info_free.inverse();
    MatrixXd vcov = expand_free_covariance(Xmm.cols(), fixed_spec, cov_free, true);
	res.ssq_b_j = (j > 0 && j <= Xmm.cols()) ? vcov(j - 1, j - 1) : NA_REAL;
	res.ssq_b_2 = (Xmm.cols() >= 2) ? vcov(1, 1) : NA_REAL;

	return List::create(
		Named("b") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("mu") = res.mu,
		Named("converged") = res.converged,
        Named("vcov") = vcov
	);
}

// [[Rcpp::export]]
List fast_quasipoisson_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
												   const Eigen::VectorXd& y,
												   int j = 2,
												   int maxit = 100,
												   double tol = 1e-8,
                                                   Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                   Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                   std::string optimization_alg = "irls") {
	ModelResult res = fast_poisson_internal(Xmm, y, Eigen::VectorXd(), maxit, tol, fixed_idx, fixed_values, optimization_alg);
    FixedParamSpec fixed_spec = make_fixed_param_spec(Xmm.cols(), fixed_idx, fixed_values);
    MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    MatrixXd cov_free = info_free.inverse();
    MatrixXd vcov = expand_free_covariance(Xmm.cols(), fixed_spec, cov_free, true);

	const int df_resid = Xmm.rows() - Xmm.cols();
	if (df_resid > 0) {
		ArrayXd pearson_terms = ((y - res.mu).array().square()) / res.mu.array();
		res.dispersion = pearson_terms.sum() / static_cast<double>(df_resid);
		if (std::isfinite(res.dispersion) && res.dispersion > 0) {
            vcov *= res.dispersion;
			res.ssq_b_j = (j > 0 && j <= Xmm.cols()) ? vcov(j - 1, j - 1) : NA_REAL;
			if (Xmm.cols() >= 2) res.ssq_b_2 = vcov(1, 1);
		}
	}

	return List::create(
		Named("b") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("dispersion") = res.dispersion,
		Named("mu") = res.mu,
		Named("converged") = res.converged,
        Named("vcov") = vcov
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

		ModelResult res = fast_poisson_internal(X_full, y_shifted, Eigen::VectorXd());
		if (res.converged && p_full >= 2 && std::isfinite(res.b[1])) {
			results[b] = res.b[1];
		}
	}

	return wrap(results);
}
