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
        eta = eta.cwiseMin(700.0);
        VectorXd mu = eta.array().exp().matrix();
        ArrayXd wt = ArrayXd::Ones(m_n);
        if (m_use_weights) wt = m_weights.array();
        double neg_ll = (wt * (mu.array() - m_y.array() * eta.array())).sum();
        grad = m_X.transpose() * (wt * (mu - m_y).array()).matrix();
        return neg_ll;
    }

    MatrixXd hessian(const VectorXd& beta) {
        VectorXd eta = m_X * beta;
        eta = eta.cwiseMin(700.0);
        VectorXd w = eta.array().exp().matrix();
        if (m_use_weights) w = w.cwiseProduct(m_weights);
        return weighted_crossprod(m_X, w);
    }
};

ModelResult fast_poisson_internal(const Eigen::MatrixXd& X,
							 const Eigen::VectorXd& y,
                             const Eigen::VectorXd& weights = Eigen::VectorXd(),
							 int maxit = 100,
							 double tol = 1e-8,
                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                             std::string optimization_alg = "lbfgs") {
	const int n = X.rows();
	const int p = X.cols();
    bool use_weights = (weights.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    std::string alg = normalize_optimizer_algorithm(optimization_alg, "lbfgs", true);

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
        eta = eta.cwiseMin(700.0);
        res.mu = eta.array().exp().matrix();
        res.XtWX = fun.hessian(res.b);
        res.converged = fit.converged;
        return res;
    }

    const int p_free = fixed_spec.free_idx.size();
    RowMajorMatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
    VectorXd beta_free = VectorXd::Zero(p_free);
    VectorXd eta_fixed = VectorXd::Zero(n);
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        eta_fixed.noalias() += X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }
    ModelResult res;

	res.b = VectorXd::Zero(p);
    VectorXd mu = (y.array() + 0.1).matrix();
    VectorXd eta = VectorXd::Zero(n);
	VectorXd XtWz(p_free);
    VectorXd w(n);
    VectorXd z(n);
    VectorXd z_adj(n);
    VectorXd beta_new(p_free);
    MatrixXd XtWX_free(p_free, p_free);

	for (int iter = 0; iter < maxit; ++iter) {
        eta.noalias() = X_free * beta_free;
        eta += eta_fixed;
		eta = eta.cwiseMin(700.0);
		mu = eta.array().exp().matrix();
		mu = mu.array().max(1e-10);
        
        if (use_weights) {
            w.array() = mu.array() * weights.array();
        } else {
            w = mu;
		}

		z.noalias() = y - mu;
        z.array() /= mu.array();
        z += eta;
        z_adj = z;
        z_adj -= eta_fixed;
		XtWX_free = weighted_crossprod(X_free, w);
		XtWz.noalias() = weighted_crossprod_rhs(X_free, w, z_adj);

		beta_new = XtWX_free.ldlt().solve(XtWz);
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
	eta = eta.cwiseMin(700.0);
	res.mu = eta.array().exp().matrix();
	res.mu = res.mu.array().max(1e-10);
    
    if (use_weights) {
        w = res.mu.cwiseProduct(weights);
    } else {
        w = res.mu;
    }
    MatrixXd info_free = weighted_crossprod(X_free, w);
	res.XtWX = expand_free_covariance(p, fixed_spec, info_free, false);

	return res;
}

} // namespace

ModelResult fast_poisson_regression_internal(const Eigen::MatrixXd& X,
                                             const Eigen::VectorXd& y,
                                             const Eigen::VectorXd& weights = Eigen::VectorXd(),
                                             int maxit = 100,
                                             double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "lbfgs") {
    return fast_poisson_internal(X, y, weights, maxit, tol, fixed_idx, fixed_values, optimization_alg);
}

//' @title Compute Poisson Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a Poisson regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the score.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd get_poisson_regression_score_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return X.transpose() * (y - mu);
}

//' @title Compute Poisson Regression Hessian
//' @description Calculates the Hessian matrix (second derivatives of the log-likelihood) for a Poisson regression model.
//' @param X A numeric matrix of predictors.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the Hessian.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_poisson_regression_hessian_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return -weighted_crossprod(X, mu);
}

//' @title Compute Weighted Poisson Regression Score
//' @description Calculates the score vector for a weighted Poisson regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param weights A numeric vector of weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the weighted score.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd get_poisson_regression_weighted_score_cpp(const Eigen::MatrixXd& X,
                                                          const Eigen::VectorXd& y,
                                                          const Eigen::VectorXd& weights,
                                                          const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    return X.transpose() * weights.cwiseProduct(y - mu);
}

//' @title Compute Weighted Poisson Regression Hessian
//' @description Calculates the Hessian matrix for a weighted Poisson regression model.
//' @param X A numeric matrix of predictors.
//' @param weights A numeric vector of weights.
//' @param beta A numeric vector of coefficients.
//' @return A numeric matrix representing the weighted Hessian.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd get_poisson_regression_weighted_hessian_cpp(const Eigen::MatrixXd& X,
                                                            const Eigen::VectorXd& weights,
                                                            const Eigen::VectorXd& beta) {
    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu = eta.array().exp().matrix();
    Eigen::VectorXd w = weights.cwiseProduct(mu);
    return -weighted_crossprod(X, w);
}

//' @title Fast Poisson Regression (C++)
//' @description High-performance Poisson regression fitting using IRLS or L-BFGS.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm ("lbfgs" or "irls").
//' @return A list containing coefficients, fitted values, and information matrix.
//' @export
// [[Rcpp::export]]
List fast_poisson_regression_cpp(const Eigen::MatrixXd& X,
									 const Eigen::VectorXd& y,
									 int maxit = 100,
									 double tol = 1e-8,
                                     Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "lbfgs") {
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), maxit, tol, fixed_idx, fixed_values, optimization_alg);
    return List::create(
		Named("b") = res.b,
		Named("mu") = res.mu,
		Named("XtWX") = res.XtWX,
		Named("converged") = res.converged
	);
}

//' @title Fast Weighted Poisson Regression (C++)
//' @description High-performance weighted Poisson regression fitting.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param weights A numeric vector of weights.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, fitted values, and information matrix.
//' @export
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

//' @title Fast Poisson Regression with Variance (C++)
//' @description Poisson regression with variance-covariance matrix and score calculation.
//' @param Xmm A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, vcov, score, and likelihood statistics.
//' @export
// [[Rcpp::export]]
List fast_poisson_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
											  const Eigen::VectorXd& y,
											  int j = 2,
											  int maxit = 100,
											  double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "lbfgs") {
	ModelResult res = fast_poisson_internal(Xmm, y, Eigen::VectorXd(), maxit, tol, fixed_idx, fixed_values, optimization_alg);
    FixedParamSpec fixed_spec = make_fixed_param_spec(Xmm.cols(), fixed_idx, fixed_values);
    MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);
    MatrixXd cov_free = info_free.inverse();
    MatrixXd vcov = expand_free_covariance(Xmm.cols(), fixed_spec, cov_free, true);
	res.ssq_b_j = (j > 0 && j <= Xmm.cols()) ? vcov(j - 1, j - 1) : NA_REAL;
	res.ssq_b_2 = (Xmm.cols() >= 2) ? vcov(1, 1) : NA_REAL;
	VectorXd score = get_poisson_regression_score_cpp(Xmm, y, res.b);
	MatrixXd information = res.XtWX;
	Eigen::ArrayXd eta = (Xmm * res.b).array().max(-30.0).min(30.0);
	double neg_loglik = (eta.exp() - y.array() * eta).sum();

	return List::create(
		Named("b") = res.b,
		Named("params") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("mu") = res.mu,
		Named("converged") = res.converged,
        Named("vcov") = vcov,
		Named("score") = score,
		Named("observed_information") = information,
		Named("fisher_information") = information,
		Named("information") = information,
		Named("information_type") = "fisher",
		Named("hessian") = -information,
		Named("neg_loglik") = neg_loglik,
		Named("neg_ll") = neg_loglik,
		Named("loglik") = R_finite(neg_loglik) ? -neg_loglik : NA_REAL
	);
}

//' @title Fast Quasi-Poisson Regression with Variance (C++)
//' @description Quasi-Poisson regression with dispersion-adjusted variance-covariance matrix.
//' @param Xmm A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @return A list containing coefficients, vcov, and dispersion estimate.
//' @export
// [[Rcpp::export]]
List fast_quasipoisson_regression_with_var_cpp(const Eigen::MatrixXd& Xmm,
												   const Eigen::VectorXd& y,
												   int j = 2,
												   int maxit = 100,
												   double tol = 1e-8,
                                                   Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                   Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                   std::string optimization_alg = "lbfgs") {
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

	// Pre-allocate one workspace per thread; fill invariant columns (intercept
	// and covariates) once so the inner loop only overwrites the treatment column.
	struct PoissonWorkspace {
		Eigen::MatrixXd X_full;
		Eigen::VectorXd y_shifted;
	};
	std::vector<PoissonWorkspace> ws(num_cores);
	for (int t = 0; t < num_cores; ++t) {
		ws[t].X_full.resize(n, p_full);
		ws[t].y_shifted.resize(n);
		ws[t].X_full.col(0).setOnes();
		for (int k = 0; k < p_covars; ++k) {
			ws[t].X_full.col(2 + k) = X_covars.col(k);
		}
	}

#pragma omp parallel for schedule(static)
	for (int b = 0; b < nsim; ++b) {
		int tid = 0;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		const int* w_col = w_ptr + (size_t)b * n;
		PoissonWorkspace& W = ws[tid];

		for (int i = 0; i < n; ++i) {
			W.X_full(i, 1) = (double)w_col[i];
			const bool treated = (w_col[i] == 1);
			W.y_shifted[i] = treated ? (log_transform ? y[i] * exp_delta : y[i] + delta) : y[i];
		}

		ModelResult res = fast_poisson_internal(W.X_full, W.y_shifted, Eigen::VectorXd());
		if (res.converged && p_full >= 2 && std::isfinite(res.b[1])) {
			results[b] = res.b[1];
		}
	}

	return wrap(results);
}
