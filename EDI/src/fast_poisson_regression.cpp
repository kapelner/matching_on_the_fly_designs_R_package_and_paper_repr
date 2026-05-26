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
        double neg_ll = 0.0;
        VectorXd diff(m_n);
        
        if (m_use_weights) {
            for (int i = 0; i < m_n; ++i) {
                double ei = eta[i];
                if (ei > 700.0) ei = 700.0;
                double mui = std::exp(ei);
                double wi = m_weights[i];
                neg_ll += wi * (mui - m_y[i] * ei);
                diff[i] = wi * (mui - m_y[i]);
            }
        } else {
            for (int i = 0; i < m_n; ++i) {
                double ei = eta[i];
                if (ei > 700.0) ei = 700.0;
                double mui = std::exp(ei);
                neg_ll += (mui - m_y[i] * ei);
                diff[i] = (mui - m_y[i]);
            }
        }
        grad.noalias() = m_X.transpose() * diff;
        return neg_ll;
    }

    MatrixXd hessian(const VectorXd& beta) {
        VectorXd eta = m_X * beta;
        VectorXd w(m_n);
        if (m_use_weights) {
            for (int i = 0; i < m_n; ++i) {
                double ei = eta[i];
                if (ei > 700.0) ei = 700.0;
                w[i] = std::exp(ei) * m_weights[i];
            }
        } else {
            for (int i = 0; i < m_n; ++i) {
                double ei = eta[i];
                if (ei > 700.0) ei = 700.0;
                w[i] = std::exp(ei);
            }
        }
        return weighted_crossprod(m_X, w);
    }
};

ModelResult fast_poisson_internal(const Eigen::MatrixXd& X,
							 const Eigen::VectorXd& y,
                             const Eigen::VectorXd& weights = Eigen::VectorXd(),
                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                             bool smart_cold_start = false,
							 int maxit = 100,
							 double tol = 1e-8,
                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                             std::string optimization_alg = "lbfgs",
                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                             Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                             bool estimate_only = false) {
	const int n = X.rows();
	const int p = X.cols();
    bool use_weights = (weights.size() == n);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    std::string alg = normalize_optimizer_algorithm(optimization_alg, "lbfgs", true);
    VectorXd beta_start = VectorXd::Zero(p);
    if (warm_start_beta.isNotNull()) {
        beta_start = as<Eigen::VectorXd>(Rcpp::NumericVector(warm_start_beta));
        if (beta_start.size() != p) Rcpp::stop("warm_start_beta must have length equal to ncol(X)");
    } else if (smart_cold_start) {
        beta_start = edi_opt::poisson_smart_cold_start(X, y);
    }
    beta_start = apply_fixed_values(beta_start, fixed_spec);

    if (alg != "irls") {
        VectorXd beta = beta_start;
        PoissonNegLogLik fun(X, y, weights);
        
        Eigen::MatrixXd H_start_val;
        const Eigen::MatrixXd* h_ptr = nullptr;
        if (warm_start_fisher_info.isNotNull()) {
            H_start_val = as<Eigen::MatrixXd>(warm_start_fisher_info);
            h_ptr = &H_start_val;
        } else if (smart_cold_start) {
            H_start_val = edi_opt::poisson_smart_hessian(X, beta_start);
            h_ptr = &H_start_val;
        }

        LikelihoodFitResult fit = (alg == "lbfgs") ?
            optimize_fixed_likelihood_lbfgs(fun, beta, fixed_spec, maxit, tol) :
            optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, alg, "newton_raphson", 0, h_ptr);

        ModelResult res;
        res.b = fit.params;
        if (!estimate_only) {
            VectorXd eta = X * res.b;
            eta = eta.cwiseMin(700.0);
            res.mu = eta.array().exp().matrix();
            res.XtWX = fun.hessian(res.b);
        }
        res.iterations = fit.niter;
        res.converged = fit.converged;
        return res;
    }

    const int p_free = fixed_spec.free_idx.size();
    RowMajorMatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
    VectorXd beta_free(p_free);
    for (int j = 0; j < p_free; ++j) beta_free[j] = beta_start[fixed_spec.free_idx[j]];
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
        res.iterations = iter + 1;
        eta.noalias() = X_free * beta_free;
        eta += eta_fixed;
        
        // Combined clamping and exp to avoid multiple passes and temps
        for (int i = 0; i < n; ++i) {
            double ei = eta[i];
            if (ei > 700.0) ei = 700.0;
            double mui = std::exp(ei);
            if (mui < 1e-10) mui = 1e-10;
            mu[i] = mui;
            
            // Compute z and w in the same pass
            double wi = mui;
            if (use_weights) wi *= weights[i];
            w[i] = wi;
            
            z_adj[i] = (y[i] - mui) / mui + ei - eta_fixed[i];
        }

        if (iter == 0 && warm_start_weights.isNotNull()) {
            Eigen::VectorXd ww = as<Eigen::VectorXd>(warm_start_weights);
            if (ww.size() == n) w = ww;
        }

        if (iter == 0 && warm_start_fisher_info.isNotNull()) {
            Eigen::MatrixXd info_full = as<Eigen::MatrixXd>(warm_start_fisher_info);
            if (info_full.rows() == p && info_full.cols() == p) {
                XtWX_free = subset_matrix(info_full, fixed_spec.free_idx, fixed_spec.free_idx);
            } else {
                XtWX_free = weighted_crossprod(X_free, w);
            }
        } else if (iter == 0 && smart_cold_start && warm_start_beta.isNull()) {
            Eigen::MatrixXd H_full = edi_opt::poisson_smart_hessian(X, beta_start);
            XtWX_free = subset_matrix(H_full, fixed_spec.free_idx, fixed_spec.free_idx);
        } else {
            XtWX_free = weighted_crossprod(X_free, w);
        }
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

    if (!estimate_only) {
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
    }

	return res;
}

} // namespace

ModelResult fast_poisson_regression_internal(const Eigen::MatrixXd& X,
                                             const Eigen::VectorXd& y,
                                             const Eigen::VectorXd& weights = Eigen::VectorXd(),
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                             bool smart_cold_start = false,
                                             int maxit = 100,
                                             double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "irls",
                                             Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                             bool estimate_only = false) {
    return fast_poisson_internal(X, y, weights, warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
}

//' @title Compute Poisson Regression Score
//' @description Calculates the score vector (gradient of the log-likelihood) for a Poisson regression model.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param beta A numeric vector of coefficients.
//' @return A numeric vector representing the score.
//' @export
//' @keywords internal
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
//' @keywords internal
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
//' @keywords internal
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
//' @keywords internal
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
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm ("lbfgs" or "irls").
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, fitted values, and information matrix.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rpois(10, 2)
//' fast_poisson_regression_cpp(X, y)
// [[Rcpp::export]]
List fast_poisson_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y,
                                     Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                     bool smart_cold_start = false,
									 int maxit = 100,
									 double tol = 1e-8,
                                     Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                     std::string optimization_alg = "lbfgs",
                                     Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
                                     bool estimate_only = false) {
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
	if (estimate_only) {
		return List::create(
			Named("b") = res.b,
			Named("converged") = res.converged,
			Named("iterations") = res.iterations
		);
	}
	return List::create(
		Named("b") = res.b,
		Named("mu") = res.mu,
		Named("XtWX") = res.XtWX,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations
	);
}

//' @title Fast Weighted Poisson Regression (C++)
//' @description High-performance weighted Poisson regression fitting.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param weights A numeric vector of weights.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, fitted values, and information matrix.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_poisson_regression_weighted_cpp(const Eigen::MatrixXd& X,
                                          const Eigen::VectorXd& y,
                                          const Eigen::VectorXd& weights,
                                          Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                          bool smart_cold_start = false,
                                          int maxit = 100,
                                          double tol = 1e-8,
                                          Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                          Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                          std::string optimization_alg = "irls",
                                          Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                          Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
	ModelResult res = fast_poisson_internal(X, y, weights, warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
	return List::create(
		Named("b") = res.b,
		Named("mu") = res.mu,
		Named("XtWX") = res.XtWX,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations
	);
}

//' @title Fast Poisson Regression with Variance (C++)
//' @description Poisson regression with variance-covariance matrix and score calculation.
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, vcov, score, and likelihood statistics.
//' @export
//' @keywords internal
//' @examples
//' X = matrix(rnorm(100), 10, 10)
//' y = rpois(10, 2)
//' fast_poisson_regression_with_var_cpp(X, y)
// [[Rcpp::export]]
List fast_poisson_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int j = 2,
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                              bool smart_cold_start = false,
											  int maxit = 100,
											  double tol = 1e-8,
                                              Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                              std::string optimization_alg = "irls",
                                              Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                              Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int k) -> int {
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == k) return jj + 1;
        return -1;
    };
    int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
    res.ssq_b_j = (free_j > 0) ? compute_diagonal_inverse_entry(info_free, free_j) : NA_REAL;
    int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;
    res.ssq_b_2 = (free_2 > 0) ? compute_diagonal_inverse_entry(info_free, free_2) : NA_REAL;

	VectorXd score = get_poisson_regression_score_cpp(X, y, res.b);
	MatrixXd information = res.XtWX;
	Eigen::ArrayXd eta = (X * res.b).array().max(-30.0).min(30.0);
	double neg_loglik = (eta.exp() - y.array() * eta).sum();

	return List::create(
		Named("b") = res.b,
		Named("params") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("mu") = res.mu,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations,
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
//' @param X A numeric matrix of predictors.
//' @param y A numeric vector of responses (counts).
//' @param j 1-based index of the parameter for which to return specific variance.
//' @param warm_start_beta Optional starting values for coefficients. If provided, \code{smart_cold_start} is ignored.
//' @param smart_cold_start Logical. If TRUE, use an initial OLS-based guess when starting from scratch (a "cold start") with no prior knowledge. This is ignored if a warm start is provided.
//' @param maxit Maximum number of iterations.
//' @param tol Convergence tolerance.
//' @param fixed_idx Optional indices of fixed parameters.
//' @param fixed_values Optional values for fixed parameters.
//' @param optimization_alg Optimization algorithm.
//' @param warm_start_weights Optional initial working weights for the first IRLS iteration.
//' @param warm_start_fisher_info Optional initial Fisher Information matrix for the first IRLS iteration.
//' @return A list containing coefficients, vcov, and dispersion estimate.
//' @export
//' @keywords internal
// [[Rcpp::export]]
List fast_quasipoisson_regression_with_var_cpp(const Eigen::MatrixXd& X,
												   const Eigen::VectorXd& y,
												   int j = 2,
                                                   Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
                                                   bool smart_cold_start = false,
												   int maxit = 100,
												   double tol = 1e-8,
                                                   Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                   Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                   std::string optimization_alg = "irls",
                                                   Rcpp::Nullable<Rcpp::NumericVector> warm_start_weights = R_NilValue,
                                                   Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
	ModelResult res = fast_poisson_internal(X, y, Eigen::VectorXd(), warm_start_beta, smart_cold_start, maxit, tol, fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    MatrixXd info_free = subset_matrix(res.XtWX, fixed_spec.free_idx, fixed_spec.free_idx);

    auto free_idx_of = [&](int k) -> int {
        for (int jj = 0; jj < (int)fixed_spec.free_idx.size(); ++jj)
            if (fixed_spec.free_idx[jj] == k) return jj + 1;
        return -1;
    };

	const int df_resid = X.rows() - X.cols();
	if (df_resid > 0) {
        double pearson_sum = 0.0;
        int n = X.rows();
        for (int i = 0; i < n; ++i) {
            double diff = y[i] - res.mu[i];
            pearson_sum += (diff * diff) / res.mu[i];
        }
		res.dispersion = pearson_sum / static_cast<double>(df_resid);

		if (std::isfinite(res.dispersion) && res.dispersion > 0) {
            int free_j = (j > 0 && j <= X.cols()) ? free_idx_of(j - 1) : -1;
            int free_2 = (X.cols() >= 2) ? free_idx_of(1) : -1;

            if (free_j > 0) {
                res.ssq_b_j = res.dispersion * compute_diagonal_inverse_entry(info_free, free_j);
            } else {
                res.ssq_b_j = NA_REAL;
            }

            if (free_2 > 0) {
                if (free_2 == free_j) {
                    res.ssq_b_2 = res.ssq_b_j;
                } else {
                    res.ssq_b_2 = res.dispersion * compute_diagonal_inverse_entry(info_free, free_2);
                }
            }
		}
	}

	return List::create(
		Named("b") = res.b,
		Named("ssq_b_j") = res.ssq_b_j,
		Named("ssq_b_2") = res.ssq_b_2,
		Named("dispersion") = res.dispersion,
		Named("mu") = res.mu,
		Named("converged") = res.converged,
		Named("iterations") = res.iterations
	);
}

//' Parallel Poisson Randomization Distribution
//'
//' @param X Matrix of covariates (without intercept or treatment).
//' @param y Numeric vector of response values (pre-null-shifted for treated).
//' @param w_mat Integer matrix of permuted treatment assignments (n x nsim).
//' @param delta Null treatment effect shift.
//' @param log_transform If TRUE, apply multiplicative delta shift (exp scale); otherwise additive.
//' @param num_cores Number of OpenMP threads.
//' @return Numeric vector of length nsim with treatment coefficients.
// [[Rcpp::export]]
NumericVector compute_poisson_distr_parallel_cpp(
	const Eigen::MatrixXd& X,
	const Eigen::VectorXd& y,
	const Rcpp::IntegerMatrix& w_mat,
	double delta,
	bool log_transform,
	int num_cores
) {
	int nsim = w_mat.cols();
	int n = y.size();
	int p_covars = X.cols();
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
			ws[t].X_full.col(2 + k) = X.col(k);
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

		try {
			ModelResult res = fast_poisson_internal(W.X_full, W.y_shifted, Eigen::VectorXd());
			if (res.converged && p_full >= 2 && std::isfinite(res.b[1])) {
				results[b] = res.b[1];
			}
		} catch (...) {
			// Leave results[b] = NA_REAL on any optimizer failure
		}
	}

	return wrap(results);
}
