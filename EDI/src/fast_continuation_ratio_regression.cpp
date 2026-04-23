#include "_helper_functions.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

struct ContinuationRatioObjective {
    const MatrixXd& X_aug;
    const VectorXd& z;

    ContinuationRatioObjective(const MatrixXd& X_aug, const VectorXd& z) :
        X_aug(X_aug), z(z) {}

    double operator()(const VectorXd& beta, VectorXd& grad) const {
        VectorXd eta = X_aug * beta;
        VectorXd mu = eta.unaryExpr([](double x) {
            if (x > 20) return 1.0;
            if (x < -20) return 0.0;
            return 1.0 / (1.0 + std::exp(-x));
        });
        grad = X_aug.transpose() * (mu - z); // Negative log-likelihood gradient
        
        double nll = 0.0;
        for (int i = 0; i < z.size(); ++i) {
            if (z[i] > 0.5) {
                nll -= std::log(std::max(1e-12, mu[i]));
            } else {
                nll -= std::log(std::max(1e-12, 1.0 - mu[i]));
            }
        }
        return nll;
    }

    MatrixXd hessian(const VectorXd& beta) const {
        VectorXd eta = X_aug * beta;
        VectorXd mu = eta.unaryExpr([](double x) {
            if (x > 20) return 1.0;
            if (x < -20) return 0.0;
            return 1.0 / (1.0 + std::exp(-x));
        });
        VectorXd w = mu.array() * (1.0 - mu.array());
        return X_aug.transpose() * w.asDiagonal() * X_aug;
    }
};

static List build_continuation_ratio_augmented_data(const Eigen::MatrixXd& X,
													const Eigen::VectorXd& y) {
	int n = X.rows();
	int p = X.cols();

	std::vector<double> levels;
	for (int i = 0; i < y.size(); ++i) {
		if (std::find(levels.begin(), levels.end(), y[i]) == levels.end()) {
			levels.push_back(y[i]);
		}
	}
	std::sort(levels.begin(), levels.end());
	int K = levels.size();
	if (K < 2) {
		return List::create(Named("X_aug") = MatrixXd(0, p), Named("z") = VectorXd(0), Named("n_alpha") = 0);
	}
	int n_alpha = K - 1;

	std::vector<double> z_vals;
	std::vector<VectorXd> x_aug_vecs;
	for (int i = 0; i < n; ++i) {
		int yi_level = -1;
		for (int k = 0; k < K; ++k) {
			if (y[i] == levels[k]) {
				yi_level = k;
				break;
			}
		}
		for (int j = 0; j < std::min(yi_level + 1, n_alpha); ++j) {
			VectorXd x_row(n_alpha + p);
			x_row.setZero();
			x_row[j] = 1.0;
			if (p > 0) x_row.tail(p) = X.row(i);
			x_aug_vecs.push_back(x_row);
			z_vals.push_back((yi_level == j) ? 1.0 : 0.0);
		}
	}

	MatrixXd X_aug(z_vals.size(), n_alpha + p);
	VectorXd z(z_vals.size());
	for (int i = 0; i < X_aug.rows(); ++i) {
		X_aug.row(i) = x_aug_vecs[i];
		z[i] = z_vals[i];
	}
	return List::create(Named("X_aug") = X_aug, Named("z") = z, Named("n_alpha") = n_alpha);
}

} // namespace

// [[Rcpp::export]]
Eigen::VectorXd get_continuation_ratio_regression_score_cpp(const Eigen::MatrixXd& X,
															const Eigen::VectorXd& y,
															const Eigen::VectorXd& params) {
	List aug = build_continuation_ratio_augmented_data(X, y);
	MatrixXd X_aug = aug["X_aug"];
	VectorXd z = aug["z"];
	if (X_aug.rows() == 0) return VectorXd::Zero(params.size());
	VectorXd eta = X_aug * params;
	VectorXd mu = eta.unaryExpr([](double x) {
		if (x > 20) return 1.0;
		if (x < -20) return 0.0;
		return 1.0 / (1.0 + std::exp(-x));
	});
	return X_aug.transpose() * (z - mu);
}

// [[Rcpp::export]]
Eigen::MatrixXd get_continuation_ratio_regression_hessian_cpp(const Eigen::MatrixXd& X,
															  const Eigen::VectorXd& y,
															  const Eigen::VectorXd& params) {
	List aug = build_continuation_ratio_augmented_data(X, y);
	MatrixXd X_aug = aug["X_aug"];
	if (X_aug.rows() == 0) return MatrixXd::Zero(params.size(), params.size());
	VectorXd eta = X_aug * params;
	VectorXd mu = eta.unaryExpr([](double x) {
		if (x > 20) return 1.0;
		if (x < -20) return 0.0;
		return 1.0 / (1.0 + std::exp(-x));
	});
	VectorXd w = mu.array() * (1.0 - mu.array());
	return -(X_aug.transpose() * w.asDiagonal() * X_aug);
}

// [[Rcpp::export]]
List fast_continuation_ratio_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "newton_raphson") {
    int p = X.cols();
    List aug = build_continuation_ratio_augmented_data(X, y);
    MatrixXd X_aug = aug["X_aug"];
    VectorXd z = aug["z"];
    int n_alpha = aug["n_alpha"];
    if (n_alpha == 0) {
        return List::create(Named("b") = VectorXd::Zero(p), Named("alpha") = VectorXd::Zero(0));
    }
    
    int p_aug = n_alpha + p;
    ContinuationRatioObjective fun(X_aug, z);
    VectorXd beta = VectorXd::Zero(p_aug);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p_aug, fixed_idx, fixed_values);
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, optimization_alg, "newton_raphson");
    
    return List::create(
        Named("b") = fit.params.tail(p),
        Named("alpha") = fit.params.head(n_alpha),
        Named("beta_full") = fit.params,
        Named("X_aug") = X_aug,
        Named("z") = z,
        Named("converged") = fit.converged
    );
}

// [[Rcpp::export]]
List fast_continuation_ratio_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8,
                                                      Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                                      Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                                      std::string optimization_alg = "newton_raphson") {
    int p = X.cols();
    List aug = build_continuation_ratio_augmented_data(X, y);
    MatrixXd X_aug = aug["X_aug"];
    VectorXd z = aug["z"];
    int n_alpha = aug["n_alpha"];
    if (n_alpha == 0) {
         return List::create(Named("b") = NumericVector::create(NA_REAL), Named("ssq_b_j") = NA_REAL, Named("converged") = false);
    }
    
    int p_aug = n_alpha + p;
    ContinuationRatioObjective fun(X_aug, z);
    VectorXd beta = VectorXd::Zero(p_aug);
    FixedParamSpec fixed_spec = make_fixed_param_spec(p_aug, fixed_idx, fixed_values);
    
    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, optimization_alg, "newton_raphson");
    
    MatrixXd info = fun.hessian(fit.params);
    MatrixXd info_free = subset_matrix(info, fixed_spec.free_idx, fixed_spec.free_idx);
    FullPivLU<MatrixXd> lu(info_free);
    
    double ssq_b_j = NA_REAL;
    MatrixXd vcov = MatrixXd::Constant(p_aug, p_aug, NA_REAL);
    if (lu.isInvertible()) {
        MatrixXd inv_free = lu.inverse();
        vcov = expand_free_covariance(p_aug, fixed_spec, inv_free, true);
        if (p >= 1) {
            ssq_b_j = vcov(n_alpha, n_alpha);
        }
    }
    
    return List::create(
        Named("b") = fit.params.tail(p),
        Named("ssq_b_j") = ssq_b_j,
        Named("vcov") = vcov,
        Named("converged") = fit.converged
    );
}
