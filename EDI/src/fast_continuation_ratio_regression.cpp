#include "_helper_functions.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

namespace {

inline Eigen::ArrayXd plogis_array_clamped(const Eigen::ArrayXd& eta) {
    const Eigen::ArrayXd eta_clamped = eta.max(-20.0).min(20.0);
    return 1.0 / (1.0 + (-eta_clamped).exp());
}

struct ContinuationRatioObjective {
    const MatrixXd& X_aug;
    const VectorXd& z;

    ContinuationRatioObjective(const MatrixXd& X_aug, const VectorXd& z) :
        X_aug(X_aug), z(z) {}

    double operator()(const VectorXd& beta, VectorXd& grad) const {
        VectorXd eta = X_aug * beta;
        VectorXd mu = plogis_array_clamped(eta.array()).matrix();
        grad = X_aug.transpose() * (mu - z); // Negative log-likelihood gradient

        const Eigen::ArrayXd mu_arr = mu.array();
        const Eigen::ArrayXd log_mu = mu_arr.max(1e-12).log();
        const Eigen::ArrayXd log_one_minus_mu = (1.0 - mu_arr).max(1e-12).log();
        return -(z.array() * log_mu + (1.0 - z.array()) * log_one_minus_mu).sum();
    }

    MatrixXd hessian(const VectorXd& beta) const {
        VectorXd eta = X_aug * beta;
        VectorXd mu = plogis_array_clamped(eta.array()).matrix();
        VectorXd w = mu.array() * (1.0 - mu.array());
        return weighted_crossprod(X_aug, w);
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
	VectorXd mu = plogis_array_clamped(eta.array()).matrix();
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
	VectorXd mu = plogis_array_clamped(eta.array()).matrix();
	VectorXd w = mu.array() * (1.0 - mu.array());
	return -weighted_crossprod(X_aug, w);
}

// [[Rcpp::export]]
List fast_continuation_ratio_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8,
                                             Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
                                             Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue,
                                             std::string optimization_alg = "newton_raphson",
                                             Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
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
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, optimization_alg, "newton_raphson", 0, info_start_ptr);
    
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
                                                      std::string optimization_alg = "newton_raphson",
                                                      Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue) {
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
    
    Eigen::MatrixXd info_start;
    const Eigen::MatrixXd* info_start_ptr = nullptr;
    if (warm_start_fisher_info.isNotNull()) {
        info_start = as<Eigen::MatrixXd>(warm_start_fisher_info);
        info_start_ptr = &info_start;
    }

    LikelihoodFitResult fit = optimize_fixed_likelihood(fun, beta, fixed_spec, maxit, tol, optimization_alg, "newton_raphson", 0, info_start_ptr);
    
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
