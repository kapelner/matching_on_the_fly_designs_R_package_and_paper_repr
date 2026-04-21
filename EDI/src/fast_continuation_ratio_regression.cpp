#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

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
List fast_continuation_ratio_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8) {
    int n = X.rows();
    int p = X.cols();
    
    // Determine unique levels of y
    std::vector<double> levels;
    for (int i = 0; i < y.size(); ++i) {
        if (std::find(levels.begin(), levels.end(), y[i]) == levels.end()) {
            levels.push_back(y[i]);
        }
    }
    std::sort(levels.begin(), levels.end());
    int K = levels.size();
    if (K < 2) {
        return List::create(Named("b") = VectorXd::Zero(p), Named("alpha") = VectorXd::Zero(0));
    }
    int n_alpha = K - 1; // Number of continuation trials
    
    // Data augmentation
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
        
        // Trial j: 1 if they stopped at category levels[j], 0 if they went further
        // If they reached levels[yi_level], they had to pass through 0, 1, ..., yi_level-1
        for (int j = 0; j < std::min(yi_level + 1, n_alpha); ++j) {
            VectorXd x_row(n_alpha + p);
            x_row.setZero();
            x_row[j] = 1.0; // Intercept for category j
            if (p > 0) {
                x_row.tail(p) = X.row(i);
            }
            
            x_aug_vecs.push_back(x_row);
            z_vals.push_back((yi_level == j) ? 1.0 : 0.0);
        }
    }
    
    int n_aug = z_vals.size();
    int p_aug = n_alpha + p;
    
    MatrixXd X_aug(n_aug, p_aug);
    VectorXd z(n_aug);
    for (int i = 0; i < n_aug; ++i) {
        X_aug.row(i) = x_aug_vecs[i];
        z[i] = z_vals[i];
    }
    
    // Fit logistic regression on augmented data using IRWLS
    VectorXd beta = VectorXd::Zero(p_aug);
    bool converged = false;
    
    for (int iter = 0; iter < maxit; ++iter) {
        VectorXd eta = X_aug * beta;
        VectorXd mu = eta.unaryExpr([](double x) { 
            if (x > 20) return 1.0;
            if (x < -20) return 0.0;
            return 1.0 / (1.0 + std::exp(-x)); 
        });
        VectorXd w = mu.array() * (1.0 - mu.array());
        
        for (int i = 0; i < n_aug; ++i) if (w[i] < 1e-10) w[i] = 1e-10;
        
        VectorXd grad = X_aug.transpose() * (z - mu);
        MatrixXd hess = X_aug.transpose() * w.asDiagonal() * X_aug;
        
        FullPivLU<MatrixXd> lu(hess);
        if (!lu.isInvertible()) break;
        
        VectorXd step = lu.solve(grad);
        beta += step;
        
        if (step.norm() < tol) {
            converged = true;
            break;
        }
    }
    
    return List::create(
        Named("b") = beta.tail(p),
        Named("alpha") = beta.head(n_alpha),
        Named("beta_full") = beta,
        Named("X_aug") = X_aug,
        Named("z") = z,
        Named("converged") = converged
    );
}

// [[Rcpp::export]]
List fast_continuation_ratio_regression_with_var_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) {
    List res = fast_continuation_ratio_regression_cpp(X, y);
    if (!res.containsElementNamed("beta_full")) {
         return List::create(Named("b") = NumericVector::create(NA_REAL), Named("ssq_b_2") = NA_REAL, Named("converged") = false);
    }
    VectorXd beta_full = res["beta_full"];
    MatrixXd X_aug = res["X_aug"];
    bool converged = res["converged"];
    
    VectorXd eta = X_aug * beta_full;
    VectorXd mu = eta.unaryExpr([](double x) { 
        if (x > 20) return 1.0;
        if (x < -20) return 0.0;
        return 1.0 / (1.0 + std::exp(-x)); 
    });
    VectorXd w = mu.array() * (1.0 - mu.array());
    for (int i = 0; i < w.size(); ++i) if (w[i] < 1e-10) w[i] = 1e-10;
    
    MatrixXd XtWX = X_aug.transpose() * w.asDiagonal() * X_aug;
    FullPivLU<MatrixXd> lu(XtWX);
    
    double ssq_b_j = NA_REAL;
    if (lu.isInvertible()) {
        MatrixXd inv = lu.inverse();
        int n_alpha = (as<VectorXd>(res["alpha"])).size();
        if (X.cols() >= 1) {
            ssq_b_j = inv(n_alpha, n_alpha);
        }
    }
    
    return List::create(
        Named("b") = res["b"],
        Named("ssq_b_2") = ssq_b_j,
        Named("converged") = converged
    );
}
