#include "_helper_functions.h"
#include <RcppEigen.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// --- Weight Functions ---

// Huber weight function
double huber_w(double r, double c) {
    double abs_r = std::abs(r);
    if (abs_r <= c) return 1.0;
    return c / abs_r;
}

// Tukey's Bisquare (Biweight) weight function
double bisquare_w(double r, double c) {
    double abs_r = std::abs(r);
    if (abs_r <= c) {
        double tmp = 1.0 - (r / c) * (r / c);
        return tmp * tmp;
    }
    return 0.0;
}

// --- Internal IRLS Logic ---

struct RobustModelResult {
    Eigen::VectorXd b;
    Eigen::VectorXd w;
    Eigen::MatrixXd XtWX;
    double scale;
    int iterations;
    bool converged;
    double ssq_b_j;

    RobustModelResult() : scale(NA_REAL), iterations(0), converged(false), ssq_b_j(NA_REAL) {}
};

RobustModelResult fast_robust_regression_internal(
    const Eigen::MatrixXd& X, 
    const Eigen::VectorXd& y, 
    Nullable<NumericVector> start_beta = R_NilValue,
    std::string method = "MM",
    double c = 1.345, // Huber constant
    double c_bisquare = 4.685, // Bisquare constant
    int maxit = 50,
    double tol = 1e-7,
    double scale_est = -1.0, // If negative, compute MAD
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
    int n = X.rows();
    int p = X.cols();
    FixedParamSpec fixed_spec = make_fixed_param_spec(p, fixed_idx, fixed_values);
    const int p_free = fixed_spec.free_idx.size();
    Eigen::MatrixXd X_free(n, p_free);
    for (int j = 0; j < p_free; ++j) X_free.col(j) = X.col(fixed_spec.free_idx[j]);
    Eigen::VectorXd y_adj = y;
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) {
        y_adj.noalias() -= X.col(fixed_spec.fixed_idx[j]) * fixed_spec.fixed_values[j];
    }
    RobustModelResult res;

    // 1. Initial estimate
    Eigen::VectorXd b_free;
    if (start_beta.isNotNull()) {
        Eigen::VectorXd b_start = as<Eigen::VectorXd>(NumericVector(start_beta));
        b_start = apply_fixed_values(b_start, fixed_spec);
        b_free = subset_vector(b_start, fixed_spec.free_idx);
    } else {
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(X_free);
        b_free = cod.solve(y_adj);
    }
    res.b = Eigen::VectorXd::Zero(p);
    for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = b_free[j];
    for (int j = 0; j < fixed_spec.fixed_idx.size(); ++j) res.b[fixed_spec.fixed_idx[j]] = fixed_spec.fixed_values[j];

    Eigen::VectorXd r = y - X * res.b;
    
    // 2. Scale estimation (MAD of residuals)
    if (scale_est < 0) {
        std::vector<double> abs_r(n);
        for (int i = 0; i < n; ++i) abs_r[i] = std::abs(r[i]);
        std::sort(abs_r.begin(), abs_r.end());
        double median_abs_r = (n % 2 == 0) ? (abs_r[n / 2 - 1] + abs_r[n / 2]) / 2.0 : abs_r[n / 2];
        res.scale = median_abs_r / 0.6745;
    } else {
        res.scale = scale_est;
    }

    if (res.scale < 1e-10) res.scale = 1e-10;

    // 3. IRLS loop
    res.w = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd b_old = b_free;
    
    for (int iter = 1; iter <= maxit; ++iter) {
        res.iterations = iter;
        
        // Update weights
        for (int i = 0; i < n; ++i) {
            double u = r[i] / res.scale;
            if (method == "M") {
                res.w[i] = huber_w(u, c);
            } else {
                res.w[i] = bisquare_w(u, c_bisquare);
            }
        }

        // Solve Weighted Least Squares
        Eigen::MatrixXd XtWX = X_free.transpose() * res.w.asDiagonal() * X_free;
        Eigen::VectorXd XtWy = X_free.transpose() * res.w.asDiagonal() * y_adj;
        
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
        if (ldlt.info() != Eigen::Success) break; // Numerical failure
        
        b_free = ldlt.solve(XtWy);
        for (int j = 0; j < p_free; ++j) res.b[fixed_spec.free_idx[j]] = b_free[j];
        r = y - X * res.b;

        // Check convergence
        if ((b_free - b_old).norm() / (b_free.norm() + 1e-10) < tol) {
            res.converged = true;
            res.XtWX = expand_free_covariance(p, fixed_spec, XtWX, false);
            break;
        }
        b_old = b_free;
    }

    return res;
}

// [[Rcpp::export]]
List fast_robust_regression_cpp(
    const Eigen::MatrixXd& X, 
    const Eigen::VectorXd& y, 
    Nullable<NumericVector> start_beta = R_NilValue,
    std::string method = "MM",
    int j = 2,
    double c = 1.345,
    int maxit = 50,
    double tol = 1e-7,
    Rcpp::Nullable<Rcpp::IntegerVector> fixed_idx = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> fixed_values = R_NilValue
) {
    RobustModelResult res = fast_robust_regression_internal(X, y, start_beta, method, c, 4.685, maxit, tol, -1.0, fixed_idx, fixed_values);
    FixedParamSpec fixed_spec = make_fixed_param_spec(X.cols(), fixed_idx, fixed_values);
    
    if (!res.converged && res.XtWX.rows() == 0) {
        Eigen::MatrixXd X_free(X.rows(), fixed_spec.free_idx.size());
        for (int jj = 0; jj < fixed_spec.free_idx.size(); ++jj) X_free.col(jj) = X.col(fixed_spec.free_idx[jj]);
        Eigen::MatrixXd XtWX_free = X_free.transpose() * res.w.asDiagonal() * X_free;
        res.XtWX = expand_free_covariance(X.cols(), fixed_spec, XtWX_free, false);
    }

    int n = X.rows();
    int p = X.cols();
    Eigen::VectorXd r = y - X * res.b;
    
    double ssq_j = NA_REAL;
    if (res.converged || res.iterations == maxit) {
        Eigen::VectorXd psi_r(n);
        double sum_psi_prime = 0;
        
        for (int i = 0; i < n; ++i) {
            double u = r[i] / res.scale;
            double abs_u = std::abs(u);
            if (method == "M") {
                if (abs_u <= c) {
                    psi_r[i] = r[i];
                    sum_psi_prime += 1.0;
                } else {
                    psi_r[i] = c * res.scale * (u > 0 ? 1.0 : -1.0);
                    sum_psi_prime += 0.0;
                }
            } else {
                double c_b = 4.685;
                if (abs_u <= c_b) {
                    double tmp = 1.0 - (u/c_b)*(u/c_b);
                    psi_r[i] = r[i] * tmp * tmp;
                    sum_psi_prime += tmp * (1.0 - 5.0 * (u/c_b)*(u/c_b));
                } else {
                    psi_r[i] = 0;
                    sum_psi_prime += 0;
                }
            }
        }
        
        double m = sum_psi_prime / n;
        Eigen::MatrixXd X_free(n, fixed_spec.free_idx.size());
        for (int jj = 0; jj < fixed_spec.free_idx.size(); ++jj) X_free.col(jj) = X.col(fixed_spec.free_idx[jj]);
        Eigen::MatrixXd XtX = X_free.transpose() * X_free;
        Eigen::MatrixXd invXtX_free = XtX.inverse();
        Eigen::MatrixXd vcov_free;
        
        double sum_psi_sq = psi_r.squaredNorm();
        double factor = (n / (double(n - fixed_spec.free_idx.size()))) * sum_psi_sq / (n * m * m);
        vcov_free = factor * invXtX_free;
        Eigen::MatrixXd vcov = expand_free_covariance(p, fixed_spec, vcov_free, true);
        
        if (j > 0 && j <= p) {
            ssq_j = vcov(j-1, j-1);
        }
    }

    return List::create(
        Named("coefficients") = res.b,
        Named("scale") = res.scale,
        Named("converged") = res.converged,
        Named("iterations") = res.iterations,
        Named("ssq_b_j") = ssq_j
    );
}
