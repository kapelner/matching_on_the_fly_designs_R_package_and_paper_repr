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

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        Eigen::VectorXi w = w_mat.col(b);
        
        // Construct the design matrix
        Eigen::MatrixXd X_full(n, p_full);
        X_full.col(0).setOnes();
        X_full.col(1) = w.cast<double>();
        if (p_covars > 0) {
            X_full.rightCols(p_covars) = X_covars;
        }
        
        // OLS via LDLT
        Eigen::MatrixXd XtX = X_full.transpose() * X_full;
        Eigen::VectorXd Xty = X_full.transpose() * y;
        
        Eigen::VectorXd beta = XtX.ldlt().solve(Xty);
        results[b] = beta[1]; // Treatment effect is the 2nd coefficient
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
