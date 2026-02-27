#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector compute_simple_mean_diff_parallel_cpp(
    const Eigen::MatrixXd& y_mat,
    const Eigen::MatrixXi& w_mat,
    int num_cores) {
    
    int nsim = w_mat.cols();
    int n = w_mat.rows();
    NumericVector results(nsim);
    
#ifdef _OPENMP
    if (num_cores > 1) {
        omp_set_num_threads(num_cores);
    }
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        Eigen::VectorXd y = y_mat.col(b);
        Eigen::VectorXi w = w_mat.col(b);
        
        double sum_T = 0, sum_C = 0;
        int n_T = 0, n_C = 0;
        
        for (int i = 0; i < n; ++i) {
            if (w[i] == 1) {
                sum_T += y[i];
                n_T++;
            } else {
                sum_C += y[i];
                n_C++;
            }
        }
        
        if (n_T == 0 || n_C == 0) {
            results[b] = NA_REAL;
        } else {
            results[b] = (sum_T / n_T) - (sum_C / n_C);
        }
    }
    
    return results;
}
