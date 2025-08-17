#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd compute_proportional_mahal_distances_cpp(
    const Eigen::VectorXd& xt_prev,              // fixed vector
    const Eigen::MatrixXd& X_prev,               // full data matrix
    const Eigen::VectorXi& reservoir_indices,    // indices (1-based from R)
    const Eigen::MatrixXd& S_xs_inv)             // inverse covariance matrix
{
    int n_R = reservoir_indices.size();
    Eigen::VectorXd sqd_distances_times_two(n_R);

    for (int r = 0; r < n_R; r++) {
        // Convert from R's 1-based to C++'s 0-based indexing
        int idx = reservoir_indices[r] - 1;
        
        // Difference vector
        Eigen::VectorXd diff = xt_prev - X_prev.row(idx).transpose();
        
        // Quadratic form: diff' * S_xs_inv * diff
        sqd_distances_times_two[r] = diff.transpose() * S_xs_inv * diff;
    }
    
    return sqd_distances_times_two;
}