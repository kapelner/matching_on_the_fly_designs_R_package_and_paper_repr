#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd scale_columns_cpp(const Eigen::MatrixXd& X) {
  Eigen::Index n = X.rows();
  Eigen::Index p = X.cols();
  
  Eigen::MatrixXd X_scaled(n, p);
  
  for (Eigen::Index j = 0; j < p; ++j) {
    double mean = X.col(j).mean();
    double sd = std::sqrt((X.col(j).array() - mean).square().sum() / (n - 1));
    
    if (sd > 0) {
      X_scaled.col(j) = (X.col(j).array() - mean) / sd;
    } else {
      // If standard deviation is 0, just subtract mean (constant column)
      X_scaled.col(j) = X.col(j).array() - mean;
    }
  }

  return X_scaled;
}
