#include "_helper_functions.h"
#include <RcppEigen.h>

using namespace Rcpp;

// [[Rcpp::export]]
List fast_logistic_regression_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, int maxit = 100, double tol = 1e-8) {
  int n = X.rows();
  int p = X.cols();
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd prob(n);
  Eigen::VectorXd w(n);
  Eigen::VectorXd z(n);
  
  for (int iter = 0; iter < maxit; iter++) {
    Eigen::VectorXd eta = X * beta;
    // prob = 1 / (1 + exp(-eta))
    // Avoid overflow/underflow
    for(int i=0; i<n; i++) {
        double val = 1.0 / (1.0 + std::exp(-eta[i]));
        prob[i] = val;
        // Weights w = p * (1-p)
        w[i] = std::max(val * (1.0 - val), 1e-10); 
    }
    
    // IRLS step
    // z = eta + (y - prob) / w
    z = eta + (y - prob).cwiseQuotient(w);
    
    // beta_new = (X'WX)^-1 X'W z
    Eigen::MatrixXd XtW = X.transpose() * w.asDiagonal();
    Eigen::MatrixXd XtWX = XtW * X;
    Eigen::VectorXd XtWz = XtW * z;
    
    Eigen::VectorXd beta_new = XtWX.ldlt().solve(XtWz);
    
    if ((beta - beta_new).norm() < tol) {
      beta = beta_new;
      break;
    }
    beta = beta_new;
  }
  
  // Recalculate weights for final output if needed, but we have them from last iter.
  // Actually, we should make sure w matches beta.
  Eigen::VectorXd eta = X * beta;
  for(int i=0; i<n; i++) {
      double val = 1.0 / (1.0 + std::exp(-eta[i]));
      w[i] = std::max(val * (1.0 - val), 1e-10);
  }

  return List::create(
    Named("b") = beta,
    Named("w") = w
  );
}

//' Fast Logistic Regression with variance using Rcpp and IRLS
//'
//' @param Xmm Design matrix.
//' @param y Response vector.
//' @return A list with coefficients and specific treatment variance.
//' @export
// [[Rcpp::export]]
List fast_logistic_regression_with_var_cpp(const Eigen::MatrixXd& Xmm, const Eigen::VectorXd& y) {
  List mod = fast_logistic_regression_cpp(Xmm, y);
  Eigen::VectorXd w = mod["w"];
  Eigen::VectorXd b = mod["b"];
  
  // Xmm.transpose() * w.asDiagonal() * Xmm
  Eigen::MatrixXd XtWX = Xmm.transpose() * w.asDiagonal() * Xmm;
  
  double ssq_b_2 = eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp(XtWX, 2);
  
  return List::create(
    Named("b") = b,
    Named("ssq_b_2") = ssq_b_2
  );
}
