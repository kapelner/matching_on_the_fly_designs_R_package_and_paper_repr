#include <RcppEigen.h>
#include <cmath>    // for std::log, std::exp, std::lgamma

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double neg_loglik_nb_cpp(double theta,
							const Eigen::VectorXd &beta,
			                const Eigen::MatrixXd &X,
			                const Eigen::VectorXi &y
	                 	) {
  
  Eigen::VectorXd eta = X * beta;
  Eigen::VectorXd mu = eta.array().exp();
  
  double ll = 0.0;
  int n = y.size();
  
  for (int i = 0; i < n; i++) {
    double yi = static_cast<double>(y[i]);
    double mui = mu[i];
    
    ll += std::lgamma(yi + theta) 
        - std::lgamma(theta) 
        - std::lgamma(yi + 1.0)
        + theta * std::log(theta)
        + yi * std::log(mui)
        - (yi + theta) * std::log(mui + theta);
  }
  
  return -ll;
}