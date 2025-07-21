// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
List fast_logistic_regression_cpp(const arma::mat& X, const arma::vec& y, 
                   int max_iter = 25, double tol = 1e-8) {
  int n = X.n_rows;
  int p = X.n_cols;

  arma::vec beta = arma::zeros(p);
  arma::vec eta(n), mu(n), W_diag(n), z(n);
  arma::mat X_weighted(n, p);
  arma::mat XtWX;
  bool converged = false;

  for (int iter = 0; iter < max_iter; ++iter) {
    eta = X * beta;
    mu = 1.0 / (1.0 + arma::exp(-eta));
    mu = arma::clamp(mu, 1e-8, 1 - 1e-8);

    W_diag = mu % (1 - mu);
    z = eta + (y - mu) / W_diag;

    // X_weighted = each column of X multiplied by sqrt(W_diag)
    for (int j = 0; j < p; ++j)
      X_weighted.col(j) = X.col(j) % arma::sqrt(W_diag);

    arma::vec z_weighted = z % arma::sqrt(W_diag);

    // Solve using Cholesky
    XtWX = X_weighted.t() * X_weighted;
    arma::vec Xtz = X_weighted.t() * z_weighted;
    arma::vec beta_new = arma::solve(XtWX, Xtz, arma::solve_opts::fast);

    if (arma::norm(beta_new - beta, 2) < tol) {
      beta = beta_new;
      converged = true;
      break;
    }

    beta = beta_new;
  }

  return List::create(
    Named("b")  			= beta,
    Named("converged")     	= converged,
    Named("iterations")    	= max_iter,
    Named("XtWX") 			= XtWX
  );
}

// [[Rcpp::export]]
List fast_logistic_regression_with_sd_cpp(const arma::mat& X, const arma::vec& y, 
                   int max_iter = 25, double tol = 1e-8) {
  List mod = fast_logistic_regression_cpp(X, y, max_iter, tol);

  // Standard errors from inverse Hessian (XtWX)
  arma::mat XtWX = mod["XtWX"];
  arma::mat cov = arma::inv_sympd(XtWX);
  mod["s_b"] = arma::sqrt(cov.diag());

  return mod;
}
