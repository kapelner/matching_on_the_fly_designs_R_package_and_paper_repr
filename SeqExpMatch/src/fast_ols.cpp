// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
List fast_ols_cpp(const arma::mat& X, const arma::vec& y) {
  arma::mat XtX = X.t() * X;
  arma::vec Xty = X.t() * y;
  return List::create(
    Named("b") = solve(XtX, Xty),
    Named("XtX") = XtX
  );
}

// [[Rcpp::export]]
List fast_ols_with_sd_cpp(const arma::mat& X, const arma::vec& y) {
  int n = X.n_rows;
  int p = X.n_cols;
  List mod = fast_ols_cpp(X, y);
  arma::vec beta = mod["b"];

  // Residuals and SSE
  arma::vec residuals = y - X * beta;
  double rss = arma::dot(residuals, residuals); // Residual Sum of Squares
  double sigma2_hat = rss / (n - p);            // Residual variance estimate

  // Variance-covariance matrix of beta
  arma::mat XtX = mod["XtX"];
  arma::mat XtX_inv = inv_sympd(XtX);           // Use inv_sympd for speed/safety
  mod["s_b"] = arma::sqrt(sigma2_hat * XtX_inv.diag());

  return mod;
}
