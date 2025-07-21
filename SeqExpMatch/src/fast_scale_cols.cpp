// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat scale_columns_cpp(const arma::mat& X) {
  arma::mat scaled = X;
  arma::uword n_cols = X.n_cols;

  for (arma::uword j = 0; j < n_cols; ++j) {
    arma::vec col = X.col(j);

    if (!col.is_finite()) {
      stop("Column %d contains NA/NaN/Inf values", j + 1);
    }

    double mu = arma::mean(col);
    double sigma = arma::stddev(col, 0);  // sample std dev

    if (std::isfinite(sigma) && sigma > 0) {
      scaled.col(j) = (col - mu) / sigma;
    } else {
      scaled.col(j).zeros();  // constant or invalid column
    }
  }

  return scaled;
}