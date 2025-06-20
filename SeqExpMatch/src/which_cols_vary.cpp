#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector which_cols_vary_cpp(NumericMatrix X){
  int n = X.nrow();
  int p = X.ncol();
  LogicalVector which_cols_vary(p);
  for (int j = 0; j < p; j++){
	  for (int i = 1; i < n; i++){
	    if (X(i, j) != X(0, j)){
	      which_cols_vary[j] = TRUE;
	      break;
	    }
	  }
  }
  return which_cols_vary;
}
