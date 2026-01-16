#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix bootstrap_indices_cpp(int n, int B) {
  IntegerMatrix idx(B, n);
  for (int i = 0; i < B; ++i) {
    for (int j = 0; j < n; ++j) {
      idx(i, j) = 1 + static_cast<int>(R::unif_rand() * n);
    }
  }
  return idx;
}
