#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector sample_int_replace_cpp(int n, int size) {
  IntegerVector result(size);
  for (int i = 0; i < size; ++i) {
    result[i] = 1 + (int)(R::unif_rand() * n); // 1 to n
  }
  return result;
}
