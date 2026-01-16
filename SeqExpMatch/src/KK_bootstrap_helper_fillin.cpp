#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
void fill_i_b_with_matches_loop_cpp(IntegerVector& i_b,
                   const IntegerVector& match_indic,
                   const IntegerVector& ms_b,
                   int i_b_idx) {

  for (int m0 = 0; m0 < ms_b.size(); ++m0) {
    int target = ms_b[m0];
    int found = 0;

    for (int j = 0; j < match_indic.size(); ++j) {
      if (match_indic[j] == target) {
        i_b[i_b_idx++] = j + 1;  // R is 1-based
        found++;
        if (found == 2) break;
      }
    }
  }
}
