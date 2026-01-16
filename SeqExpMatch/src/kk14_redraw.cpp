#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector redraw_w_kk14_cpp(const IntegerVector& match_indic,
                                const NumericVector& w) {
  int n = w.size();
  int m = 0;
  for (int i = 0; i < n; ++i) {
    int id = match_indic[i];
    if (id > m) {
      m = id;
    }
  }

  std::vector<int> idx1(m, -1);
  std::vector<int> idx2(m, -1);
  std::vector<int> reservoir;
  reservoir.reserve(n);

  for (int i = 0; i < n; ++i) {
    int id = match_indic[i];
    if (id > 0 && id <= m) {
      int idx = id - 1;
      if (idx1[idx] == -1) {
        idx1[idx] = i;
      } else {
        idx2[idx] = i;
      }
    } else {
      reservoir.push_back(i);
    }
  }

  RNGScope scope;
  NumericVector w_out = clone(w);

  for (int i = 0; i < m; ++i) {
    if (idx1[i] >= 0 && idx2[i] >= 0) {
      double u = R::unif_rand();
      if (u < 0.5) {
        w_out[idx1[i]] = 0.0;
        w_out[idx2[i]] = 1.0;
      } else {
        w_out[idx1[i]] = 1.0;
        w_out[idx2[i]] = 0.0;
      }
    }
  }

  int rsize = static_cast<int>(reservoir.size());
  if (rsize > 1) {
    NumericVector vals(rsize);
    for (int i = 0; i < rsize; ++i) {
      vals[i] = w_out[reservoir[i]];
    }
    for (int i = rsize - 1; i > 0; --i) {
      int j = static_cast<int>(std::floor(R::unif_rand() * (i + 1)));
      double tmp = vals[i];
      vals[i] = vals[j];
      vals[j] = tmp;
    }
    for (int i = 0; i < rsize; ++i) {
      w_out[reservoir[i]] = vals[i];
    }
  }

  return w_out;
}
