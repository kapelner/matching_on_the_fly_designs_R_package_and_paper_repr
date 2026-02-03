#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector redraw_w_kk14_cpp(const IntegerVector& match_indic, const NumericVector& w) {
  (void)w; // Unused; kept for ABI consistency with generated RcppExports.
  int n = match_indic.size();
  NumericVector w_out(n);
  double sum_w;
  double sum_w_reservoir;

  RNGScope scope;

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
  int rsize = static_cast<int>(reservoir.size());

  do {
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

    for (int i = 0; i < rsize; ++i) {
        w_out[reservoir[i]] = (R::unif_rand() < 0.5) ? 0.0 : 1.0;
    }

    sum_w = 0;
    for (int i = 0; i < n; ++i) {
        sum_w += w_out[i];
    }
    
    sum_w_reservoir = 0;
    for (int i = 0; i < rsize; ++i) {
        sum_w_reservoir += w_out[reservoir[i]];
    }
    
  } while (sum_w == 0 || sum_w == n || (rsize > 1 && (sum_w_reservoir == 0 || sum_w_reservoir == rsize)));

  return w_out;
}
