#include <Rcpp.h>
#include <limits>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_pair_averages_cpp(const NumericMatrix& X,
                                        const IntegerVector& match_indic,
                                        const int m) {
  int n = X.nrow();
  int p = X.ncol();

  NumericMatrix pair_avg(m, p);
  IntegerVector counts(m);

  for (int i = 0; i < n; ++i) {
    int id = match_indic[i];
    if (id > 0 && id <= m) {
      int idx = id - 1;
      counts[idx] += 1;
      for (int j = 0; j < p; ++j) {
        pair_avg(idx, j) += X(i, j);
      }
    }
  }

  for (int i = 0; i < m; ++i) {
    if (counts[i] > 0) {
      double denom = static_cast<double>(counts[i]);
      for (int j = 0; j < p; ++j) {
        pair_avg(i, j) /= denom;
      }
    } else {
      for (int j = 0; j < p; ++j) {
        pair_avg(i, j) = NA_REAL;
      }
    }
  }

  return pair_avg;
}

// [[Rcpp::export]]
NumericMatrix compute_pair_distance_matrix_cpp(const NumericMatrix& pair_avg,
                                               const NumericVector& weights) {
  int m = pair_avg.nrow();
  int p = pair_avg.ncol();
  bool use_weights = weights.size() == p;

  NumericMatrix dist_mat(m, m);
  double inf = std::numeric_limits<double>::infinity();

  for (int i = 0; i < m; ++i) {
    dist_mat(i, i) = inf;
  }

  for (int i = 0; i < m; ++i) {
    for (int j = i + 1; j < m; ++j) {
      double sum = 0.0;
      for (int k = 0; k < p; ++k) {
        double diff = pair_avg(i, k) - pair_avg(j, k);
        double w = use_weights ? weights[k] : 1.0;
        sum += w * diff * diff;
      }
      dist_mat(i, j) = sum;
      dist_mat(j, i) = sum;
    }
  }

  return dist_mat;
}

// [[Rcpp::export]]
double compute_lambda_squ_cpp(const NumericVector& d_i,
                              const IntegerMatrix& halves) {
  int n_halves = halves.nrow();
  if (n_halves == 0) {
    return 0.0;
  }

  double sum = 0.0;
  for (int i = 0; i < n_halves; ++i) {
    int id1 = halves(i, 0) - 1;
    int id2 = halves(i, 1) - 1;
    if (id1 >= 0 && id2 >= 0 && id1 < d_i.size() && id2 < d_i.size()) {
      sum += d_i[id1] * d_i[id2];
    }
  }

  return sum / static_cast<double>(n_halves);
}
