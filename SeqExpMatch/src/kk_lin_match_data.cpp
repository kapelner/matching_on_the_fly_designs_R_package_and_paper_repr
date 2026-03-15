#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List compute_kk_lin_match_data_cpp(const IntegerVector& w,
                                   const IntegerVector& match_indic,
                                   const NumericVector& y,
                                   const NumericMatrix& X) {
  const int n = w.size();
  const int p = X.ncol();
  int m = 0;
  int n_reservoir = 0;

  for (int i = 0; i < n; ++i) {
    int match_id = match_indic[i];
    if (match_id == NA_INTEGER) {
      match_id = 0;
    }
    if (match_id <= 0) {
      ++n_reservoir;
    } else if (match_id > m) {
      m = match_id;
    }
  }

  NumericVector yTs_matched(m, NA_REAL);
  NumericVector yCs_matched(m, NA_REAL);
  NumericVector y_matched_diffs(m, NA_REAL);
  NumericMatrix X_matched_diffs_full(m, p);
  NumericMatrix X_matched_means_full(m, p);
  std::vector<int> found_t(static_cast<std::size_t>(m), 0);
  std::vector<int> found_c(static_cast<std::size_t>(m), 0);

  for (int i = 0; i < n; ++i) {
    int match_id = match_indic[i];
    if (match_id == NA_INTEGER || match_id <= 0) {
      continue;
    }

    const int pair_index = match_id - 1;
    if (w[i] == 1) {
      yTs_matched[pair_index] = y[i];
      found_t[static_cast<std::size_t>(pair_index)] = 1;
      for (int j = 0; j < p; ++j) {
        X_matched_diffs_full(pair_index, j) += X(i, j);
        X_matched_means_full(pair_index, j) += 0.5 * X(i, j);
      }
    } else {
      yCs_matched[pair_index] = y[i];
      found_c[static_cast<std::size_t>(pair_index)] = 1;
      for (int j = 0; j < p; ++j) {
        X_matched_diffs_full(pair_index, j) -= X(i, j);
        X_matched_means_full(pair_index, j) += 0.5 * X(i, j);
      }
    }
  }

  for (int pair_index = 0; pair_index < m; ++pair_index) {
    if (found_t[static_cast<std::size_t>(pair_index)] && found_c[static_cast<std::size_t>(pair_index)]) {
      y_matched_diffs[pair_index] = yTs_matched[pair_index] - yCs_matched[pair_index];
    }
  }

  NumericMatrix X_reservoir(n_reservoir, p);
  NumericVector y_reservoir(n_reservoir);
  IntegerVector w_reservoir(n_reservoir);
  int nRT = 0;
  int nRC = 0;

  for (int i = 0, reservoir_index = 0; i < n; ++i) {
    int match_id = match_indic[i];
    if (match_id == NA_INTEGER) {
      match_id = 0;
    }
    if (match_id > 0) {
      continue;
    }

    y_reservoir[reservoir_index] = y[i];
    w_reservoir[reservoir_index] = w[i];
    for (int j = 0; j < p; ++j) {
      X_reservoir(reservoir_index, j) = X(i, j);
    }
    if (w[i] == 1) {
      ++nRT;
    } else {
      ++nRC;
    }
    ++reservoir_index;
  }

  return List::create(
    _["X_matched_diffs_full"] = X_matched_diffs_full,
    _["X_matched_means_full"] = X_matched_means_full,
    _["y_matched_diffs"] = y_matched_diffs,
    _["yTs_matched"] = yTs_matched,
    _["yCs_matched"] = yCs_matched,
    _["X_reservoir"] = X_reservoir,
    _["y_reservoir"] = y_reservoir,
    _["w_reservoir"] = w_reservoir,
    _["nRT"] = nRT,
    _["nRC"] = nRC,
    _["m"] = m
  );
}
