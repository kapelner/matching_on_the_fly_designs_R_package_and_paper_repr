#include <RcppEigen.h>
#include <vector>
#include <array>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix bootstrap_match_indices_cpp(
  const IntegerVector& match_indic,
  const IntegerVector& i_reservoir,
  int n_reservoir,
  int m,
  int B
) {
  int row_length = n_reservoir + 2 * m;
  IntegerMatrix result(B, row_length);

  std::vector< std::array<int, 2> > match_pairs(m);
  std::vector<int> count(m, 0);
  for (int idx = 0; idx < match_indic.size(); ++idx) {
    int match_id = match_indic[idx];
    if (match_id > 0 && match_id <= m) {
      int pos = count[match_id - 1]++;
      if (pos < 2) {
        match_pairs[match_id - 1][pos] = idx + 1;
      }
    }
  }

  for (int row = 0; row < B; ++row) {
    if (n_reservoir > 0) {
      for (int j = 0; j < n_reservoir; ++j) {
        int idx = static_cast<int>(R::runif(0.0, 1.0) * n_reservoir);
        if (idx == n_reservoir) idx = n_reservoir - 1;
        result(row, j) = i_reservoir[idx];
      }
    }

    for (int k = 0; k < m; ++k) {
      int match_id = static_cast<int>(R::runif(0.0, 1.0) * m);
      if (match_id == m) match_id = m - 1;
      auto pair = match_pairs[match_id];
      result(row, n_reservoir + 2 * k) = pair[0];
      result(row, n_reservoir + 2 * k + 1) = pair[1];
    }
  }
  return result;
}

// forward declaration from match data helper
List match_diffs_cpp(const Eigen::VectorXi& w, const Eigen::VectorXi& match_indic, const Eigen::VectorXd& y, const Eigen::MatrixXd& X, int m);

// [[Rcpp::export]]
List match_stats_from_indices_cpp(
  const NumericVector& y,
  const NumericVector& w,
  const NumericMatrix& X,
  const IntegerVector& original_match_indic, // Changed name
  const IntegerVector& i_b,
  int m
) {
  int n_rows = i_b.size();
  int p = X.ncol();
  NumericVector y_sample(n_rows);
  NumericVector w_sample(n_rows);
  NumericMatrix X_sample(n_rows, p);
  IntegerVector match_indic_sample(n_rows); // This will hold the resampled match_indic

  for (int i = 0; i < n_rows; ++i) {
    int idx = i_b[i] - 1;
    y_sample[i] = y[idx];
    w_sample[i] = w[idx];
    match_indic_sample[i] = original_match_indic[idx]; // Corrected resampling
    for (int j = 0; j < p; ++j) {
      X_sample(i, j) = X(idx, j);
    }
  }

  Eigen::VectorXi w_eigen = as<Eigen::VectorXi>(w_sample);
  Eigen::VectorXi match_indic_eigen = as<Eigen::VectorXi>(match_indic_sample);
  Eigen::VectorXd y_eigen = as<Eigen::VectorXd>(y_sample);
  Eigen::MatrixXd X_eigen = as<Eigen::MatrixXd>(X_sample);
  List match_data = match_diffs_cpp(w_eigen, match_indic_eigen, y_eigen, X_eigen, m);

  std::vector<int> reservoir_rows;
  for (int i = 0; i < n_rows; ++i) {
    if (match_indic_sample[i] == 0) {
      reservoir_rows.push_back(i);
    }
  }

  int n_res = reservoir_rows.size();
  NumericMatrix X_reservoir(n_res, p);
  NumericVector y_reservoir(n_res);
  NumericVector w_reservoir(n_res);
  for (int k = 0; k < n_res; ++k) {
    int row = reservoir_rows[k];
    y_reservoir[k] = y_sample[row];
    w_reservoir[k] = w_sample[row];
    for (int j = 0; j < p; ++j) {
      X_reservoir(k, j) = X_sample(row, j);
    }
  }

  int nRT = 0, nRC = 0;
  for (int i = 0; i < n_res; ++i) {
    if (w_reservoir[i] == 1) nRT++;
    else nRC++;
  }

  NumericMatrix X_matched_diffs = match_data["X_matched_diffs"];
  NumericVector yTs_matched = match_data["yTs_matched"];
  NumericVector yCs_matched = match_data["yCs_matched"];

  return List::create(
    _["X_matched_diffs"] = X_matched_diffs,
    _["yTs_matched"] = yTs_matched,
    _["yCs_matched"] = yCs_matched,
    _["y_matched_diffs"] = yTs_matched - yCs_matched,
    _["X_reservoir"] = X_reservoir,
    _["y_reservoir"] = y_reservoir,
    _["w_reservoir"] = w_reservoir,
    _["nRT"] = nRT,
    _["nRC"] = nRC,
    _["m"] = m
  );
}
