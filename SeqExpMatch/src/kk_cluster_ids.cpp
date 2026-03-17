#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List compute_kk_grouping_cpp(const IntegerVector& match_indic) {
  const int n = match_indic.size();
  IntegerVector cluster_id(n);
  NumericVector pair_active(n);
  NumericVector reservoir_ind(n);

  int max_match = 0;
  for (int i = 0; i < n; ++i) {
    int match_id = match_indic[i];
    if (match_id == NA_INTEGER) {
      match_id = 0;
    }
    if (match_id > max_match) {
      max_match = match_id;
    }
  }

  int next_id = max_match;
  for (int i = 0; i < n; ++i) {
    int match_id = match_indic[i];
    if (match_id == NA_INTEGER) {
      match_id = 0;
    }
    if (match_id > 0) {
      cluster_id[i] = match_id;
      pair_active[i] = 1.0;
      reservoir_ind[i] = 0.0;
    } else {
      ++next_id;
      cluster_id[i] = next_id;
      pair_active[i] = 0.0;
      reservoir_ind[i] = 1.0;
    }
  }

  return List::create(
    _["cluster_id"] = cluster_id,
    _["pair_active"] = pair_active,
    _["reservoir_ind"] = reservoir_ind
  );
}

// [[Rcpp::export]]
IntegerVector compute_kk_cluster_ids_cpp(const IntegerVector& match_indic) {
  List grouping = compute_kk_grouping_cpp(match_indic);
  return grouping["cluster_id"];
}
