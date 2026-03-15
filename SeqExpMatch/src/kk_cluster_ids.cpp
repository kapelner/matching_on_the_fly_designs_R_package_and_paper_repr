#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector compute_kk_cluster_ids_cpp(const IntegerVector& match_indic) {
  const int n = match_indic.size();
  IntegerVector cluster_id(n);

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
    } else {
      ++next_id;
      cluster_id[i] = next_id;
    }
  }

  return cluster_id;
}
