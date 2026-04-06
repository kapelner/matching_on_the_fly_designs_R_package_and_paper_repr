// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <algorithm>
#include <map>
#include <vector>

using namespace Rcpp;

namespace {

struct CandidateColumn {
  int index;
  int num_levels;
};

inline bool candidate_less(const CandidateColumn& a, const CandidateColumn& b) {
  if (a.num_levels < b.num_levels) return true;
  if (a.num_levels > b.num_levels) return false;
  return a.index < b.index;
}

}

//' Compute automatic survival strata IDs from low-cardinality covariates
//'
//' Selects numeric covariate columns with a small number of observed levels and
//' combines them into a single all-subject stratum identifier. This is used to
//' support automatic stratified Cox models when the design object stores observed
//' covariates but no explicit stratum variable.
//'
//' @param X Numeric covariate matrix.
//' @param max_unique_per_col Maximum number of unique values allowed for a column
//'   to be considered a stratification candidate.
//' @param max_strata_cols Maximum number of candidate columns to combine.
//' @param min_count_per_level Minimum frequency required for every level in a
//'   candidate column.
//' @return A list with `strata_id`, `selected_cols`, and `num_strata`.
// [[Rcpp::export]]
List compute_survival_strata_ids_cpp(const Eigen::MatrixXd& X,
                                     int max_unique_per_col = 4,
                                     int max_strata_cols = 4,
                                     int min_count_per_level = 2) {
  const int n = X.rows();
  const int p = X.cols();
  IntegerVector default_ids(n, 1);

  if (n == 0 || p == 0 || max_unique_per_col < 2 || max_strata_cols < 1) {
    return List::create(
      _["strata_id"] = default_ids,
      _["selected_cols"] = IntegerVector(0),
      _["num_strata"] = 1
    );
  }

  std::vector< std::map<double, int> > level_maps(static_cast<std::size_t>(p));
  std::vector<CandidateColumn> candidates;
  candidates.reserve(static_cast<std::size_t>(p));

  for (int j = 0; j < p; ++j) {
    std::map<double, int> counts;
    bool ok = true;
    for (int i = 0; i < n; ++i) {
      const double xij = X(i, j);
      if (!R_finite(xij)) {
        ok = false;
        break;
      }
      counts[xij] += 1;
      if (static_cast<int>(counts.size()) > max_unique_per_col) {
        ok = false;
        break;
      }
    }
    if (!ok) continue;
    if (static_cast<int>(counts.size()) < 2) continue;

    for (std::map<double, int>::const_iterator it = counts.begin(); it != counts.end(); ++it) {
      if (it->second < min_count_per_level) {
        ok = false;
        break;
      }
    }
    if (!ok) continue;

    level_maps[static_cast<std::size_t>(j)] = counts;
    candidates.push_back(CandidateColumn{j, static_cast<int>(counts.size())});
  }

  if (candidates.empty()) {
    return List::create(
      _["strata_id"] = default_ids,
      _["selected_cols"] = IntegerVector(0),
      _["num_strata"] = 1
    );
  }

  std::sort(candidates.begin(), candidates.end(), candidate_less);
  if (static_cast<int>(candidates.size()) > max_strata_cols) {
    candidates.resize(static_cast<std::size_t>(max_strata_cols));
  }

  IntegerVector selected_cols(static_cast<int>(candidates.size()));
  std::vector< std::map<double, int> > selected_level_codes(candidates.size());
  std::vector<int> multipliers(candidates.size(), 1);

  for (std::size_t k = 0; k < candidates.size(); ++k) {
    const int j = candidates[k].index;
    selected_cols[static_cast<int>(k)] = j + 1;

    int code = 1;
    for (std::map<double, int>::const_iterator it = level_maps[static_cast<std::size_t>(j)].begin();
         it != level_maps[static_cast<std::size_t>(j)].end(); ++it) {
      selected_level_codes[k][it->first] = code;
      ++code;
    }
    if (k > 0) {
      multipliers[k] = multipliers[k - 1] * (candidates[k - 1].num_levels + 1);
    }
  }

  IntegerVector raw_ids(n);
  std::map<int, int> id_lookup;
  int next_id = 1;
  for (int i = 0; i < n; ++i) {
    int combined = 0;
    for (std::size_t k = 0; k < candidates.size(); ++k) {
      const double xij = X(i, candidates[k].index);
      const int code = selected_level_codes[k][xij];
      combined += code * multipliers[k];
    }
    std::map<int, int>::iterator it = id_lookup.find(combined);
    if (it == id_lookup.end()) {
      id_lookup[combined] = next_id;
      raw_ids[i] = next_id;
      ++next_id;
    } else {
      raw_ids[i] = it->second;
    }
  }

  return List::create(
    _["strata_id"] = raw_ids,
    _["selected_cols"] = selected_cols,
    _["num_strata"] = next_id - 1
  );
}
