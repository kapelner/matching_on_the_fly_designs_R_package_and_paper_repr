// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <algorithm>
#include <map>
#include <vector>

using namespace Rcpp;

namespace {

double log_choose_int(int n, int k) {
  if (k < 0 || k > n) return R_NegInf;
  return R::lchoose(static_cast<double>(n), static_cast<double>(k));
}

int compute_jt_statistic2_from_counts(const std::vector<int>& treat_counts,
                                      const std::vector<int>& total_counts) {
  const int K = static_cast<int>(total_counts.size());
  int treated_seen = 0;
  int total_seen = 0;
  int stat2 = 0;

  for (int k = 0; k < K; ++k) {
    const int tk = treat_counts[static_cast<std::size_t>(k)];
    const int nk = total_counts[static_cast<std::size_t>(k)];
    const int lower_controls = total_seen - treated_seen;
    stat2 += tk * (2 * lower_controls + (nk - tk));
    treated_seen += tk;
    total_seen += nk;
  }
  return stat2;
}

void recurse_jt_distribution(int idx,
                             int remaining_treated,
                             const std::vector<int>& total_counts,
                             std::vector<int>& treat_counts,
                             double log_weight,
                             std::map<int, double>& stat_prob,
                             double log_norm) {
  const int K = static_cast<int>(total_counts.size());
  if (idx == K - 1) {
    const int nk = total_counts[static_cast<std::size_t>(idx)];
    if (remaining_treated < 0 || remaining_treated > nk) return;
    treat_counts[static_cast<std::size_t>(idx)] = remaining_treated;
    const double lw = log_weight + log_choose_int(nk, remaining_treated) - log_norm;
    const int stat2 = compute_jt_statistic2_from_counts(treat_counts, total_counts);
    stat_prob[stat2] += std::exp(lw);
    return;
  }

  const int nk = total_counts[static_cast<std::size_t>(idx)];
  const int max_take = std::min(nk, remaining_treated);
  for (int tk = 0; tk <= max_take; ++tk) {
    treat_counts[static_cast<std::size_t>(idx)] = tk;
    recurse_jt_distribution(
      idx + 1,
      remaining_treated - tk,
      total_counts,
      treat_counts,
      log_weight + log_choose_int(nk, tk),
      stat_prob,
      log_norm
    );
  }
}

}  // namespace

// [[Rcpp::export]]
List exact_jonckheere_terpstra_pval_cpp(const IntegerVector& y,
                                        const IntegerVector& w) {
  const int n = y.size();
  if (w.size() != n) stop("dimension mismatch in exact_jonckheere_terpstra_pval_cpp");
  if (n == 0) stop("empty input in exact_jonckheere_terpstra_pval_cpp");

  const int* y_ptr = y.begin();
  const int* w_ptr = w.begin();

  int n_treat = 0;
  int n_control = 0;
  for (int i = 0; i < n; ++i) {
    if (IntegerVector::is_na(y_ptr[i]) || IntegerVector::is_na(w_ptr[i])) {
      stop("missing values are not allowed in exact_jonckheere_terpstra_pval_cpp");
    }
    if (w_ptr[i] == 1) {
      ++n_treat;
    } else if (w_ptr[i] == 0) {
      ++n_control;
    } else {
      stop("treatment assignments must be 0/1 in exact_jonckheere_terpstra_pval_cpp");
    }
  }
  if (n_treat == 0 || n_control == 0) {
    stop("both treatment arms must be present in exact_jonckheere_terpstra_pval_cpp");
  }

  std::vector<int> levels(n);
  for(int i=0; i<n; ++i) levels[i] = y_ptr[i];
  std::sort(levels.begin(), levels.end());
  levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
  const int K = static_cast<int>(levels.size());

  std::vector<int> total_counts(static_cast<std::size_t>(K), 0);
  std::vector<int> treat_counts_obs(static_cast<std::size_t>(K), 0);

  for (int i = 0; i < n; ++i) {
    const int yi = y_ptr[i];
    const int k = static_cast<int>(std::lower_bound(levels.begin(), levels.end(), yi) - levels.begin());
    ++total_counts[static_cast<std::size_t>(k)];
    if (w_ptr[i] == 1) ++treat_counts_obs[static_cast<std::size_t>(k)];
  }

  const int stat2_obs = compute_jt_statistic2_from_counts(treat_counts_obs, total_counts);
  const double log_norm = log_choose_int(n, n_treat);
  std::vector<int> treat_counts(static_cast<std::size_t>(K), 0);
  std::map<int, double> stat_prob;
  recurse_jt_distribution(0, n_treat, total_counts, treat_counts, 0.0, stat_prob, log_norm);

  double p_lower = 0.0;
  double p_upper = 0.0;
  for (auto const& [stat2, prob] : stat_prob) {
    if (stat2 <= stat2_obs) p_lower += prob;
    if (stat2 >= stat2_obs) p_upper += prob;
  }

  const double p_exact = std::min(1.0, 2.0 * std::min(p_lower, p_upper));
  const double superiority = static_cast<double>(stat2_obs) / (2.0 * n_treat * n_control);

  return List::create(
    _["stat2"] = stat2_obs,
    _["n_treat"] = n_treat,
    _["n_control"] = n_control,
    _["superiority"] = superiority,
    _["p_lower"] = p_lower,
    _["p_upper"] = p_upper,
    _["p_exact"] = p_exact
  );
}
