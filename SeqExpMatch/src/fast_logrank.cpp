#include "_helper_functions.h"
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

struct SubjectRecord {
  double time;
  int dead;
  int w;
};

inline bool record_less(const SubjectRecord& a, const SubjectRecord& b) {
  if (a.time < b.time) return true;
  if (a.time > b.time) return false;
  if (a.dead > b.dead) return true;
  if (a.dead < b.dead) return false;
  return a.w < b.w;
}

ModelResult fast_logrank_internal(const Eigen::VectorXd& time,
                                const std::vector<int>& dead,
                                const std::vector<int>& w) {
  const int n = time.size();
  ModelResult res;
  if (n == 0) return res;

  std::vector<SubjectRecord> recs;
  recs.reserve(n);
  int n_treat = 0;
  int n_control = 0;
  for (int i = 0; i < n; ++i) {
    recs.push_back(SubjectRecord{time[i], dead[i], w[i]});
    if (w[i] == 1) ++n_treat;
    else ++n_control;
  }
  std::sort(recs.begin(), recs.end(), record_less);

  std::vector<double> martingale(n, NA_REAL);
  double score = 0.0;
  double var_score = 0.0;
  double cum_hazard = 0.0;
  int risk_all = n;
  int risk_treat = n_treat;

  int start = 0;
  while (start < n) {
    const double curr_time = recs[start].time;
    int end = start;
    int d_all = 0;
    int d_treat = 0;
    int remove_treat = 0;

    while (end < n && recs[end].time == curr_time) {
      const SubjectRecord& rec = recs[end];
      if (rec.dead == 1) {
        ++d_all;
        if (rec.w == 1) ++d_treat;
      }
      if (rec.w == 1) ++remove_treat;
      ++end;
    }

    if (d_all > 0 && risk_all > 0) {
      const double expected_treat = static_cast<double>(d_all) * static_cast<double>(risk_treat) / static_cast<double>(risk_all);
      score += static_cast<double>(d_treat) - expected_treat;
      if (risk_all > 1) {
        const double frac_treat = static_cast<double>(risk_treat) / static_cast<double>(risk_all);
        var_score += static_cast<double>(d_all) * frac_treat * (1.0 - frac_treat) * (static_cast<double>(risk_all - d_all) / static_cast<double>(risk_all - 1));
      }
      cum_hazard += static_cast<double>(d_all) / static_cast<double>(risk_all);
    }

    for (int i = start; i < end; ++i) martingale[i] = static_cast<double>(recs[i].dead) - cum_hazard;
    risk_all -= (end - start);
    risk_treat -= remove_treat;
    start = end;
  }

  double mean_treat = 0.0, mean_control = 0.0;
  for (int i = 0; i < n; ++i) {
    if (recs[i].w == 1) mean_treat += martingale[i];
    else mean_control += martingale[i];
  }
  if (n_treat > 0) mean_treat /= static_cast<double>(n_treat);
  if (n_control > 0) mean_control /= static_cast<double>(n_control);
  res.b = Eigen::VectorXd(1);
  res.b[0] = (n_treat > 0 && n_control > 0) ? (mean_treat - mean_control) : NA_REAL;

  double var_treat = 0.0, var_control = 0.0;
  if (n_treat > 1) {
    double ss = 0.0;
    for (int i = 0; i < n; ++i) if (recs[i].w == 1) { double diff = martingale[i] - mean_treat; ss += diff * diff; }
    var_treat = ss / static_cast<double>(n_treat - 1);
  }
  if (n_control > 1) {
    double ss = 0.0;
    for (int i = 0; i < n; ++i) if (recs[i].w == 0) { double diff = martingale[i] - mean_control; ss += diff * diff; }
    var_control = ss / static_cast<double>(n_control - 1);
  }

  if (std::isfinite(var_treat) && std::isfinite(var_control)) {
    const double se_sq = var_treat / static_cast<double>(n_treat) + var_control / static_cast<double>(n_control);
    if (std::isfinite(se_sq) && se_sq > 0.0) res.ssq_b_2 = std::sqrt(se_sq); // repurposing ssq_b_2 for SE
  }
  res.dispersion = score; // repurposing dispersion for score
  res.sigma2_hat = var_score; // repurposing sigma2_hat for var_score
  return res;
}

} // namespace

// [[Rcpp::export]]
List fast_logrank_stats_cpp(const Eigen::VectorXd& time,
                            const IntegerVector& dead,
                            const IntegerVector& w) {
  int n = time.size();
  std::vector<int> dead_std(n), w_std(n);
  for(int i=0; i<n; ++i) { dead_std[i] = dead[i]; w_std[i] = w[i]; }

  ModelResult res = fast_logrank_internal(time, dead_std, w_std);
  int n_treat = 0;
  for(int val : w_std) if (val == 1) n_treat++;

  return List::create(
    _["score"] = res.dispersion,
    _["var_score"] = (res.sigma2_hat > 0.0) ? res.sigma2_hat : NA_REAL,
    _["beta_hat"] = res.b[0],
    _["se_beta_hat"] = res.ssq_b_2,
    _["n_treat"] = n_treat,
    _["n_control"] = n - n_treat
  );
}
