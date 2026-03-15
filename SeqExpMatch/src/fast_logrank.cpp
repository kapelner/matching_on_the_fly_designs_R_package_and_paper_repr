// [[Rcpp::depends(RcppEigen)]]
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

}  // namespace

//' Fast log-rank summary statistics in C++
//'
//' Computes the standard two-sample log-rank score statistic, its hypergeometric
//' tie-adjusted variance, and a score-style treatment estimate based on the
//' difference in mean martingale residuals under the null pooled hazard.
//'
//' @param time Observed follow-up times.
//' @param dead Event indicators (`1` = event, `0` = censored).
//' @param w Treatment indicators (`1` = treatment, `0` = control).
//' @return A list with the log-rank score, its variance, a martingale-residual
//'   mean-difference estimate, and its standard error.
// [[Rcpp::export]]
List fast_logrank_stats_cpp(const Eigen::VectorXd& time,
                            const IntegerVector& dead,
                            const IntegerVector& w) {
  const int n = time.size();
  if (dead.size() != n || w.size() != n) {
    stop("dimension mismatch in fast_logrank_stats_cpp");
  }
  if (n == 0) {
    return List::create(
      _["score"] = NA_REAL,
      _["var_score"] = NA_REAL,
      _["beta_hat"] = NA_REAL,
      _["se_beta_hat"] = NA_REAL,
      _["n_treat"] = 0,
      _["n_control"] = 0
    );
  }

  std::vector<SubjectRecord> recs;
  recs.reserve(static_cast<std::size_t>(n));
  int n_treat = 0;
  int n_control = 0;
  for (int i = 0; i < n; ++i) {
    const double ti = time[i];
    if (!R_finite(ti) || ti < 0.0) {
      stop("times must be finite and nonnegative");
    }
    if (dead[i] != 0 && dead[i] != 1) {
      stop("dead must contain only 0/1");
    }
    if (w[i] != 0 && w[i] != 1) {
      stop("w must contain only 0/1");
    }
    recs.push_back(SubjectRecord{ti, dead[i], w[i]});
    if (w[i] == 1) {
      ++n_treat;
    } else {
      ++n_control;
    }
  }
  std::sort(recs.begin(), recs.end(), record_less);

  std::vector<double> martingale(static_cast<std::size_t>(n), NA_REAL);
  double score = 0.0;
  double var_score = 0.0;
  double cum_hazard = 0.0;
  int risk_all = n;
  int risk_treat = n_treat;

  int start = 0;
  while (start < n) {
    const double curr_time = recs[static_cast<std::size_t>(start)].time;
    int end = start;
    int d_all = 0;
    int d_treat = 0;
    int remove_treat = 0;

    while (end < n && recs[static_cast<std::size_t>(end)].time == curr_time) {
      const SubjectRecord& rec = recs[static_cast<std::size_t>(end)];
      if (rec.dead == 1) {
        ++d_all;
        if (rec.w == 1) ++d_treat;
      }
      if (rec.w == 1) ++remove_treat;
      ++end;
    }

    if (d_all > 0 && risk_all > 0) {
      const double expected_treat = static_cast<double>(d_all) *
        static_cast<double>(risk_treat) / static_cast<double>(risk_all);
      score += static_cast<double>(d_treat) - expected_treat;

      if (risk_all > 1) {
        const double frac_treat = static_cast<double>(risk_treat) / static_cast<double>(risk_all);
        var_score += static_cast<double>(d_all) * frac_treat * (1.0 - frac_treat) *
          (static_cast<double>(risk_all - d_all) / static_cast<double>(risk_all - 1));
      }

      cum_hazard += static_cast<double>(d_all) / static_cast<double>(risk_all);
    }

    for (int i = start; i < end; ++i) {
      martingale[static_cast<std::size_t>(i)] =
        static_cast<double>(recs[static_cast<std::size_t>(i)].dead) - cum_hazard;
    }

    risk_all -= (end - start);
    risk_treat -= remove_treat;
    start = end;
  }

  double mean_treat = 0.0;
  double mean_control = 0.0;
  for (int i = 0; i < n; ++i) {
    if (recs[static_cast<std::size_t>(i)].w == 1) {
      mean_treat += martingale[static_cast<std::size_t>(i)];
    } else {
      mean_control += martingale[static_cast<std::size_t>(i)];
    }
  }
  if (n_treat > 0) mean_treat /= static_cast<double>(n_treat);
  if (n_control > 0) mean_control /= static_cast<double>(n_control);
  const double beta_hat = (n_treat > 0 && n_control > 0) ? (mean_treat - mean_control) : NA_REAL;

  double var_treat = NA_REAL;
  double var_control = NA_REAL;
  if (n_treat > 1) {
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
      if (recs[static_cast<std::size_t>(i)].w == 1) {
        const double diff = martingale[static_cast<std::size_t>(i)] - mean_treat;
        ss += diff * diff;
      }
    }
    var_treat = ss / static_cast<double>(n_treat - 1);
  }
  if (n_control > 1) {
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
      if (recs[static_cast<std::size_t>(i)].w == 0) {
        const double diff = martingale[static_cast<std::size_t>(i)] - mean_control;
        ss += diff * diff;
      }
    }
    var_control = ss / static_cast<double>(n_control - 1);
  }

  double se_beta_hat = NA_REAL;
  if (R_finite(var_treat) && R_finite(var_control)) {
    const double se_sq = var_treat / static_cast<double>(n_treat) +
      var_control / static_cast<double>(n_control);
    if (R_finite(se_sq) && se_sq > 0.0) {
      se_beta_hat = std::sqrt(se_sq);
    }
  }

  return List::create(
    _["score"] = score,
    _["var_score"] = (var_score > 0.0) ? var_score : NA_REAL,
    _["beta_hat"] = beta_hat,
    _["se_beta_hat"] = se_beta_hat,
    _["n_treat"] = n_treat,
    _["n_control"] = n_control
  );
}
