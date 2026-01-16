#include <Rcpp.h>

using namespace Rcpp;

namespace {

double sample_variance(const NumericVector& x) {
  int n = x.size();
  if (n <= 1) {
    return NA_REAL;
  }
  double sum = 0.0;
  double sumsq = 0.0;
  for (int i = 0; i < n; ++i) {
    double v = x[i];
    sum += v;
    sumsq += v * v;
  }
  double mean = sum / n;
  double var = (sumsq - n * mean * mean) / (n - 1);
  return var;
}

double mean_or_na(const NumericVector& x) {
  int n = x.size();
  if (n == 0) {
    return NA_REAL;
  }
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    sum += x[i];
  }
  return sum / n;
}

} // namespace

// [[Rcpp::export]]
List compute_kk_reservoir_stats_cpp(
  const NumericVector& y_matched_diffs,
  const NumericVector& y_reservoir,
  const NumericVector& w_reservoir
) {
  int m = y_matched_diffs.size();
  int nR = w_reservoir.size();
  int nRT = 0;
  int nRC = 0;
  NumericVector y_reservoir_T;
  NumericVector y_reservoir_C;
  for (int i = 0; i < nR; ++i) {
    if (w_reservoir[i] == 1) {
      y_reservoir_T.push_back(y_reservoir[i]);
      ++nRT;
    } else {
      y_reservoir_C.push_back(y_reservoir[i]);
      ++nRC;
    }
  }

  double ssqD_bar = NA_REAL;
  if (m > 1) {
    double varD = sample_variance(y_matched_diffs);
    if (!ISNA(varD)) {
      ssqD_bar = varD / m;
    }
  }
  double r_bar = NA_REAL;
  if (nRT > 0 && nRC > 0) {
    double meanT = mean_or_na(y_reservoir_T);
    double meanC = mean_or_na(y_reservoir_C);
    r_bar = meanT - meanC;
  }

  double ssqR = NA_REAL;
  if (nRT > 1 && nRC > 1 && (nRT + nRC) > 2) {
    double varT = sample_variance(y_reservoir_T);
    double varC = sample_variance(y_reservoir_C);
    if (!ISNA(varT) && !ISNA(varC)) {
      ssqR = (varT * (nRT - 1) + varC * (nRC - 1)) / (nRT + nRC - 2) * (1.0 / nRT + 1.0 / nRC);
    }
  }

  double d_bar = NA_REAL;
  if (m > 0) {
    d_bar = mean_or_na(y_matched_diffs);
  }
  double w_star = NA_REAL;
  if (!ISNA(ssqR) && !ISNA(ssqD_bar)) {
    double denom = ssqR + ssqD_bar;
    if (denom != 0) {
      w_star = ssqR / denom;
    }
  }

  return List::create(
    _["d_bar"] = d_bar,
    _["ssqD_bar"] = ssqD_bar,
    _["r_bar"] = r_bar,
    _["ssqR"] = ssqR,
    _["w_star"] = w_star,
    _["nRT"] = nRT,
    _["nRC"] = nRC
  );
}
