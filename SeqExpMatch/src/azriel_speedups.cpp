#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

namespace {

struct BlockAccumulator {
  int n = 0;
  double sum = 0.0;
  double sum_sq = 0.0;
};

}  // namespace

// [[Rcpp::export]]
double compute_azriel_block_se_cpp(const NumericVector& y,
                                   const IntegerVector& m_vec,
                                   int n_total) {
  if (y.size() != m_vec.size()) {
    stop("compute_azriel_block_se_cpp: y and m_vec must have the same length.");
  }
  if (n_total <= 0) {
    return NA_REAL;
  }

  std::unordered_map<int, BlockAccumulator> blocks;
  blocks.reserve(static_cast<std::size_t>(y.size()));

  for (int i = 0; i < y.size(); ++i) {
    const int match_id = m_vec[i];
    if (match_id == NA_INTEGER || match_id <= 0 || !R_finite(y[i])) {
      continue;
    }
    BlockAccumulator& block = blocks[match_id];
    block.n += 1;
    block.sum += y[i];
    block.sum_sq += y[i] * y[i];
  }

  double var_cmh = 0.0;
  for (const auto& entry : blocks) {
    const BlockAccumulator& block = entry.second;
    if (block.n <= 1) {
      return NA_REAL;
    }
    const double ss = block.sum_sq - (block.sum * block.sum) / static_cast<double>(block.n);
    var_cmh += (static_cast<double>(block.n) / (block.n - 1.0)) * ss;
  }

  return (2.0 / static_cast<double>(n_total)) * std::sqrt(var_cmh);
}

// [[Rcpp::export]]
double compute_extended_robins_block_se_cpp(const NumericVector& y,
                                            const NumericVector& w,
                                            const IntegerVector& m_vec,
                                            int n_total) {
  if (y.size() != m_vec.size() || y.size() != w.size()) {
    stop("compute_extended_robins_block_se_cpp: y, w, and m_vec must have the same length.");
  }
  if (n_total <= 0) {
    return NA_REAL;
  }

  struct BlockAccumulator {
    int n = 0;
    double sum_t = 0.0;
    double sum_c = 0.0;
  };

  std::unordered_map<int, BlockAccumulator> blocks;
  blocks.reserve(static_cast<std::size_t>(y.size()));

  for (int i = 0; i < y.size(); ++i) {
    const int match_id = m_vec[i];
    if (match_id == NA_INTEGER || match_id <= 0) {
      continue;
    }
    if (!R_finite(y[i]) || !R_finite(w[i]) || (w[i] != 0.0 && w[i] != 1.0)) {
      return NA_REAL;
    }

    BlockAccumulator& block = blocks[match_id];
    block.n += 1;
    if (w[i] == 1.0) {
      block.sum_t += y[i];
    } else {
      block.sum_c += y[i];
    }
  }

  if (blocks.empty()) {
    return NA_REAL;
  }

  double variance_tot = 0.0;
  for (const auto& entry : blocks) {
    const BlockAccumulator& block = entry.second;
    if (block.n <= 1) {
      return NA_REAL;
    }

    const double n_b = static_cast<double>(block.n);
    const double n_b_over_two = n_b / 2.0;
    const double p_hat_T_b = block.sum_t / n_b_over_two;
    const double p_hat_C_b = block.sum_c / n_b_over_two;
    const double m_1_b = std::max(p_hat_T_b, p_hat_C_b);
    const double m_0_b = std::min(p_hat_T_b, p_hat_C_b);

    variance_tot +=
      m_1_b * (1.0 - m_1_b) / n_b_over_two +
      m_0_b * (1.0 - m_0_b) / n_b_over_two +
      ((2.0 * m_0_b - m_1_b) * (1.0 - m_1_b) - m_0_b * (1.0 - m_0_b)) / n_b;
  }

  return (1.0 / static_cast<double>(blocks.size())) * std::sqrt(variance_tot);
}
