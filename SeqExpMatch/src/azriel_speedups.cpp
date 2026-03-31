#include <Rcpp.h>
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
