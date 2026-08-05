#ifdef EDI_CORE_ONLY
#include <Eigen/Dense>
#include <limits>
constexpr double NA_REAL = std::numeric_limits<double>::quiet_NaN();
#else
#include <RcppEigen.h>
#endif
#include <algorithm>
#include <cmath>
#include <vector>

#ifndef EDI_CORE_ONLY
using namespace Rcpp;
#endif
using namespace Eigen;

namespace {

struct RiditScoreData {
    std::vector<double> scores;
    std::vector<int> levels;
    std::vector<double> ref_p;
};

RiditScoreData fast_ridit_scores_from_ref_indices(const int* y_ptr,
                                                  int n,
                                                  const std::vector<int>& ref_idx_zero_based) {
    const int n_ref = static_cast<int>(ref_idx_zero_based.size());

    // 1. Get unique levels and their counts in the reference group
    std::vector<int> levels;
    std::vector<int> counts;
    levels.reserve(static_cast<std::size_t>(std::min(n_ref, 16)));
    counts.reserve(levels.capacity());
    for (int idx : ref_idx_zero_based) {
        const int yi = y_ptr[idx];
        auto level_it = std::lower_bound(levels.begin(), levels.end(), yi);
        const std::size_t pos = static_cast<std::size_t>(level_it - levels.begin());
        if (level_it != levels.end() && *level_it == yi) {
            counts[pos]++;
        } else {
            levels.insert(level_it, yi);
            counts.insert(counts.begin() + static_cast<std::ptrdiff_t>(pos), 1);
        }
    }

    const int K = static_cast<int>(levels.size());
    std::vector<double> ref_p(static_cast<std::size_t>(K));
    std::vector<double> ridit_values(static_cast<std::size_t>(K));

    // 2. Calculate Ridit scores for each category
    // R_k = sum_{j < k} p_j + 0.5 * p_k
    double cumulative_p = 0.0;
    for (int k = 0; k < K; ++k) {
        const double pk = static_cast<double>(counts[static_cast<std::size_t>(k)]) / n_ref;
        ref_p[static_cast<std::size_t>(k)] = pk;
        const double ridit = cumulative_p + 0.5 * pk;
        ridit_values[static_cast<std::size_t>(k)] = ridit;
        cumulative_p += pk;
    }

    // 3. Assign scores to all subjects
    std::vector<double> scores(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const int yi = y_ptr[i];
        auto level_it = std::lower_bound(levels.begin(), levels.end(), yi);
        if (level_it != levels.end() && *level_it == yi) {
            scores[static_cast<std::size_t>(i)] = ridit_values[static_cast<std::size_t>(level_it - levels.begin())];
        } else {
            // If a level wasn't in the reference group, find its place
            if (level_it == levels.begin()) {
                scores[static_cast<std::size_t>(i)] = 0.0; // Extremely low
            } else if (level_it == levels.end()) {
                scores[static_cast<std::size_t>(i)] = 1.0; // Extremely high
            } else {
                // Average of the ridits of the categories it falls between
                const int idx = static_cast<int>(std::distance(levels.begin(), level_it));
                scores[static_cast<std::size_t>(i)] = 0.5 * (
                    ridit_values[static_cast<std::size_t>(idx - 1)] +
                    ridit_values[static_cast<std::size_t>(idx)]
                );
            }
        }
    }

    return RiditScoreData{scores, levels, ref_p};
}

}  // namespace

struct RiditAnalysisResult {
    double mean_ridit_t;
    double mean_ridit_c;
    double estimate;
    double se;
    std::vector<double> scores;
    std::vector<int> levels;
    std::vector<double> ref_p;
};

// Portable (EDI_CORE_ONLY-safe) sibling of fast_ridit_analysis_cpp below:
// identical logic, on plain Eigen::Ref inputs / a plain struct return
// instead of SEXP/Rcpp::List, so a separate Python binding translation
// unit can call it. fast_ridit_scores_from_ref_indices above already used
// only std::vector/plain-int-pointer types internally -- the Rcpp coupling
// was solely in RiditScoreData's NumericVector fields (now std::vector<double>)
// and this function's own SEXP boundary. An empty scores vector on return
// is the sentinel for "no reference group" (mirrors fast_ridit_analysis_cpp's
// List::create() early return for that case).
RiditAnalysisResult fast_ridit_analysis_result(const Eigen::Ref<const Eigen::VectorXi>& w,
                                               const Eigen::Ref<const Eigen::VectorXi>& y,
                                               const std::string& reference) {
    int n = static_cast<int>(y.size());
    std::vector<int> ref_idx;
    ref_idx.reserve(static_cast<std::size_t>(n));

    if (reference == "control") {
        for (int i = 0; i < n; ++i) if (w[i] == 0) ref_idx.push_back(i);
    } else if (reference == "treatment") {
        for (int i = 0; i < n; ++i) if (w[i] == 1) ref_idx.push_back(i);
    } else { // pooled
        for (int i = 0; i < n; ++i) ref_idx.push_back(i);
    }

    if (ref_idx.empty()) {
        return RiditAnalysisResult{NA_REAL, NA_REAL, NA_REAL, NA_REAL, {}, {}, {}};
    }

    RiditScoreData ridit_data = fast_ridit_scores_from_ref_indices(y.data(), n, ref_idx);
    const std::vector<double>& scores = ridit_data.scores;

    double sum_t = 0.0, sum_c = 0.0;
    int n_t = 0, n_c = 0;

    for (int i = 0; i < n; ++i) {
        if (w[i] == 1) {
            sum_t += scores[static_cast<std::size_t>(i)];
            n_t++;
        } else {
            sum_c += scores[static_cast<std::size_t>(i)];
            n_c++;
        }
    }

    double mean_ridit_t = (n_t > 0) ? sum_t / n_t : NA_REAL;
    double mean_ridit_c = (n_c > 0) ? sum_c / n_c : NA_REAL;

    double var_t = 0.0;
    if (n_t > 1) {
        for (int i = 0; i < n; ++i) {
            if (w[i] == 1) {
                const double diff = scores[static_cast<std::size_t>(i)] - mean_ridit_t;
                var_t += diff * diff;
            }
        }
        var_t /= (n_t - 1);
    }

    return RiditAnalysisResult{
        mean_ridit_t, mean_ridit_c, mean_ridit_t - 0.5, std::sqrt(var_t / n_t),
        ridit_data.scores, ridit_data.levels, ridit_data.ref_p
    };
}

#ifndef EDI_CORE_ONLY
// [[Rcpp::export]]
List fast_ridit_scores_cpp(SEXP y_sexp, SEXP ref_idx_sexp) {
    IntegerVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXi> y(y_vec.begin(), y_vec.size());
    IntegerVector ref_idx_vec(ref_idx_sexp);
    Eigen::Map<const Eigen::VectorXi> ref_idx(ref_idx_vec.begin(), ref_idx_vec.size());
    int n = y.size();
    int n_ref = ref_idx.size();

    std::vector<int> ref_idx_zero_based;
    ref_idx_zero_based.reserve(static_cast<std::size_t>(n_ref));
    for (int i = 0; i < n_ref; ++i) {
        ref_idx_zero_based.push_back(ref_idx[i] - 1);
    }

    RiditScoreData ridit_data = fast_ridit_scores_from_ref_indices(y.data(), n, ref_idx_zero_based);
    return List::create(
        Named("scores") = wrap(ridit_data.scores),
        Named("levels") = wrap(ridit_data.levels),
        Named("ref_p") = wrap(ridit_data.ref_p)
    );
}

// [[Rcpp::export]]
List fast_ridit_analysis_cpp(SEXP w_sexp, SEXP y_sexp, const std::string& reference = "control") {
    IntegerVector w_vec(w_sexp);
    Eigen::Map<const Eigen::VectorXi> w(w_vec.begin(), w_vec.size());
    IntegerVector y_vec(y_sexp);
    Eigen::Map<const Eigen::VectorXi> y(y_vec.begin(), y_vec.size());

    RiditAnalysisResult res = fast_ridit_analysis_result(w, y, reference);
    if (res.scores.empty()) {
        return List::create();
    }

    return List::create(
        Named("mean_ridit_t") = res.mean_ridit_t,
        Named("mean_ridit_c") = res.mean_ridit_c,
        Named("estimate") = res.estimate,
        Named("se") = res.se,
        Named("scores") = wrap(res.scores),
        Named("levels") = wrap(res.levels),
        Named("ref_p") = wrap(res.ref_p)
    );
}
#endif // EDI_CORE_ONLY
