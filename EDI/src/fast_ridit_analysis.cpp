#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <map>

using namespace Rcpp;

// [[Rcpp::export]]
List fast_ridit_scores_cpp(const IntegerVector& y, const IntegerVector& ref_idx) {
    int n = y.size();
    int n_ref = ref_idx.size();
    
    // 1. Get unique levels and their counts in the reference group
    std::map<int, int> counts;
    for (int i = 0; i < n_ref; ++i) {
        counts[y[ref_idx[i] - 1]]++;
    }
    
    std::vector<int> levels;
    for (auto const& [level, count] : counts) {
        levels.push_back(level);
    }
    std::sort(levels.begin(), levels.end());
    
    int K = levels.size();
    std::vector<double> p(K);
    for (int k = 0; k < K; ++k) {
        p[k] = static_cast<double>(counts[levels[k]]) / n_ref;
    }
    
    // 2. Calculate Ridit scores for each category
    // R_k = sum_{j < k} p_j + 0.5 * p_k
    std::map<int, double> ridit_map;
    double cumulative_p = 0.0;
    for (int k = 0; k < K; ++k) {
        ridit_map[levels[k]] = cumulative_p + 0.5 * p[k];
        cumulative_p += p[k];
    }
    
    // 3. Assign scores to all subjects
    NumericVector scores(n);
    for (int i = 0; i < n; ++i) {
        if (ridit_map.count(y[i])) {
            scores[i] = ridit_map[y[i]];
        } else {
            // If a level wasn't in the reference group, find its place
            auto it = std::lower_bound(levels.begin(), levels.end(), y[i]);
            if (it == levels.begin()) {
                scores[i] = 0.0; // Extremely low
            } else if (it == levels.end()) {
                scores[i] = 1.0; // Extremely high
            } else {
                // Average of the ridits of the categories it falls between
                int idx = std::distance(levels.begin(), it);
                // This is a heuristic for unseen categories
                scores[i] = (ridit_map[levels[idx-1]] + ridit_map[levels[idx]]) / 2.0;
            }
        }
    }
    
    NumericVector ref_p(K);
    for(int k = 0; k < K; ++k) ref_p[k] = p[k];
    
    return List::create(
        Named("scores") = scores,
        Named("levels") = wrap(levels),
        Named("ref_p") = ref_p
    );
}

// [[Rcpp::export]]
List fast_ridit_analysis_cpp(const IntegerVector& y, const IntegerVector& w, const std::string& reference = "control") {
    int n = y.size();
    std::vector<int> ref_idx;
    
    if (reference == "control") {
        for (int i = 0; i < n; ++i) if (w[i] == 0) ref_idx.push_back(i + 1);
    } else if (reference == "treatment") {
        for (int i = 0; i < n; ++i) if (w[i] == 1) ref_idx.push_back(i + 1);
    } else { // pooled
        for (int i = 0; i < n; ++i) ref_idx.push_back(i + 1);
    }
    
    if (ref_idx.empty()) return List::create();
    
    List ridit_data = fast_ridit_scores_cpp(y, wrap(ref_idx));
    NumericVector scores = ridit_data["scores"];
    
    double sum_t = 0.0, sum_c = 0.0;
    int n_t = 0, n_c = 0;
    std::vector<double> scores_t, scores_c;
    
    for (int i = 0; i < n; ++i) {
        if (w[i] == 1) {
            sum_t += scores[i];
            n_t++;
            scores_t.push_back(scores[i]);
        } else {
            sum_c += scores[i];
            n_c++;
            scores_c.push_back(scores[i]);
        }
    }
    
    double mean_ridit_t = (n_t > 0) ? sum_t / n_t : NA_REAL;
    double mean_ridit_c = (n_c > 0) ? sum_c / n_c : NA_REAL;
    
    // Variance of mean ridit (Bross 1958 approximation)
    // Var(W) = 1 / (12 * n_t) if reference is control and n_c is large
    // More generally, we can use the sample variance of the scores
    double var_t = 0.0;
    if (n_t > 1) {
        for (double s : scores_t) var_t += std::pow(s - mean_ridit_t, 2);
        var_t /= (n_t - 1);
    }
    
    return List::create(
        Named("mean_ridit_t") = mean_ridit_t,
        Named("mean_ridit_c") = mean_ridit_c,
        Named("estimate") = mean_ridit_t - 0.5, // Centered at 0
        Named("se") = std::sqrt(var_t / n_t),
        Named("scores") = scores,
        Named("levels") = ridit_data["levels"],
        Named("ref_p") = ridit_data["ref_p"]
    );
}
