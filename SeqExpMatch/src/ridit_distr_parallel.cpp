#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <map>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// Helper to compute ridit scores and levels from a reference set
void get_ridit_map(const std::vector<int>& y_ref, 
                  std::map<int, double>& ridit_map, 
                  std::vector<int>& levels) {
    int n_ref = y_ref.size();
    if (n_ref == 0) return;

    std::map<int, int> counts;
    for (int val : y_ref) {
        counts[val]++;
    }
    
    for (auto const& [level, count] : counts) {
        levels.push_back(level);
    }
    // std::map keys are already sorted
    
    double cumulative_p = 0.0;
    for (int level : levels) {
        double p_k = static_cast<double>(counts[level]) / n_ref;
        ridit_map[level] = cumulative_p + 0.5 * p_k;
        cumulative_p += p_k;
    }
}

double compute_mean_ridit_with_map(const std::vector<int>& y_target,
                                  const std::map<int, double>& ridit_map,
                                  const std::vector<int>& levels) {
    if (y_target.empty()) return NA_REAL;
    if (ridit_map.empty()) return NA_REAL;

    double sum_t = 0.0;
    for (int val : y_target) {
        if (ridit_map.count(val)) {
            sum_t += ridit_map.at(val);
        } else {
            auto it = std::lower_bound(levels.begin(), levels.end(), val);
            if (it == levels.begin()) {
                // sum_t += 0.0
            } else if (it == levels.end()) {
                sum_t += 1.0;
            } else {
                int idx = std::distance(levels.begin(), it);
                sum_t += (ridit_map.at(levels[idx-1]) + ridit_map.at(levels[idx])) / 2.0;
            }
        }
    }
    return sum_t / y_target.size();
}

// Full analysis for one realization
double compute_single_ridit_estimate(const std::vector<int>& y_b, 
                                   const std::vector<int>& w_b, 
                                   const std::string& reference) {
    int n = y_b.size();
    std::vector<int> y_ref;
    std::vector<int> y_t;
    
    for (int i = 0; i < n; ++i) {
        if (w_b[i] == 1) y_t.push_back(y_b[i]);
        
        if (reference == "control") {
            if (w_b[i] == 0) y_ref.push_back(y_b[i]);
        } else if (reference == "treatment") {
            if (w_b[i] == 1) y_ref.push_back(y_b[i]);
        } else { // pooled
            y_ref.push_back(y_b[i]);
        }
    }
    
    if (y_ref.empty() || y_t.empty()) return NA_REAL;
    
    std::map<int, double> ridit_map;
    std::vector<int> levels;
    get_ridit_map(y_ref, ridit_map, levels);
    
    return compute_mean_ridit_with_map(y_t, ridit_map, levels) - 0.5;
}

// [[Rcpp::export]]
NumericVector compute_ridit_distr_parallel_cpp(const IntegerVector& y, 
                                             const IntegerMatrix& w_mat, 
                                             std::string reference, 
                                             int num_cores) {
    int nsim = w_mat.cols();
    int n = y.size();
    NumericVector results(nsim);
    
    std::vector<int> y_std(n);
    for(int i=0; i<n; ++i) y_std[i] = y[i];

    // Precompute if pooled
    bool is_pooled = (reference == "pooled");
    std::map<int, double> global_ridit_map;
    std::vector<int> global_levels;
    if (is_pooled) {
        get_ridit_map(y_std, global_ridit_map, global_levels);
    }

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < nsim; ++b) {
        std::vector<int> y_t;
        std::vector<int> w_b(n);
        for(int i=0; i<n; ++i) {
            int wb = w_mat(i, b);
            w_b[i] = wb;
            if (wb == 1) y_t.push_back(y_std[i]);
        }
        
        if (is_pooled) {
            if (y_t.empty()) results[b] = NA_REAL;
            else results[b] = compute_mean_ridit_with_map(y_t, global_ridit_map, global_levels) - 0.5;
        } else {
            results[b] = compute_single_ridit_estimate(y_std, w_b, reference);
        }
    }
    
    return results;
}

// [[Rcpp::export]]
NumericVector compute_ridit_bootstrap_parallel_cpp(const IntegerVector& y, 
                                                 const IntegerVector& w, 
                                                 const IntegerMatrix& indices_mat, 
                                                 std::string reference, 
                                                 int num_cores) {
    int B = indices_mat.cols();
    int n = y.size();
    NumericVector results(B);

#ifdef _OPENMP
    if (num_cores > 1) omp_set_num_threads(num_cores);
#endif

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < B; ++b) {
        if (indices_mat(0, b) == -1) {
            results[b] = NA_REAL;
            continue;
        }
        
        std::vector<int> y_b(n);
        std::vector<int> w_b(n);
        for (int i = 0; i < n; ++i) {
            int idx = indices_mat(i, b) - 1; // 1-indexed to 0-indexed
            y_b[i] = y[idx];
            w_b[i] = w[idx];
        }
        results[b] = compute_single_ridit_estimate(y_b, w_b, reference);
    }
    
    return results;
}
