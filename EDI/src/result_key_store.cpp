#include <Rcpp.h>
#include <unordered_set>
#include <string>
#include <sstream>
#include <iomanip>

using namespace Rcpp;

// Using a static pointer to persist the store across calls within the same R session.
// Note: In a package, this is safe as long as we don't need multiple independent stores
// or we manage them via Rcpp Modules. For this framework, one global store for the 
// current simulation run is sufficient.
static std::unordered_set<std::string>* global_key_store = nullptr;

// Helper to format double to string consistently
inline std::string format_double(double d) {
    std::ostringstream oss;
    oss << std::setprecision(15) << d; // Match a reasonable precision
    return oss.str();
}

// [[Rcpp::export]]
void clear_result_key_store_cpp() {
    if (global_key_store) {
        delete global_key_store;
        global_key_store = nullptr;
    }
}

// [[Rcpp::export]]
void init_result_key_store_cpp(int expected_size) {
    clear_result_key_store_cpp();
    global_key_store = new std::unordered_set<std::string>();
    if (expected_size > 0) {
        global_key_store->reserve(expected_size);
    }
}

// [[Rcpp::export]]
void add_to_result_key_store_cpp(CharacterVector response_type,
                               CharacterVector cond_exp_func_model,
                               IntegerVector n,
                               IntegerVector p,
                               NumericVector betaT,
                               IntegerVector rep,
                               CharacterVector design,
                               CharacterVector inference,
                               CharacterVector inference_type,
                               int start = 1,
                               int end = -1) {
    if (!global_key_store) return;
    
    int n_rows = response_type.size();
    int i_start = std::max(0, start - 1);
    int i_end = (end < 0) ? n_rows : std::min(n_rows, end);
    
    for (int i = i_start; i < i_end; ++i) {
        std::string key;
        key.reserve(128); // Pre-allocate to reduce reallocations
        key += response_type[i];
        key += "|";
        key += cond_exp_func_model[i];
        key += "|";
        key += std::to_string(n[i]);
        key += "|";
        key += std::to_string(p[i]);
        key += "|";
        key += format_double(betaT[i]);
        key += "|";
        key += std::to_string(rep[i]);
        key += "|";
        key += design[i];
        key += "|";
        key += inference[i];
        key += "|";
        key += inference_type[i];
        global_key_store->insert(key);
    }
}

// [[Rcpp::export]]
LogicalVector check_in_result_key_store_cpp(CharacterVector response_type,
                                          CharacterVector cond_exp_func_model,
                                          IntegerVector n,
                                          IntegerVector p,
                                          NumericVector betaT,
                                          IntegerVector rep,
                                          CharacterVector design,
                                          CharacterVector inference,
                                          CharacterVector inference_type,
                                          int start = 1,
                                          int end = -1) {
    int n_rows = response_type.size();
    int i_start = std::max(0, start - 1);
    int i_end = (end < 0) ? n_rows : std::min(n_rows, end);
    
    LogicalVector results(i_end - i_start);
    if (!global_key_store) return results;
    
    for (int i = i_start; i < i_end; ++i) {
        std::string key;
        key.reserve(128);
        key += response_type[i];
        key += "|";
        key += cond_exp_func_model[i];
        key += "|";
        key += std::to_string(n[i]);
        key += "|";
        key += std::to_string(p[i]);
        key += "|";
        key += format_double(betaT[i]);
        key += "|";
        key += std::to_string(rep[i]);
        key += "|";
        key += design[i];
        key += "|";
        key += inference[i];
        key += "|";
        key += inference_type[i];
        results[i - i_start] = global_key_store->find(key) != global_key_store->end();
    }
    return results;
}

// [[Rcpp::export]]
int result_key_store_size_cpp() {
    return global_key_store ? global_key_store->size() : 0;
}
