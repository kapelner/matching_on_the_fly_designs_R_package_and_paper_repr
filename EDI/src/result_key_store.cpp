#include <Rcpp.h>
#include <unordered_set>
#include <cstdint>
#include <cstdio>

using namespace Rcpp;

// uint64_t keys instead of std::string: no per-key heap allocation, cache-friendly lookup.
// FNV-1a 64-bit collision probability is ~N^2/2^65 — negligible for any realistic Nrep*n_cells.
static std::unordered_set<uint64_t>* global_key_store = nullptr;

static const uint64_t FNV_OFFSET = 14695981039346656037ULL;
static const uint64_t FNV_PRIME  = 1099511628211ULL;

inline uint64_t fnv1a_str(uint64_t h, const char* s) {
    while (*s) { h ^= static_cast<uint8_t>(*s++); h *= FNV_PRIME; }
    return h;
}
inline uint64_t fnv1a_sep(uint64_t h) { h ^= static_cast<uint8_t>('|'); return h * FNV_PRIME; }

// CHAR(STRING_ELT(...)) gives a direct pointer into R's string cache — zero allocation.
// snprintf onto a stack buffer for numeric fields — zero heap allocation.
inline uint64_t compute_key(
    CharacterVector& response_type, CharacterVector& cond_exp_func_model,
    IntegerVector& n, IntegerVector& p, NumericVector& betaT, IntegerVector& rep,
    CharacterVector& design, CharacterVector& inference, CharacterVector& inference_type,
    int i)
{
    char buf[64];
    uint64_t h = FNV_OFFSET;
    h = fnv1a_str(h, CHAR(STRING_ELT(response_type,        i))); h = fnv1a_sep(h);
    h = fnv1a_str(h, CHAR(STRING_ELT(cond_exp_func_model,  i))); h = fnv1a_sep(h);
    std::snprintf(buf, sizeof(buf), "%d",    n[i]);    h = fnv1a_str(h, buf); h = fnv1a_sep(h);
    std::snprintf(buf, sizeof(buf), "%d",    p[i]);    h = fnv1a_str(h, buf); h = fnv1a_sep(h);
    std::snprintf(buf, sizeof(buf), "%.15g", betaT[i]); h = fnv1a_str(h, buf); h = fnv1a_sep(h);
    std::snprintf(buf, sizeof(buf), "%d",    rep[i]);  h = fnv1a_str(h, buf); h = fnv1a_sep(h);
    h = fnv1a_str(h, CHAR(STRING_ELT(design,          i))); h = fnv1a_sep(h);
    h = fnv1a_str(h, CHAR(STRING_ELT(inference,       i))); h = fnv1a_sep(h);
    h = fnv1a_str(h, CHAR(STRING_ELT(inference_type,  i)));
    return h;
}

// [[Rcpp::export]]
void clear_result_key_store_cpp() {
    delete global_key_store;
    global_key_store = nullptr;
}

// [[Rcpp::export]]
void init_result_key_store_cpp(int expected_size) {
    clear_result_key_store_cpp();
    global_key_store = new std::unordered_set<uint64_t>();
    if (expected_size > 0) {
        global_key_store->reserve(static_cast<size_t>(expected_size));
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
        global_key_store->insert(compute_key(
            response_type, cond_exp_func_model, n, p, betaT, rep,
            design, inference, inference_type, i));
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
    LogicalVector results(i_end - i_start, false);
    if (!global_key_store) return results;
    for (int i = i_start; i < i_end; ++i) {
        results[i - i_start] = global_key_store->count(compute_key(
            response_type, cond_exp_func_model, n, p, betaT, rep,
            design, inference, inference_type, i)) > 0;
    }
    return results;
}

// [[Rcpp::export]]
int result_key_store_size_cpp() {
    return global_key_store ? static_cast<int>(global_key_store->size()) : 0;
}
