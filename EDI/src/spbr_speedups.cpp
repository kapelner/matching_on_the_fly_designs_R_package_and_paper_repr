#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <vector>
#include <string>
#include <map>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector spbr_redraw_w_cpp(SEXP strata_keys_sexp, SEXP block_size_sexp, SEXP prob_T_sexp) {
    CharacterVector strata_keys = as<CharacterVector>(strata_keys_sexp);
    int block_size = as<int>(block_size_sexp);
    double prob_T = as<double>(prob_T_sexp);
    
    int n = strata_keys.size();
    NumericVector w(n);
    std::map<std::string, std::vector<double>> strata_blocks;
    
    int n_T_per_block = (int)std::round((double)block_size * prob_T);
    int n_C_per_block = block_size - n_T_per_block;

    static std::random_device rd;
    static std::mt19937 g(rd());

    for (int i = 0; i < n; ++i) {
        std::string key = as<std::string>(strata_keys[i]);
        
        if (strata_blocks[key].empty()) {
            std::vector<double> new_block;
            new_block.reserve(block_size);
            for (int j = 0; j < n_T_per_block; ++j) new_block.push_back(1.0);
            for (int j = 0; j < n_C_per_block; ++j) new_block.push_back(0.0);
            
            std::shuffle(new_block.begin(), new_block.end(), g);
            strata_blocks[key] = new_block;
        }
        
        w[i] = strata_blocks[key].back();
        strata_blocks[key].pop_back();
    }
    
    return w;
}

// [[Rcpp::export]]
IntegerVector stratified_bootstrap_indices_cpp(SEXP strata_keys_sexp) {
    CharacterVector strata_keys = as<CharacterVector>(strata_keys_sexp);
    int n = strata_keys.size();
    
    std::map<std::string, std::vector<int>> strata_map;
    for (int i = 0; i < n; ++i) {
        std::string key = as<std::string>(strata_keys[i]);
        strata_map[key].push_back(i + 1);
    }
    
    static std::random_device rd;
    static std::mt19937 gen(rd());
    
    IntegerVector indices(n);
    int current_idx = 0;
    for (auto const& [key, members] : strata_map) {
        int m = members.size();
        std::uniform_int_distribution<> dis(0, m - 1);
        for (int j = 0; j < m; ++j) {
            indices[current_idx++] = members[dis(gen)];
        }
    }
    
    std::shuffle(indices.begin(), indices.end(), gen);
    return indices;
}
