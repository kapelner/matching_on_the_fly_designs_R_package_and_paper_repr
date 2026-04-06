#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <vector>
#include <string>
#include <map>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector random_block_size_redraw_w_cpp(SEXP strata_keys_sexp, SEXP block_sizes_sexp, SEXP prob_T_sexp) {
    CharacterVector strata_keys = as<CharacterVector>(strata_keys_sexp);
    IntegerVector block_sizes = as<IntegerVector>(block_sizes_sexp);
    double prob_T = as<double>(prob_T_sexp);
    
    int n = strata_keys.size();
    NumericVector w(n);
    std::map<std::string, std::vector<double>> strata_blocks;
    
    static std::random_device rd;
    static std::mt19937 g(rd());
    std::uniform_int_distribution<> block_dist(0, block_sizes.size() - 1);

    for (int i = 0; i < n; ++i) {
        std::string key = as<std::string>(strata_keys[i]);
        
        if (strata_blocks[key].empty()) {
            // Randomly pick a block size
            int current_block_size = block_sizes[block_dist(g)];
            
            int n_T_per_block = (int)std::round((double)current_block_size * prob_T);
            int n_C_per_block = current_block_size - n_T_per_block;

            std::vector<double> new_block;
            new_block.reserve(current_block_size);
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
