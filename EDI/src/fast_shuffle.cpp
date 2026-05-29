#include <RcppEigen.h>
#include <random>
#include <chrono>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector shuffle_cpp(SEXP w_sexp) {
        NumericVector w_vec(w_sexp);
        Eigen::Map<Eigen::VectorXd> w(w_vec.begin(), w_vec.size());
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::shuffle(w.data(), w.data() + w.size(), std::default_random_engine(seed));
        return wrap(w);
}

