#include <Rcpp.h>
#include <random>
#include <chrono>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector shuffle_cpp(NumericVector w) {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(w.begin(), w.end(), std::default_random_engine(seed));
	return w;
}
