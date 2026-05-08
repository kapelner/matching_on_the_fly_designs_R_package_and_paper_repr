#include <RcppEigen.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Check for MKL
#ifdef EIGEN_USE_MKL_ALL
#include <mkl.h>
#endif

using namespace Rcpp;

//' Set the number of threads for OpenMP, Eigen, and MKL
//'
//' @param n_threads Integer.
//' @keywords internal
// [[Rcpp::export]]
void set_omp_num_threads_cpp(int n_threads) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_nested(0);
    omp_set_num_threads(n_threads);
#endif
    Eigen::setNbThreads(n_threads);
#ifdef EIGEN_USE_MKL_ALL
    mkl_set_num_threads(n_threads);
#endif
}

//' Get the maximum number of threads for OpenMP
//'
//' @return Integer.
//' @keywords internal
// [[Rcpp::export]]
int get_omp_max_threads_cpp() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}
