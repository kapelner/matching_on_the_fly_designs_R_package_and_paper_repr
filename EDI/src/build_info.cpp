#include <Rcpp.h>
#include "edi_build_flags.h"

//' Return EDI Build Information
//'
//' Returns the compiler and package build metadata that was compiled into the
//' loaded EDI shared object. This is intended for benchmark reports and
//' reproducibility audits where the installed binary's build context matters.
//'
//' @return A named list with build timestamp, R compiler configuration, EDI
//'   build environment variables, package-level compiler and linker flags, and
//'   selected compiler preprocessor macro indicators.
//' @export
//' @examples
//' info = edi_build_info_cpp()
//' info$pkg_cxxflags
// [[Rcpp::export]]
Rcpp::List edi_build_info_cpp() {
    return Rcpp::List::create(
        Rcpp::Named("capture_method") = EDI_BUILD_CAPTURE_METHOD,
        Rcpp::Named("build_timestamp") = EDI_BUILD_TIMESTAMP,
        Rcpp::Named("build_host") = EDI_BUILD_HOST,
        Rcpp::Named("r_home") = EDI_BUILD_R_HOME,
        Rcpp::Named("r_version") = EDI_BUILD_R_VERSION,
        Rcpp::Named("r_cxx20") = EDI_BUILD_R_CXX20,
        Rcpp::Named("r_cxx20std") = EDI_BUILD_R_CXX20STD,
        Rcpp::Named("r_cxx20flags") = EDI_BUILD_R_CXX20FLAGS,
        Rcpp::Named("r_shlib_openmp_cxxflags") = EDI_BUILD_R_SHLIB_OPENMP_CXXFLAGS,
        Rcpp::Named("env_edi_portable") = EDI_BUILD_ENV_EDI_PORTABLE,
        Rcpp::Named("env_edi_disable_vectorization") = EDI_BUILD_ENV_EDI_DISABLE_VECTORIZATION,
        Rcpp::Named("env_edi_native_speed") = EDI_BUILD_ENV_EDI_NATIVE_SPEED,
        Rcpp::Named("env_edi_native_lto") = EDI_BUILD_ENV_EDI_NATIVE_LTO,
        Rcpp::Named("pkg_cppflags") = EDI_BUILD_PKG_CPPFLAGS,
        Rcpp::Named("pkg_cxxflags") = EDI_BUILD_PKG_CXXFLAGS,
        Rcpp::Named("pkg_libs") = EDI_BUILD_PKG_LIBS,
        Rcpp::Named("compiler") = __VERSION__,
#ifdef __OPTIMIZE__
        Rcpp::Named("compiler_optimize_macro") = true,
#else
        Rcpp::Named("compiler_optimize_macro") = false,
#endif
#ifdef __FAST_MATH__
        Rcpp::Named("compiler_fast_math_macro") = true,
#else
        Rcpp::Named("compiler_fast_math_macro") = false,
#endif
#ifdef EIGEN_DONT_VECTORIZE
        Rcpp::Named("eigen_dont_vectorize_macro") = true
#else
        Rcpp::Named("eigen_dont_vectorize_macro") = false
#endif
    );
}
