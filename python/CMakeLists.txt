cmake_minimum_required(VERSION 3.18)
project(edi_kernels LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# EDI/src is the single source of truth for every model-fitting kernel and
# shared header -- nothing under python/ ever copies a .cpp/.h file from
# there. This include path is what lets python/cpp/*.cpp #include headers
# (and, once a given EDI/src/*.cpp has been given the EDI_CORE_ONLY
# treatment, call its *_internal functions) directly out of the R package's
# own source tree. See package_metadata/python_bindings_package_spec.md.
set(EDI_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../EDI/src)

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(OpenMP)
# symmetric_crossprod's EDI_CORE_ONLY branch (_helper_functions_core.h) calls
# cblas_dsyrk directly (the ABI-stable CBLAS entry point, not R's raw Fortran
# dsyrk_/F77_CALL convention -- see that header's own comment). Any standard
# BLAS provider (reference libblas, OpenBLAS, MKL) exports this symbol.
find_package(BLAS REQUIRED)
include(FetchContent)

# Prefer a real Eigen3 CMake package (apt's libeigen3-dev, conda-forge, a
# vendored checkout, etc. -- widely packaged, so this is the common case).
# Falls back to fetching Eigen directly from its own upstream so a from-
# scratch build never needs anything R-specific.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found via find_package(); fetching from upstream (gitlab.com/libeigen/eigen).")
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE
    )
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
endif()

# _helper_functions_core.h pulls in LBFGSpp's <optimization/LBFGS.h>.
# LBFGSpp (github.com/yixuan/LBFGSpp, MIT license) is itself a standalone,
# header-only C++ library with zero R/Rcpp dependency -- RcppNumerical just
# vendors a copy of it for R's use. A Python-only user has no R installed,
# so this fetches LBFGSpp directly from its own upstream instead of relying
# on RcppNumerical's copy (which is what an earlier version of this file
# did, and which only worked on machines that happened to have R + that
# package installed). Pin to a tagged release for reproducibility.
FetchContent_Declare(
    lbfgspp
    GIT_REPOSITORY https://github.com/yixuan/LBFGSpp.git
    GIT_TAG v0.4.0
    GIT_SHALLOW TRUE
)
set(LBFGSPP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(lbfgspp)

# EDI/src/_helper_functions_core.h does #include <optimization/LBFGS.h> --
# that "optimization/" prefix is RcppNumerical's own vendoring convention
# (inst/include/optimization/LBFGS.h), not upstream LBFGSpp's layout (which
# is just include/LBFGS.h). Rather than touch EDI/src to special-case the
# Python build's directory layout, shim it here: a symlink so the same
# #include line resolves against the fetched upstream copy too.
set(EDI_LBFGSPP_SHIM_DIR ${CMAKE_BINARY_DIR}/lbfgspp_shim)
file(MAKE_DIRECTORY ${EDI_LBFGSPP_SHIM_DIR})
if(NOT EXISTS ${EDI_LBFGSPP_SHIM_DIR}/optimization)
    file(CREATE_LINK ${lbfgspp_SOURCE_DIR}/include ${EDI_LBFGSPP_SHIM_DIR}/optimization SYMBOLIC)
endif()

option(EDI_PY_PORTABLE "Drop -march=native for portable/cross-compiled wheels" OFF)

# Building every EDI/src kernel with EDI_CORE_ONLY defined swaps out the
# Rcpp/R-specific include branches (RcppEigen.h, Rmath.h, R_ext/BLAS.h) for
# vanilla Eigen + a CBLAS declaration -- see _helper_functions_core.h and
# fast_gamma_functions.h in EDI/src. No R installation or Rmath/libRmath
# linkage is required for this build.
add_compile_definitions(EDI_CORE_ONLY)

# EDI/src/*.cpp sources compiled directly into this extension (not copied).
# Add one line here per EDI/src file whose *_internal function gets a
# Python binding, in step with the corresponding python/cpp/bindings_*.cpp
# file that declares/calls it.
set(EDI_KERNEL_SOURCES
    ${EDI_SRC_DIR}/fast_poisson_glmm.cpp
    ${EDI_SRC_DIR}/fast_ols.cpp
    ${EDI_SRC_DIR}/fast_robust_regression.cpp
    ${EDI_SRC_DIR}/fast_logistic_regression.cpp
    ${EDI_SRC_DIR}/fast_probit_regression.cpp
    ${EDI_SRC_DIR}/fast_poisson_regression.cpp
    ${EDI_SRC_DIR}/fast_negbin_regression.cpp
    ${EDI_SRC_DIR}/fast_zinb.cpp
    ${EDI_SRC_DIR}/fast_zero_augmented_poisson.cpp
    ${EDI_SRC_DIR}/fast_beta_regression.cpp
    ${EDI_SRC_DIR}/fast_zero_one_inflated_beta.cpp
    ${EDI_SRC_DIR}/fast_adjacent_category_logit.cpp
    ${EDI_SRC_DIR}/fast_continuation_ratio_regression.cpp
    ${EDI_SRC_DIR}/fast_gee.cpp
    ${EDI_SRC_DIR}/fast_coxph_regression.cpp
    ${EDI_SRC_DIR}/fast_weibull_regression.cpp
    ${EDI_SRC_DIR}/fast_gaussian_lmm.cpp
    ${EDI_SRC_DIR}/fast_logistic_glmm.cpp
    ${EDI_SRC_DIR}/fast_ordinal_regression.cpp
    ${EDI_SRC_DIR}/fast_ordinal_probit_regression.cpp
    ${EDI_SRC_DIR}/fast_ordinal_cauchit_regression.cpp
    ${EDI_SRC_DIR}/fast_ordinal_cloglog_regression.cpp
    ${EDI_SRC_DIR}/fast_ordinal_clmm.cpp
    ${EDI_SRC_DIR}/fast_stereotype_logit.cpp
    ${EDI_SRC_DIR}/fast_ordinal_glmm.cpp
    ${EDI_SRC_DIR}/fast_cpoisson_combined.cpp
    ${EDI_SRC_DIR}/fast_weibull_frailty.cpp
    ${EDI_SRC_DIR}/fast_clogit_plus_glmm.cpp
    ${EDI_SRC_DIR}/fast_survival_models_optim.cpp
    ${EDI_SRC_DIR}/fast_hurdle_poisson_glmm.cpp
    ${EDI_SRC_DIR}/fast_hurdle_negbin.cpp
    ${EDI_SRC_DIR}/fast_log_binomial_regression.cpp
)

pybind11_add_module(_core
    cpp/bindings_module.cpp
    cpp/bindings_fast_math.cpp
    cpp/bindings_glmm.cpp
    cpp/bindings_continuous.cpp
    cpp/bindings_binary.cpp
    cpp/bindings_count.cpp
    cpp/bindings_proportion.cpp
    cpp/bindings_ordinal.cpp
    cpp/bindings_incidence.cpp
    cpp/bindings_survival.cpp
    ${EDI_KERNEL_SOURCES}
)

target_include_directories(_core PRIVATE ${EDI_SRC_DIR} cpp ${EDI_LBFGSPP_SHIM_DIR} ${lbfgspp_SOURCE_DIR}/include)
if(Eigen3_FOUND)
    target_link_libraries(_core PRIVATE Eigen3::Eigen)
else()
    target_link_libraries(_core PRIVATE eigen)
endif()
if(OpenMP_CXX_FOUND)
    target_link_libraries(_core PRIVATE OpenMP::OpenMP_CXX)
endif()
target_link_libraries(_core PRIVATE ${BLAS_LIBRARIES})

target_compile_definitions(_core PRIVATE NDEBUG EIGEN_NO_DEBUG)
target_compile_options(_core PRIVATE -Wno-ignored-attributes -fno-lto)
if(NOT EDI_PY_PORTABLE)
    target_compile_options(_core PRIVATE -march=native -mtune=native)
endif()
target_compile_options(_core PRIVATE $<$<CONFIG:Release>:-O3>)

install(TARGETS _core LIBRARY DESTINATION edi_kernels)
