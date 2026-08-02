// Compiles directly against EDI/src/fast_gamma_functions.h (via the CMake
// include path set in ../CMakeLists.txt) -- no file under python/ is a copy
// of anything in EDI/src. EDI_CORE_ONLY (defined globally by CMakeLists.txt)
// makes this header pull in vanilla Eigen instead of RcppEigen.h/Rmath.h, so
// this translation unit has zero R/Rcpp dependency.

#include <pybind11/pybind11.h>
#include "fast_gamma_functions.h"

namespace py = pybind11;

void bind_fast_math(py::module_& m) {
    m.def("fast_pchisq_upper", &fast_pchisq_upper,
          py::arg("statistic"), py::arg("df"),
          "Upper-tail chi-squared p-value P(X > statistic) for X ~ chi-squared(df). "
          "Matches R's pchisq(statistic, df, lower.tail=FALSE) / scipy.stats.chi2.sf(statistic, df).");
}
