// Binds GLMM-family model-fitting kernels. fast_poisson_glmm_internal is
// *defined* in EDI/src/fast_poisson_glmm.cpp, which is compiled as its own
// source in this same CMake target (see ../CMakeLists.txt) -- it is
// declared, not redefined, here. Nothing under python/ is a copy of
// anything in EDI/src.

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "result_map_pybind.h"
#include <Eigen/Dense>
#include <optional>
#include <string>

namespace py = pybind11;

edi::ResultMap fast_poisson_glmm_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXi>& group_id,
    int j_T,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    bool estimate_only,
    int n_gh,
    int maxit,
    double eps_g,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    std::optional<Eigen::VectorXd> row_weights
);

void bind_glmm(py::module_& m) {
    m.def("fast_poisson_glmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                   const Eigen::Ref<const Eigen::VectorXd>& y,
                                   const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                   int j_T,
                                   std::optional<Eigen::VectorXd> warm_start_params,
                                   bool smart_cold_start,
                                   bool estimate_only,
                                   int n_gh,
                                   int maxit,
                                   double eps_g,
                                   std::optional<Eigen::VectorXi> fixed_idx,
                                   std::optional<Eigen::VectorXd> fixed_values,
                                   std::string optimization_alg,
                                   std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                   std::optional<Eigen::VectorXd> row_weights) {
        return edi::to_py_dict(fast_poisson_glmm_internal(
            X, y, group_id, j_T, warm_start_params, smart_cold_start,
            estimate_only, n_gh, maxit, eps_g, fixed_idx, fixed_values,
            optimization_alg, warm_start_fisher_info, row_weights));
    },
    py::arg("X"), py::arg("y"), py::arg("group_id"), py::arg("j_T"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("n_gh") = 20,
    py::arg("maxit") = 300,
    py::arg("eps_g") = 1e-6,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("row_weights") = py::none(),
    "Fit a Poisson GLMM with a random intercept via Gauss-Hermite quadrature + L-BFGS.");
}
