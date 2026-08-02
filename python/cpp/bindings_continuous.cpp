// Binds continuous-outcome model-fitting kernels (OLS, robust regression).
// Each *_internal function is defined in its own EDI/src/fast_*.cpp file,
// compiled as its own source in this same CMake target (see
// ../CMakeLists.txt) -- declared, not redefined, here.

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "_helper_functions_core.h"
#include <optional>

namespace py = pybind11;

ModelResult fast_ols_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    bool estimate_only
);

// Matches EDI/src/fast_robust_regression.cpp's file-local RobustModelResult
// exactly (no shared header exists for it) -- keep in sync if that struct's
// field list ever changes.
struct RobustModelResult {
    Eigen::VectorXd b;
    Eigen::VectorXd w;
    Eigen::MatrixXd XtWX;
    Eigen::MatrixXd X_free;
    double XtX_inv_diag_j;
    double scale;
    int iterations;
    bool converged;
    double ssq_b_j;
};

RobustModelResult fast_robust_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    std::string method,
    double c,
    double c_bisquare,
    int maxit,
    double tol,
    double scale_est,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::VectorXd> warm_start_weights,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only,
    int variance_j
);

void bind_continuous(py::module_& m) {
    m.def("fast_ols", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                          const Eigen::Ref<const Eigen::VectorXd>& y,
                          std::optional<Eigen::VectorXi> fixed_idx,
                          std::optional<Eigen::VectorXd> fixed_values) {
        ModelResult res = fast_ols_internal(X, y, fixed_idx, fixed_values, true);
        py::dict out;
        out["b"] = res.b;
        return out;
    },
    py::arg("X"), py::arg("y"),
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    "Fast OLS regression (coefficients only).");

    m.def("fast_robust_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                        const Eigen::Ref<const Eigen::VectorXd>& y,
                                        std::optional<Eigen::VectorXd> warm_start_beta,
                                        bool smart_cold_start,
                                        std::string method,
                                        int j,
                                        double c,
                                        int maxit,
                                        double tol,
                                        std::optional<Eigen::VectorXi> fixed_idx,
                                        std::optional<Eigen::VectorXd> fixed_values,
                                        std::optional<Eigen::VectorXd> warm_start_weights,
                                        std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                        bool estimate_only) {
        RobustModelResult res = fast_robust_regression_internal(
            X, y, warm_start_beta, smart_cold_start, method, c, 4.685, maxit, tol, -1.0,
            fixed_idx, fixed_values, warm_start_weights, warm_start_fisher_info,
            estimate_only, j);
        py::dict out;
        out["coefficients"] = res.b;
        out["scale"] = res.scale;
        out["converged"] = res.converged;
        out["iterations"] = res.iterations;
        return out;
    },
    py::arg("X"), py::arg("y"),
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("method") = "MM",
    py::arg("j") = 2,
    py::arg("c") = 1.345,
    py::arg("maxit") = 50,
    py::arg("tol") = 1e-7,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_weights") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Robust regression via IRLS (Huber/Tukey bisquare M/MM-estimation).");
}
