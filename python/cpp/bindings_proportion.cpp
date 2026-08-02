// Binds proportion-outcome model-fitting kernels (beta regression, zero-one-inflated beta).

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "_helper_functions_core.h"
#include <optional>
#include <string>

namespace py = pybind11;

ModelResult fast_beta_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::VectorXd* weights,
    const Eigen::VectorXd* warm_start_beta,
    bool smart_cold_start,
    double start_phi,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

LikelihoodFitResult fast_zero_one_inflated_beta_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::MatrixXd>& X_zero_one,
    const Eigen::Ref<const Eigen::VectorXd>& y_eigen,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

void bind_proportion(py::module_& m) {
    m.def("fast_beta_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                      const Eigen::Ref<const Eigen::VectorXd>& y,
                                      std::optional<Eigen::VectorXd> weights,
                                      std::optional<Eigen::VectorXd> warm_start_beta,
                                      bool smart_cold_start,
                                      double start_phi,
                                      std::optional<Eigen::VectorXi> fixed_idx,
                                      std::optional<Eigen::VectorXd> fixed_values,
                                      std::string optimization_alg,
                                      std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                      bool estimate_only) {
        ModelResult res = fast_beta_regression_internal(
            X, y,
            weights.has_value() ? &(*weights) : nullptr,
            warm_start_beta.has_value() ? &(*warm_start_beta) : nullptr,
            smart_cold_start, start_phi,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);
        py::dict out;
        out["b"] = res.b;
        out["mu"] = res.mu;
        out["XtWX"] = res.XtWX;
        out["converged"] = res.converged;
        out["neg_loglik"] = res.neg_ll;
        out["dispersion"] = res.dispersion;
        out["iterations"] = res.iterations;
        out["gradient_norm"] = res.gradient_norm;
        return out;
    },
    py::arg("X"), py::arg("y"),
    py::arg("weights") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("start_phi") = 10.0,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Fast beta regression via L-BFGS. Returned 'b' is [beta, log_phi].");

    m.def("fast_zero_one_inflated_beta", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                             const Eigen::Ref<const Eigen::MatrixXd>& X_zero_one,
                                             const Eigen::Ref<const Eigen::VectorXd>& y,
                                             std::optional<Eigen::VectorXd> warm_start_params,
                                             bool smart_cold_start,
                                             std::optional<Eigen::VectorXi> fixed_idx,
                                             std::optional<Eigen::VectorXd> fixed_values,
                                             std::string optimization_alg,
                                             std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        LikelihoodFitResult fit = fast_zero_one_inflated_beta_internal(
            X, X_zero_one, y, warm_start_params, smart_cold_start,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        py::dict out;
        out["params"] = fit.params;
        out["neg_loglik"] = fit.value;
        out["iterations"] = fit.niter;
        out["converged"] = fit.converged;
        out["gradient_norm"] = fit.gradient_norm;
        return out;
    },
    py::arg("X"), py::arg("X_zero_one"), py::arg("y"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast zero-one-inflated beta regression via L-BFGS. Returned 'params' packs "
    "[beta, log_phi, zero_one_inflation_params]; see EDI/src/fast_zero_one_inflated_beta.cpp.");
}
