// Binds survival-outcome model-fitting kernels (Cox PH, Weibull AFT).

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "result_map_pybind.h"
#include <optional>
#include <string>

namespace py = pybind11;

edi::ResultMap fast_coxph_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& dead,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    bool estimate_only,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_stratified_coxph_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& dead,
    const Eigen::Ref<const Eigen::VectorXi>& strata,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    bool estimate_only,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_weibull_frailty_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& dead,
    const Eigen::Ref<const Eigen::VectorXi>& group_id,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool estimate_only,
    int n_gh,
    double max_abs_log_sigma,
    int maxit,
    double eps_g,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_clayton_weibull_aft_optim_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& dead,
    const Eigen::Ref<const Eigen::MatrixXi>& pair_idx,
    const Eigen::Ref<const Eigen::VectorXi>& singleton_rows,
    const Eigen::Ref<const Eigen::VectorXd>& warm_start_params,
    bool estimate_only,
    int maxit,
    double reltol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_dep_cens_transform_optim_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& dead,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    bool estimate_only,
    int maxit,
    double reltol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_weibull_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& dead,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    bool estimate_only,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

void bind_survival(py::module_& m) {
    m.def("fast_coxph_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                       const Eigen::Ref<const Eigen::VectorXd>& y,
                                       const Eigen::Ref<const Eigen::VectorXd>& dead,
                                       std::optional<Eigen::VectorXd> warm_start_beta,
                                       bool smart_cold_start,
                                       bool estimate_only,
                                       int maxit,
                                       double tol,
                                       std::optional<Eigen::VectorXi> fixed_idx,
                                       std::optional<Eigen::VectorXd> fixed_values,
                                       std::string optimization_alg,
                                       std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_coxph_regression_internal(
            X, y, dead, warm_start_beta, smart_cold_start, estimate_only, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("dead"),
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("maxit") = 20,
    py::arg("tol") = 1e-9,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "newton_raphson",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast Cox proportional-hazards regression via Newton-Raphson (unstratified, "
    "no cluster-robust vcov).");

    m.def("fast_stratified_coxph_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                  const Eigen::Ref<const Eigen::VectorXd>& y,
                                                  const Eigen::Ref<const Eigen::VectorXd>& dead,
                                                  const Eigen::Ref<const Eigen::VectorXi>& strata,
                                                  std::optional<Eigen::VectorXd> warm_start_beta,
                                                  bool smart_cold_start,
                                                  bool estimate_only,
                                                  int maxit,
                                                  double tol,
                                                  std::optional<Eigen::VectorXi> fixed_idx,
                                                  std::optional<Eigen::VectorXd> fixed_values,
                                                  std::string optimization_alg,
                                                  std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_stratified_coxph_regression_internal(
            X, y, dead, strata, warm_start_beta, smart_cold_start, estimate_only, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("dead"), py::arg("strata"),
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("maxit") = 20,
    py::arg("tol") = 1e-9,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "newton_raphson",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast stratified Cox proportional-hazards regression via Newton-Raphson: one shared "
    "beta fit jointly across strata-specific baseline hazards. strata: integer group "
    "labels, any values (grouped internally, not required to be 0-based/contiguous).");

    m.def("fast_weibull_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                        const Eigen::Ref<const Eigen::VectorXd>& y,
                                        const Eigen::Ref<const Eigen::VectorXd>& dead,
                                        std::optional<Eigen::VectorXd> warm_start_params,
                                        bool smart_cold_start,
                                        bool estimate_only,
                                        int maxit,
                                        double tol,
                                        std::optional<Eigen::VectorXi> fixed_idx,
                                        std::optional<Eigen::VectorXd> fixed_values,
                                        std::string optimization_alg,
                                        std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_weibull_regression_internal(
            X, y, dead, warm_start_params, smart_cold_start, estimate_only, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("dead"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast Weibull AFT regression via L-BFGS. Returned 'params' = [beta, log_sigma].");

    m.def("fast_weibull_frailty", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                      const Eigen::Ref<const Eigen::VectorXd>& y,
                                      const Eigen::Ref<const Eigen::VectorXd>& dead,
                                      const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                      std::optional<Eigen::VectorXd> warm_start_params,
                                      std::optional<Eigen::VectorXd> warm_start_beta,
                                      bool estimate_only,
                                      int n_gh,
                                      double max_abs_log_sigma,
                                      int maxit,
                                      double eps_g,
                                      std::optional<Eigen::VectorXi> fixed_idx,
                                      std::optional<Eigen::VectorXd> fixed_values,
                                      std::string optimization_alg,
                                      std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_weibull_frailty_internal(
            X, y, dead, group_id, warm_start_params, warm_start_beta, estimate_only,
            n_gh, max_abs_log_sigma, maxit, eps_g, fixed_idx, fixed_values,
            optimization_alg, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("dead"), py::arg("group_id"),
    py::arg("warm_start_params") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("estimate_only") = false,
    py::arg("n_gh") = 20,
    py::arg("max_abs_log_sigma") = 8.0,
    py::arg("maxit") = 300,
    py::arg("eps_g") = 1e-6,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast Weibull AFT frailty regression (random intercept, Gauss-Hermite quadrature) "
    "via L-BFGS.");

    m.def("fast_clayton_weibull_aft_optim", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                const Eigen::Ref<const Eigen::VectorXd>& y,
                                                const Eigen::Ref<const Eigen::VectorXd>& dead,
                                                const Eigen::Ref<const Eigen::MatrixXi>& pair_idx,
                                                const Eigen::Ref<const Eigen::VectorXi>& singleton_rows,
                                                const Eigen::Ref<const Eigen::VectorXd>& warm_start_params,
                                                bool estimate_only,
                                                int maxit,
                                                double reltol,
                                                std::optional<Eigen::VectorXi> fixed_idx,
                                                std::optional<Eigen::VectorXd> fixed_values,
                                                std::string optimization_alg,
                                                std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_clayton_weibull_aft_optim_internal(
            X, y, dead, pair_idx, singleton_rows, warm_start_params, estimate_only,
            maxit, reltol, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("dead"), py::arg("pair_idx"), py::arg("singleton_rows"),
    py::arg("warm_start_params"),
    py::arg("estimate_only") = false,
    py::arg("maxit") = 2000,
    py::arg("reltol") = 1e-9,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast Clayton-copula Weibull AFT regression (matched pairs + singleton reservoir "
    "design) via L-BFGS. pair_idx: 0-based (n_pairs x 2) row indices into X/y/dead; "
    "singleton_rows: 0-based row indices of reservoir singletons. warm_start_params "
    "is required (no default cold start).");

    m.def("fast_dep_cens_transform_optim", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                               const Eigen::Ref<const Eigen::VectorXd>& y,
                                               const Eigen::Ref<const Eigen::VectorXd>& dead,
                                               std::optional<Eigen::VectorXd> warm_start_params,
                                               bool smart_cold_start,
                                               bool estimate_only,
                                               int maxit,
                                               double reltol,
                                               std::optional<Eigen::VectorXi> fixed_idx,
                                               std::optional<Eigen::VectorXd> fixed_values,
                                               std::string optimization_alg,
                                               std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_dep_cens_transform_optim_internal(
            X, y, dead, warm_start_params, smart_cold_start, estimate_only, maxit, reltol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("dead"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("maxit") = 2000,
    py::arg("reltol") = 1e-9,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast dependent-censoring transformation-model regression via L-BFGS.");
}
