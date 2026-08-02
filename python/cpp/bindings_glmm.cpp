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

edi::ResultMap fast_gaussian_lmm_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXi>& group_id,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool estimate_only,
    int maxit,
    double eps_g,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    std::optional<Eigen::VectorXd> weights
);

edi::ResultMap fast_hurdle_poisson_glmm_internal(
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
    std::string optimization_alg,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_clogit_plus_glmm_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X_disc,
    const Eigen::Ref<const Eigen::VectorXd>& y_disc,
    const Eigen::Ref<const Eigen::MatrixXd>& X_conc,
    const Eigen::Ref<const Eigen::VectorXd>& y_conc,
    const Eigen::Ref<const Eigen::VectorXi>& group_conc,
    bool has_discordant,
    bool has_concordant,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool estimate_only,
    double max_abs_log_sigma,
    int maxit,
    double eps_g,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    int n_gh
);

edi::ResultMap fast_logistic_glmm_internal(
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
    std::string optimization_alg,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
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

    m.def("fast_gaussian_lmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                   const Eigen::Ref<const Eigen::VectorXd>& y,
                                   const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                   std::optional<Eigen::VectorXd> warm_start_params,
                                   std::optional<Eigen::VectorXd> warm_start_beta,
                                   bool estimate_only,
                                   int maxit,
                                   double eps_g,
                                   std::optional<Eigen::VectorXi> fixed_idx,
                                   std::optional<Eigen::VectorXd> fixed_values,
                                   std::string optimization_alg,
                                   std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                   std::optional<Eigen::VectorXd> weights) {
        edi::ResultMap res = fast_gaussian_lmm_internal(
            X, y, group_id, warm_start_params, warm_start_beta, estimate_only,
            maxit, eps_g, fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, weights);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("group_id"),
    py::arg("warm_start_params") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("estimate_only") = false,
    py::arg("maxit") = 300,
    py::arg("eps_g") = 1e-6,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("weights") = py::none(),
    "Fit a Gaussian LMM with a random intercept via L-BFGS. Returned 'b' is "
    "[beta, log_sigma_e, log_sigma_b] (plain array, not R's named vector).");

    m.def("fast_logistic_glmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                    const Eigen::Ref<const Eigen::VectorXd>& y,
                                    const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                    int j_T,
                                    std::optional<Eigen::VectorXd> warm_start_params,
                                    bool smart_cold_start,
                                    bool estimate_only,
                                    int n_gh,
                                    int maxit,
                                    double eps_g,
                                    std::string optimization_alg,
                                    std::optional<Eigen::VectorXi> fixed_idx,
                                    std::optional<Eigen::VectorXd> fixed_values,
                                    std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_logistic_glmm_internal(
            X, y, group_id, j_T, warm_start_params, smart_cold_start, estimate_only,
            n_gh, maxit, eps_g, optimization_alg, fixed_idx, fixed_values, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("group_id"), py::arg("j_T"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("n_gh") = 20,
    py::arg("maxit") = 300,
    py::arg("eps_g") = 1e-6,
    py::arg("optimization_alg") = "lbfgs",
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    "Fit a Logistic GLMM with a random intercept via Gauss-Hermite quadrature + L-BFGS.");

    m.def("fast_hurdle_poisson_glmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                          const Eigen::Ref<const Eigen::VectorXd>& y,
                                          const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                          int j_T,
                                          std::optional<Eigen::VectorXd> warm_start_params,
                                          bool smart_cold_start,
                                          bool estimate_only,
                                          int n_gh,
                                          int maxit,
                                          double eps_g,
                                          std::string optimization_alg,
                                          std::optional<Eigen::VectorXi> fixed_idx,
                                          std::optional<Eigen::VectorXd> fixed_values,
                                          std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_hurdle_poisson_glmm_internal(
            X, y, group_id, j_T, warm_start_params, smart_cold_start, estimate_only,
            n_gh, maxit, eps_g, optimization_alg, fixed_idx, fixed_values, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("group_id"), py::arg("j_T"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("n_gh") = 7,
    py::arg("maxit") = 300,
    py::arg("eps_g") = 1e-6,
    py::arg("optimization_alg") = "lbfgs",
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    "Fit a Hurdle-Poisson GLMM with a random intercept via Gauss-Hermite quadrature + L-BFGS.");

    m.def("fast_clogit_plus_glmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X_disc,
                                       const Eigen::Ref<const Eigen::VectorXd>& y_disc,
                                       const Eigen::Ref<const Eigen::MatrixXd>& X_conc,
                                       const Eigen::Ref<const Eigen::VectorXd>& y_conc,
                                       const Eigen::Ref<const Eigen::VectorXi>& group_conc,
                                       bool has_discordant,
                                       bool has_concordant,
                                       std::optional<Eigen::VectorXd> warm_start_params,
                                       std::optional<Eigen::VectorXd> warm_start_beta,
                                       bool estimate_only,
                                       double max_abs_log_sigma,
                                       int maxit,
                                       double eps_g,
                                       std::optional<Eigen::VectorXi> fixed_idx,
                                       std::optional<Eigen::VectorXd> fixed_values,
                                       std::string optimization_alg,
                                       std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                       int n_gh) {
        edi::ResultMap res = fast_clogit_plus_glmm_internal(
            X_disc, y_disc, X_conc, y_conc, group_conc, has_discordant, has_concordant,
            warm_start_params, warm_start_beta, estimate_only, max_abs_log_sigma, maxit, eps_g,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, n_gh);
        return edi::to_py_dict(res);
    },
    py::arg("X_disc"), py::arg("y_disc"), py::arg("X_conc"), py::arg("y_conc"), py::arg("group_conc"),
    py::arg("has_discordant"), py::arg("has_concordant"),
    py::arg("warm_start_params") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("estimate_only") = false,
    py::arg("max_abs_log_sigma") = 8.0,
    py::arg("maxit") = 200,
    py::arg("eps_g") = 1e-5,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("n_gh") = 20,
    "Fit a combined conditional-logit (discordant pairs) + random-intercept-GLMM "
    "(concordant pairs) model via Gauss-Hermite quadrature + L-BFGS.");
}
