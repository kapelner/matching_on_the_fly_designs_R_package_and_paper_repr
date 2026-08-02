// Binds count-outcome model-fitting kernels (Poisson, NegBin, ZINB, ZAP/hurdle-Poisson).

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "_helper_functions_core.h"
#include "result_map_pybind.h"
#include <optional>
#include <string>

namespace py = pybind11;

ModelResult fast_poisson_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& weights,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::VectorXd> warm_start_weights,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

ModelResult fast_neg_bin_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXi>& y,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    int maxit,
    double eps_g,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only,
    const Eigen::VectorXd* weights
);

LikelihoodFitResult fast_zinb_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& Xc,
    const Eigen::Ref<const Eigen::MatrixXd>& Xz,
    const Eigen::Ref<const Eigen::VectorXd>& y_vec,
    std::optional<Eigen::VectorXd> warm_start_params,
    int maxit, double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    bool smart_cold_start,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

LikelihoodFitResult fast_zap_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::MatrixXd>& Xzi,
    bool is_hurdle,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_hurdle_negbin_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::MatrixXd>& X_hurdle,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    std::optional<Eigen::MatrixXd> warm_start_hurdle_fisher_info,
    bool estimate_only,
    int j
);

edi::ResultMap fast_truncated_negbin_count_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
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

edi::ResultMap fast_cpoisson_combined_internal(
    const Eigen::Ref<const Eigen::VectorXd>& yT_v,
    const Eigen::Ref<const Eigen::VectorXd>& n_k_v,
    const Eigen::Ref<const Eigen::MatrixXd>& X_diff_v,
    const Eigen::Ref<const Eigen::VectorXd>& y_r,
    const Eigen::Ref<const Eigen::VectorXd>& w_r,
    const Eigen::Ref<const Eigen::MatrixXd>& X_r,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool estimate_only
);

static py::dict likelihood_fit_result_to_dict(const LikelihoodFitResult& fit) {
    py::dict out;
    out["params"] = fit.params;
    out["neg_loglik"] = fit.value;
    out["iterations"] = fit.niter;
    out["converged"] = fit.converged;
    out["gradient_norm"] = fit.gradient_norm;
    return out;
}

void bind_count(py::module_& m) {
    m.def("fast_poisson_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                         const Eigen::Ref<const Eigen::VectorXd>& y,
                                         std::optional<Eigen::VectorXd> weights,
                                         std::optional<Eigen::VectorXd> warm_start_beta,
                                         bool smart_cold_start,
                                         int maxit,
                                         double tol,
                                         std::optional<Eigen::VectorXi> fixed_idx,
                                         std::optional<Eigen::VectorXd> fixed_values,
                                         std::string optimization_alg,
                                         std::optional<Eigen::VectorXd> warm_start_weights,
                                         std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                         bool estimate_only) {
        ModelResult res = fast_poisson_regression_internal(
            X, y, weights.value_or(Eigen::VectorXd()), warm_start_beta, smart_cold_start, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
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
    py::arg("smart_cold_start") = false,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "irls",
    py::arg("warm_start_weights") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Fast Poisson regression via IRLS or L-BFGS.");

    m.def("fast_neg_bin", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                              const Eigen::Ref<const Eigen::VectorXi>& y,
                              std::optional<Eigen::VectorXd> warm_start_params,
                              bool smart_cold_start,
                              int maxit,
                              double eps_g,
                              std::optional<Eigen::VectorXi> fixed_idx,
                              std::optional<Eigen::VectorXd> fixed_values,
                              std::string optimization_alg,
                              std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                              bool estimate_only,
                              std::optional<Eigen::VectorXd> weights) {
        ModelResult res = fast_neg_bin_internal(
            X, y, warm_start_params, smart_cold_start, maxit, eps_g,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only,
            weights.has_value() ? &(*weights) : nullptr);
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
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 1000,
    py::arg("eps_g") = 1e-6,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    py::arg("weights") = py::none(),
    "Fast negative binomial regression via L-BFGS.");

    m.def("fast_zinb", [](const Eigen::Ref<const Eigen::MatrixXd>& Xc,
                           const Eigen::Ref<const Eigen::MatrixXd>& Xz,
                           const Eigen::Ref<const Eigen::VectorXd>& y,
                           std::optional<Eigen::VectorXd> warm_start_params,
                           int maxit,
                           double tol,
                           std::optional<Eigen::VectorXi> fixed_idx,
                           std::optional<Eigen::VectorXd> fixed_values,
                           std::string optimization_alg,
                           bool smart_cold_start,
                           std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        LikelihoodFitResult fit = fast_zinb_internal(
            Xc, Xz, y, warm_start_params, maxit, tol, fixed_idx, fixed_values,
            optimization_alg, smart_cold_start, warm_start_fisher_info);
        return likelihood_fit_result_to_dict(fit);
    },
    py::arg("Xc"), py::arg("Xz"), py::arg("y"),
    py::arg("warm_start_params") = py::none(),
    py::arg("maxit") = 1000,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("smart_cold_start") = true,
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast zero-inflated negative binomial regression via L-BFGS. Returns 'params' "
    "= [beta_cond, beta_zi, log_theta]; split into components on the Python side.");

    m.def("fast_zero_augmented_poisson", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                             const Eigen::Ref<const Eigen::VectorXd>& y,
                                             const Eigen::Ref<const Eigen::MatrixXd>& Xzi,
                                             bool is_hurdle,
                                             std::optional<Eigen::VectorXd> warm_start_params,
                                             bool smart_cold_start,
                                             int maxit,
                                             double tol,
                                             std::optional<Eigen::VectorXi> fixed_idx,
                                             std::optional<Eigen::VectorXd> fixed_values,
                                             std::string optimization_alg,
                                             std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        LikelihoodFitResult fit = fast_zap_internal(
            X, y, Xzi, is_hurdle, warm_start_params, smart_cold_start, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        return likelihood_fit_result_to_dict(fit);
    },
    py::arg("X"), py::arg("y"), py::arg("Xzi"), py::arg("is_hurdle"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 1000,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast zero-inflated (is_hurdle=False) or hurdle (is_hurdle=True) Poisson "
    "regression via L-BFGS. Returns 'params' = [beta_cond, beta_zi]; split into "
    "components on the Python side.");

    m.def("fast_hurdle_negbin", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                    const Eigen::Ref<const Eigen::VectorXd>& y,
                                    const Eigen::Ref<const Eigen::MatrixXd>& X_hurdle,
                                    std::optional<Eigen::VectorXd> warm_start_params,
                                    bool smart_cold_start,
                                    int maxit,
                                    double tol,
                                    std::optional<Eigen::VectorXi> fixed_idx,
                                    std::optional<Eigen::VectorXd> fixed_values,
                                    std::string optimization_alg,
                                    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                    std::optional<Eigen::MatrixXd> warm_start_hurdle_fisher_info,
                                    bool estimate_only,
                                    int j) {
        edi::ResultMap res = fast_hurdle_negbin_internal(
            X, y, X_hurdle, warm_start_params, smart_cold_start, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info,
            warm_start_hurdle_fisher_info, estimate_only, j);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("X_hurdle"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 1000,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("warm_start_hurdle_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    py::arg("j") = 2,
    "Fast hurdle negative binomial regression via L-BFGS (count part + logistic "
    "hurdle part). Unifies the plain-fit and with-variance R exports: vcov/ssq_b_j "
    "(for parameter index j) are always computed when !estimate_only.");

    m.def("fast_truncated_negbin_count", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                             const Eigen::Ref<const Eigen::VectorXd>& y,
                                             std::optional<Eigen::VectorXd> warm_start_params,
                                             bool smart_cold_start,
                                             bool estimate_only,
                                             int maxit,
                                             double tol,
                                             std::optional<Eigen::VectorXi> fixed_idx,
                                             std::optional<Eigen::VectorXd> fixed_values,
                                             std::string optimization_alg,
                                             std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_truncated_negbin_count_internal(
            X, y, warm_start_params, smart_cold_start, estimate_only, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("maxit") = 1000,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast zero-truncated negative binomial count regression via L-BFGS.");

    m.def("fast_cpoisson_combined", [](const Eigen::Ref<const Eigen::VectorXd>& yT_v,
                                        const Eigen::Ref<const Eigen::VectorXd>& n_k_v,
                                        const Eigen::Ref<const Eigen::MatrixXd>& X_diff_v,
                                        const Eigen::Ref<const Eigen::VectorXd>& y_r,
                                        const Eigen::Ref<const Eigen::VectorXd>& w_r,
                                        const Eigen::Ref<const Eigen::MatrixXd>& X_r,
                                        int maxit,
                                        double tol,
                                        std::optional<Eigen::VectorXi> fixed_idx,
                                        std::optional<Eigen::VectorXd> fixed_values,
                                        std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                        std::optional<Eigen::VectorXd> warm_start_params,
                                        std::optional<Eigen::VectorXd> warm_start_beta,
                                        bool estimate_only) {
        edi::ResultMap res = fast_cpoisson_combined_internal(
            yT_v, n_k_v, X_diff_v, y_r, w_r, X_r, maxit, tol,
            fixed_idx, fixed_values, warm_start_fisher_info, warm_start_params,
            warm_start_beta, estimate_only);
        return edi::to_py_dict(res);
    },
    py::arg("yT_v"), py::arg("n_k_v"), py::arg("X_diff_v"),
    py::arg("y_r"), py::arg("w_r"), py::arg("X_r"),
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("warm_start_params") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("estimate_only") = false,
    "Fast combined conditional-Poisson (matched valid pairs) + Poisson (reservoir "
    "singletons) regression via Newton's method for KK designs.");
}
