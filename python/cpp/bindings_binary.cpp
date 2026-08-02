// Binds binary-outcome model-fitting kernels (logistic, probit).

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "_helper_functions_core.h"
#include "result_map_pybind.h"
#include <optional>
#include <string>

namespace py = pybind11;

enum class BinomialConstrainedLink {
  kLog = 0,
  kIdentity = 1
};

edi::ResultMap fit_constrained_binomial_cpp_impl(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    BinomialConstrainedLink link_type,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    std::optional<Eigen::VectorXd> warm_start_weights,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

edi::ResultMap fit_constrained_binomial_weighted_cpp_impl(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::VectorXd& obs_weights,
    BinomialConstrainedLink link_type,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    std::optional<Eigen::VectorXd> warm_start_weights,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

edi::ResultMap fit_constrained_binomial_with_var_cpp_impl(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    BinomialConstrainedLink link_type,
    int j,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    std::optional<Eigen::VectorXd> warm_start_weights,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

static void bind_constrained_binomial(py::module_& m, const char* name,
                                       const char* name_with_var,
                                       BinomialConstrainedLink link_type,
                                       const char* doc, const char* doc_with_var) {
    m.def(name, [link_type](const Eigen::Ref<const Eigen::MatrixXd>& X,
                             const Eigen::Ref<const Eigen::VectorXd>& y,
                             std::optional<Eigen::VectorXd> weights,
                             int maxit,
                             double tol,
                             std::optional<Eigen::VectorXi> fixed_idx,
                             std::optional<Eigen::VectorXd> fixed_values,
                             std::optional<Eigen::VectorXd> warm_start_beta,
                             bool smart_cold_start,
                             std::optional<Eigen::VectorXd> warm_start_weights,
                             std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                             bool estimate_only) {
        edi::ResultMap res = weights.has_value()
            ? fit_constrained_binomial_weighted_cpp_impl(
                  X, y, *weights, link_type, maxit, tol, fixed_idx, fixed_values,
                  warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info,
                  estimate_only)
            : fit_constrained_binomial_cpp_impl(
                  X, y, link_type, maxit, tol, fixed_idx, fixed_values,
                  warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info,
                  estimate_only);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("weights") = py::none(),
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("warm_start_weights") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    doc);

    m.def(name_with_var, [link_type](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                      const Eigen::Ref<const Eigen::VectorXd>& y,
                                      int j,
                                      int maxit,
                                      double tol,
                                      std::optional<Eigen::VectorXi> fixed_idx,
                                      std::optional<Eigen::VectorXd> fixed_values,
                                      std::optional<Eigen::VectorXd> warm_start_beta,
                                      bool smart_cold_start,
                                      std::optional<Eigen::VectorXd> warm_start_weights,
                                      std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fit_constrained_binomial_with_var_cpp_impl(
            X, y, link_type, j, maxit, tol, fixed_idx, fixed_values,
            warm_start_beta, smart_cold_start, warm_start_weights, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("j") = 2,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("warm_start_weights") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    doc_with_var);
}

ModelResult fast_logistic_regression_internal(
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

ModelResult fast_probit_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X_eigen,
    const Eigen::Ref<const Eigen::VectorXd>& y_eigen,
    const Eigen::Ref<const Eigen::VectorXd>& weights_eigen,
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

static py::dict model_result_to_dict(const ModelResult& res) {
    py::dict out;
    out["b"] = res.b;
    out["mu"] = res.mu;
    out["XtWX"] = res.XtWX;
    out["score"] = res.score;
    out["neg_loglik"] = res.neg_ll;
    out["ssq_b_j"] = res.ssq_b_j;
    out["ssq_b_2"] = res.ssq_b_2;
    out["dispersion"] = res.dispersion;
    out["sigma2_hat"] = res.sigma2_hat;
    out["iterations"] = res.iterations;
    out["converged"] = res.converged;
    out["gradient_norm"] = res.gradient_norm;
    return out;
}

void bind_binary(py::module_& m) {
    m.def("fast_logistic_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
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
        ModelResult res = fast_logistic_regression_internal(
            X, y, weights.value_or(Eigen::VectorXd()), warm_start_beta, smart_cold_start, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
        return model_result_to_dict(res);
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
    "Fast logistic regression via IRLS or L-BFGS.");

    m.def("fast_probit_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
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
        ModelResult res = fast_probit_regression_internal(
            X, y, weights.value_or(Eigen::VectorXd()), warm_start_beta, smart_cold_start, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_weights, warm_start_fisher_info, estimate_only);
        return model_result_to_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("weights") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "irls",
    py::arg("warm_start_weights") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Fast probit regression via IRLS or L-BFGS.");

    bind_constrained_binomial(m, "fast_log_binomial_regression", "fast_log_binomial_regression_with_var",
        BinomialConstrainedLink::kLog,
        "Fast log-binomial (relative-risk) regression via Fisher scoring. Pass "
        "weights to fit the weighted variant.",
        "Fast log-binomial regression with full variance-covariance matrix (ssq_b_j "
        "is the variance of parameter index j).");

    bind_constrained_binomial(m, "fast_identity_binomial_regression", "fast_identity_binomial_regression_with_var",
        BinomialConstrainedLink::kIdentity,
        "Fast identity-link (risk-difference) binomial regression via Fisher "
        "scoring. Pass weights to fit the weighted variant.",
        "Fast identity-link binomial regression with full variance-covariance "
        "matrix (ssq_b_j is the variance of parameter index j).");
}
