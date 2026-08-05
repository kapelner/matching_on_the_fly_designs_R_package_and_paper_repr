// Binds ordinal-outcome model-fitting kernels (adjacent-category logit,
// continuation-ratio, proportional-odds logit/probit/cauchit/cloglog, the
// GLMM-random-intercept proportional-odds variant, and stereotype logit).

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "_helper_functions_core.h"
#include "result_map_pybind.h"
#include <optional>
#include <string>
#include <vector>

namespace py = pybind11;

// External linkage in EDI/src/fast_adjacent_category_logit.cpp (moved out of
// that file's anonymous namespace specifically so this binding can call them).
std::vector<double> get_levels(const Eigen::Ref<const Eigen::VectorXd>& y);
std::vector<int> map_y_to_1K(const Eigen::Ref<const Eigen::VectorXd>& y, const std::vector<double>& levels);

struct RiditAnalysisResult {
    double mean_ridit_t;
    double mean_ridit_c;
    double estimate;
    double se;
    std::vector<double> scores;
    std::vector<int> levels;
    std::vector<double> ref_p;
};

RiditAnalysisResult fast_ridit_analysis_result(
    const Eigen::Ref<const Eigen::VectorXi>& w,
    const Eigen::Ref<const Eigen::VectorXi>& y,
    const std::string& reference
);

LikelihoodFitResult fast_adjacent_category_logit_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const std::vector<int>& y_mapped,
    int K,
    int maxit,
    double tol,
    bool smart_cold_start,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::optional<Eigen::VectorXd> warm_start_beta
);

struct ContinuationRatioFit {
    LikelihoodFitResult fit;
    Eigen::MatrixXd X_aug;
    Eigen::VectorXd z;
    int n_alpha;
    int p;
};

ContinuationRatioFit fast_continuation_ratio_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool smart_cold_start,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_ordinal_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    std::optional<Eigen::VectorXd> weights,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    int maxit,
    double tol,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

edi::ResultMap fast_ordinal_probit_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    int maxit,
    double tol,
    std::string optimization_alg,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

edi::ResultMap fast_ordinal_cauchit_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    int maxit,
    double tol,
    std::string optimization_alg,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

edi::ResultMap fast_ordinal_cloglog_regression_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    std::optional<Eigen::VectorXd> warm_start_params,
    bool smart_cold_start,
    int maxit,
    double tol,
    std::string optimization_alg,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    bool estimate_only
);

edi::ResultMap fast_ordinal_clmm_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXi>& y,
    const Eigen::Ref<const Eigen::VectorXi>& group_id,
    int K,
    int j_T,
    std::string link,
    bool estimate_only,
    int n_gh,
    double max_abs_log_sigma,
    int maxit,
    double eps_g,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::string optimization_alg,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_ordinal_glmm_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXi>& y,
    const Eigen::Ref<const Eigen::VectorXi>& group_id,
    int K,
    int j_T,
    bool smart_cold_start,
    bool estimate_only,
    int n_gh,
    double max_abs_log_sigma,
    int maxit,
    double eps_g,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::optional<Eigen::VectorXd> warm_start_beta,
    std::string optimization_alg,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info
);

edi::ResultMap fast_stereotype_logit_full_internal(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    int maxit,
    double tol,
    bool smart_cold_start,
    std::optional<Eigen::VectorXi> fixed_idx,
    std::optional<Eigen::VectorXd> fixed_values,
    std::string optimization_alg,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    std::optional<Eigen::VectorXd> warm_start_params,
    std::optional<Eigen::VectorXd> warm_start_beta,
    bool estimate_only
);

void bind_ordinal(py::module_& m) {
    m.def("fast_adjacent_category_logit", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                              const Eigen::Ref<const Eigen::VectorXd>& y,
                                              int maxit,
                                              double tol,
                                              bool smart_cold_start,
                                              std::optional<Eigen::VectorXi> fixed_idx,
                                              std::optional<Eigen::VectorXd> fixed_values,
                                              std::string optimization_alg,
                                              std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                              std::optional<Eigen::VectorXd> warm_start_params,
                                              std::optional<Eigen::VectorXd> warm_start_beta) {
        std::vector<double> levels = get_levels(y);
        int K = static_cast<int>(levels.size());
        if (K < 2) {
            throw std::invalid_argument("Adjacent-category logits require at least two observed outcome categories.");
        }
        std::vector<int> y_mapped = map_y_to_1K(y, levels);

        LikelihoodFitResult fit = fast_adjacent_category_logit_internal(
            X, y_mapped, K, maxit, tol, smart_cold_start,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, warm_start_params, warm_start_beta);

        py::dict out;
        out["b"] = Eigen::VectorXd(fit.params.tail(X.cols()));
        out["alpha"] = Eigen::VectorXd(fit.params.head(K - 1));
        out["params"] = fit.params;
        out["neg_loglik"] = fit.value;
        out["iterations"] = fit.niter;
        out["converged"] = fit.converged;
        out["gradient_norm"] = fit.gradient_norm;
        return out;
    },
    py::arg("X"), py::arg("y"),
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("smart_cold_start") = true,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("warm_start_params") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    "Fast adjacent-category logit ordinal regression via L-BFGS.");

    m.def("fast_continuation_ratio_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                    const Eigen::Ref<const Eigen::VectorXd>& y,
                                                    int maxit,
                                                    double tol,
                                                    std::optional<Eigen::VectorXd> warm_start_beta,
                                                    bool smart_cold_start,
                                                    std::optional<Eigen::VectorXi> fixed_idx,
                                                    std::optional<Eigen::VectorXd> fixed_values,
                                                    std::string optimization_alg,
                                                    std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        ContinuationRatioFit cr = fast_continuation_ratio_internal(
            X, y, maxit, tol, warm_start_beta, smart_cold_start,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info);
        py::dict out;
        out["b"] = Eigen::VectorXd(cr.fit.params.tail(cr.p));
        out["alpha"] = Eigen::VectorXd(cr.fit.params.head(cr.n_alpha));
        out["params"] = cr.fit.params;
        out["neg_loglik"] = cr.fit.value;
        out["iterations"] = cr.fit.niter;
        out["converged"] = cr.fit.converged;
        out["gradient_norm"] = cr.fit.gradient_norm;
        return out;
    },
    py::arg("X"), py::arg("y"),
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("warm_start_beta") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast continuation-ratio ordinal logit regression via L-BFGS.");

    m.def("fast_ordinal_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                         const Eigen::Ref<const Eigen::VectorXd>& y,
                                         std::optional<Eigen::VectorXd> weights,
                                         std::optional<Eigen::VectorXd> warm_start_params,
                                         bool smart_cold_start,
                                         int maxit,
                                         double tol,
                                         std::optional<Eigen::VectorXi> fixed_idx,
                                         std::optional<Eigen::VectorXd> fixed_values,
                                         std::string optimization_alg,
                                         std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                         bool estimate_only) {
        edi::ResultMap res = fast_ordinal_regression_internal(
            X, y, weights, warm_start_params, smart_cold_start, maxit, tol,
            fixed_idx, fixed_values, optimization_alg, warm_start_fisher_info, estimate_only);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("weights") = py::none(),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-6,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Fast proportional-odds (logit link) ordinal regression via Newton-Raphson/L-BFGS. "
    "Returns vcov/ssq_b_j (variance of the first covariate after the alphas) whenever "
    "the Hessian is invertible and !estimate_only.");

    m.def("fast_ordinal_probit_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                const Eigen::Ref<const Eigen::VectorXd>& y,
                                                std::optional<Eigen::VectorXd> warm_start_params,
                                                bool smart_cold_start,
                                                int maxit,
                                                double tol,
                                                std::string optimization_alg,
                                                std::optional<Eigen::VectorXi> fixed_idx,
                                                std::optional<Eigen::VectorXd> fixed_values,
                                                std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                                bool estimate_only) {
        edi::ResultMap res = fast_ordinal_probit_regression_internal(
            X, y, warm_start_params, smart_cold_start, maxit, tol, optimization_alg,
            fixed_idx, fixed_values, warm_start_fisher_info, estimate_only);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-6,
    py::arg("optimization_alg") = "lbfgs",
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Fast proportional-odds (probit link) ordinal regression via Newton-Raphson/L-BFGS.");

    m.def("fast_ordinal_cauchit_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                 const Eigen::Ref<const Eigen::VectorXd>& y,
                                                 std::optional<Eigen::VectorXd> warm_start_params,
                                                 bool smart_cold_start,
                                                 int maxit,
                                                 double tol,
                                                 std::string optimization_alg,
                                                 std::optional<Eigen::VectorXi> fixed_idx,
                                                 std::optional<Eigen::VectorXd> fixed_values,
                                                 std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                                 bool estimate_only) {
        edi::ResultMap res = fast_ordinal_cauchit_regression_internal(
            X, y, warm_start_params, smart_cold_start, maxit, tol, optimization_alg,
            fixed_idx, fixed_values, warm_start_fisher_info, estimate_only);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-6,
    py::arg("optimization_alg") = "lbfgs",
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Fast proportional-odds (cauchit link) ordinal regression via Newton-Raphson/L-BFGS. "
    "Returns an empty dict if the outcome has fewer than 2 observed categories.");

    m.def("fast_ordinal_cloglog_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                                 const Eigen::Ref<const Eigen::VectorXd>& y,
                                                 std::optional<Eigen::VectorXd> warm_start_params,
                                                 bool smart_cold_start,
                                                 int maxit,
                                                 double tol,
                                                 std::string optimization_alg,
                                                 std::optional<Eigen::VectorXi> fixed_idx,
                                                 std::optional<Eigen::VectorXd> fixed_values,
                                                 std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                                 bool estimate_only) {
        edi::ResultMap res = fast_ordinal_cloglog_regression_internal(
            X, y, warm_start_params, smart_cold_start, maxit, tol, optimization_alg,
            fixed_idx, fixed_values, warm_start_fisher_info, estimate_only);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("warm_start_params") = py::none(),
    py::arg("smart_cold_start") = true,
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-6,
    py::arg("optimization_alg") = "lbfgs",
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("estimate_only") = false,
    "Fast proportional-odds (cloglog link) ordinal regression via Newton-Raphson/L-BFGS. "
    "Returns an empty dict if the outcome has fewer than 2 observed categories.");

    m.def("fast_ordinal_clmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                   const Eigen::Ref<const Eigen::VectorXi>& y,
                                   const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                   int K,
                                   int j_T,
                                   std::string link,
                                   bool estimate_only,
                                   int n_gh,
                                   double max_abs_log_sigma,
                                   int maxit,
                                   double eps_g,
                                   std::optional<Eigen::VectorXd> warm_start_params,
                                   std::string optimization_alg,
                                   std::optional<Eigen::VectorXi> fixed_idx,
                                   std::optional<Eigen::VectorXd> fixed_values,
                                   std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_ordinal_clmm_internal(
            X, y, group_id, K, j_T, link, estimate_only, n_gh, max_abs_log_sigma,
            maxit, eps_g, warm_start_params, optimization_alg, fixed_idx, fixed_values,
            warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("group_id"), py::arg("K"), py::arg("j_T"),
    py::arg("link") = "logit",
    py::arg("estimate_only") = false,
    py::arg("n_gh") = 20,
    py::arg("max_abs_log_sigma") = 8.0,
    py::arg("maxit") = 300,
    py::arg("eps_g") = 1e-6,
    py::arg("warm_start_params") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast proportional-odds ordinal GLMM (random intercept, Gauss-Hermite quadrature) "
    "via L-BFGS. link is one of 'logit', 'probit', 'cauchit', 'cloglog'.");

    m.def("fast_stereotype_logit", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                       const Eigen::Ref<const Eigen::VectorXd>& y,
                                       int maxit,
                                       double tol,
                                       bool smart_cold_start,
                                       std::optional<Eigen::VectorXi> fixed_idx,
                                       std::optional<Eigen::VectorXd> fixed_values,
                                       std::string optimization_alg,
                                       std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                       std::optional<Eigen::VectorXd> warm_start_params,
                                       std::optional<Eigen::VectorXd> warm_start_beta,
                                       bool estimate_only) {
        edi::ResultMap res = fast_stereotype_logit_full_internal(
            X, y, maxit, tol, smart_cold_start, fixed_idx, fixed_values, optimization_alg,
            warm_start_fisher_info, warm_start_params, warm_start_beta, estimate_only);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"),
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    py::arg("smart_cold_start") = true,
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("optimization_alg") = "newton_raphson",
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("warm_start_params") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("estimate_only") = false,
    "Fast stereotype logit ordinal regression via Newton-Raphson. Returns ssq_b_1/"
    "ssq_b_j (variance of the first beta) via a profile-likelihood fallback when the "
    "Fisher-information diagonal entry isn't directly invertible; vcov is always None "
    "(not computed by this kernel).");

    m.def("fast_ordinal_glmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                   const Eigen::Ref<const Eigen::VectorXi>& y,
                                   const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                   int K,
                                   int j_T,
                                   bool smart_cold_start,
                                   bool estimate_only,
                                   int n_gh,
                                   double max_abs_log_sigma,
                                   int maxit,
                                   double eps_g,
                                   std::optional<Eigen::VectorXd> warm_start_params,
                                   std::optional<Eigen::VectorXd> warm_start_beta,
                                   std::string optimization_alg,
                                   std::optional<Eigen::VectorXi> fixed_idx,
                                   std::optional<Eigen::VectorXd> fixed_values,
                                   std::optional<Eigen::MatrixXd> warm_start_fisher_info) {
        edi::ResultMap res = fast_ordinal_glmm_internal(
            X, y, group_id, K, j_T, smart_cold_start, estimate_only, n_gh, max_abs_log_sigma,
            maxit, eps_g, warm_start_params, warm_start_beta, optimization_alg,
            fixed_idx, fixed_values, warm_start_fisher_info);
        return edi::to_py_dict(res);
    },
    py::arg("X"), py::arg("y"), py::arg("group_id"), py::arg("K"), py::arg("j_T"),
    py::arg("smart_cold_start") = true,
    py::arg("estimate_only") = false,
    py::arg("n_gh") = 20,
    py::arg("max_abs_log_sigma") = 8.0,
    py::arg("maxit") = 300,
    py::arg("eps_g") = 1e-6,
    py::arg("warm_start_params") = py::none(),
    py::arg("warm_start_beta") = py::none(),
    py::arg("optimization_alg") = "lbfgs",
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    "Fast proportional-odds (logit link) ordinal GLMM with a random intercept via "
    "Gauss-Hermite quadrature + L-BFGS. X must NOT include an intercept column.");

    m.def("fast_ridit_analysis", [](const Eigen::Ref<const Eigen::VectorXi>& w,
                                     const Eigen::Ref<const Eigen::VectorXi>& y,
                                     std::string reference) {
        RiditAnalysisResult res = fast_ridit_analysis_result(w, y, reference);
        py::dict out;
        out["mean_ridit_t"] = res.mean_ridit_t;
        out["mean_ridit_c"] = res.mean_ridit_c;
        out["estimate"] = res.estimate;
        out["se"] = res.se;
        out["scores"] = res.scores;
        out["levels"] = res.levels;
        out["ref_p"] = res.ref_p;
        return out;
    },
    py::arg("w"), py::arg("y"), py::arg("reference") = "control",
    "Ridit analysis (Bross 1958): assigns each subject a ridit score relative "
    "to the empirical distribution of the reference group ('control', "
    "'treatment', or 'pooled'), then compares treatment/control mean ridits. "
    "estimate = mean_ridit_t - 0.5 (centered at 0 under the null); se is the "
    "sample-variance-based SE of the treatment-arm mean ridit.");
}
