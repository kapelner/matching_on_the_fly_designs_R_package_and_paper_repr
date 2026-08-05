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

double wilcox_hl_point_estimate_result(
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXi>& w
);

struct SummarizeWithVcovResult {
    double beta_hat;
    double ssq_hat;
    double se;
    Eigen::MatrixXd vcov;
    Eigen::VectorXd std_err;
    Eigen::VectorXd z_vals;
};

SummarizeWithVcovResult ols_hc2_post_fit_result(
    const Eigen::Ref<const Eigen::MatrixXd>& X_fit,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const Eigen::Ref<const Eigen::VectorXd>& coef_hat,
    int j_treat
);

void bind_continuous(py::module_& m) {
    m.def("fast_ols", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                          const Eigen::Ref<const Eigen::VectorXd>& y,
                          std::optional<Eigen::VectorXi> fixed_idx,
                          std::optional<Eigen::VectorXd> fixed_values,
                          bool estimate_only) {
        ModelResult res = fast_ols_internal(X, y, fixed_idx, fixed_values, estimate_only);
        py::dict out;
        out["b"] = res.b;
        out["XtWX"] = res.XtWX;
        return out;
    },
    py::arg("X"), py::arg("y"),
    py::arg("fixed_idx") = py::none(),
    py::arg("fixed_values") = py::none(),
    py::arg("estimate_only") = true,
    "Fast OLS regression (coefficients only unless estimate_only=False, which also returns XtWX).");

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

    m.def("wilcox_hl_point_estimate", &wilcox_hl_point_estimate_result,
    py::arg("y"), py::arg("w"),
    "Hodges-Lehmann point estimate (median of all pairwise treatment-minus-"
    "control differences, exact for small n or via bisection selection for "
    "large n) for the two-sample Wilcoxon rank-sum problem. w is a 0/1 "
    "treatment indicator; non-finite y entries are dropped.");

    m.def("ols_hc2_post_fit", [](const Eigen::Ref<const Eigen::MatrixXd>& X_fit,
                                  const Eigen::Ref<const Eigen::VectorXd>& y,
                                  const Eigen::Ref<const Eigen::VectorXd>& coef_hat,
                                  int j_treat) {
        SummarizeWithVcovResult res = ols_hc2_post_fit_result(X_fit, y, coef_hat, j_treat);
        py::dict out;
        out["beta_hat"] = res.beta_hat;
        out["ssq_hat"] = res.ssq_hat;
        out["se"] = res.se;
        out["vcov"] = res.vcov;
        out["std_err"] = res.std_err;
        out["z_vals"] = res.z_vals;
        return out;
    },
    py::arg("X_fit"), py::arg("y"), py::arg("coef_hat"), py::arg("j_treat") = 2,
    "HC2 heteroskedasticity-robust (sandwich) standard errors for an "
    "already-fitted OLS coefficient vector. X_fit is the fitting design "
    "matrix (e.g. Lin's intercept+treatment+centered-covariates+treatment"
    "×covariate design); j_treat is the 1-indexed column whose SE/z-value "
    "is highlighted as beta_hat/se (std_err/z_vals cover every column).");
}
