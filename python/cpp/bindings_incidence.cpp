// Binds GEE (correlated-data / incidence family).

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "_helper_functions_core.h"
#include <algorithm>
#include <optional>
#include <string>
#include <vector>

namespace py = pybind11;

enum class GEEFamily { GAUSSIAN, BINOMIAL, POISSON };

struct GEEResult {
    Eigen::VectorXd beta;
    double alpha;
    Eigen::MatrixXd vcov;
    double quasi_loglik;
    Eigen::VectorXd score;
    Eigen::MatrixXd bread;
    bool converged;
    int niter;
};

GEEResult gee_pairs_singletons_cpp_impl(
    const Eigen::Ref<const Eigen::MatrixXd>& X,
    const Eigen::Ref<const Eigen::VectorXd>& y,
    const std::vector<int>& grp_start,
    const std::vector<int>& grp_size,
    GEEFamily family,
    const Eigen::Ref<const Eigen::VectorXd>& weights,
    std::optional<Eigen::VectorXd> warm_start_beta,
    std::optional<Eigen::MatrixXd> warm_start_fisher_info,
    int maxit, double tol
);

struct MNCIBounds { double lower; double upper; };
double mn_pvalue_cpp(double x_t, double n_t, double x_c, double n_c, double delta, double p_t_obs, double p_c_obs);
MNCIBounds mn_ci_internal(double x_t, double n_t, double x_c, double n_c, double p_t_obs, double p_c_obs, double alpha, double pval_epsilon);

struct NewcombeCIBounds { double lower; double upper; };
NewcombeCIBounds newcombe_independent_ci_internal(double x1, double n1, double x2, double n2, double alpha);

void bind_incidence(py::module_& m) {
    m.def("gee_pairs_singletons", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                      const Eigen::Ref<const Eigen::VectorXd>& y,
                                      const Eigen::Ref<const Eigen::VectorXi>& group_id,
                                      std::string family_str,
                                      std::optional<Eigen::VectorXd> warm_start_beta,
                                      std::optional<Eigen::MatrixXd> warm_start_fisher_info,
                                      int maxit,
                                      double tol) {
        GEEFamily family = GEEFamily::GAUSSIAN;
        if (family_str == "binomial") family = GEEFamily::BINOMIAL;
        else if (family_str == "poisson") family = GEEFamily::POISSON;

        const int n = static_cast<int>(y.size());
        std::vector<int> ord(n);
        for (int i = 0; i < n; ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](int a, int b) { return group_id[a] < group_id[b]; });

        Eigen::MatrixXd X_s(n, X.cols());
        Eigen::VectorXd y_s(n);
        for (int i = 0; i < n; ++i) {
            X_s.row(i) = X.row(ord[i]);
            y_s[i] = y[ord[i]];
        }
        std::vector<int> grp_start, grp_size;
        int prev = -1;
        for (int i = 0; i < n; ++i) {
            int g = group_id[ord[i]];
            if (g != prev) { grp_start.push_back(i); grp_size.push_back(1); prev = g; }
            else grp_size.back()++;
        }

        GEEResult res = gee_pairs_singletons_cpp_impl(
            X_s, y_s, grp_start, grp_size, family, Eigen::VectorXd(),
            warm_start_beta, warm_start_fisher_info, maxit, tol);

        py::dict out;
        out["beta"] = res.beta;
        out["alpha"] = res.alpha;
        out["vcov"] = res.vcov;
        out["quasi_loglik"] = res.quasi_loglik;
        out["score"] = res.score;
        out["fisher_information"] = res.bread;
        out["converged"] = res.converged;
        out["niter"] = res.niter;
        return out;
    },
    py::arg("X"), py::arg("y"), py::arg("group_id"), py::arg("family"),
    py::arg("warm_start_beta") = py::none(),
    py::arg("warm_start_fisher_info") = py::none(),
    py::arg("maxit") = 100,
    py::arg("tol") = 1e-8,
    "Fast GEE (singleton/pair clusters only) via Fisher scoring. "
    "family is one of 'gaussian', 'binomial', 'poisson'.");

    m.def("mn_ci", [](double x_t, double n_t, double x_c, double n_c, double p_t_obs, double p_c_obs,
                       double alpha, double pval_epsilon) {
        MNCIBounds r = mn_ci_internal(x_t, n_t, x_c, n_c, p_t_obs, p_c_obs, alpha, pval_epsilon);
        py::tuple out(2);
        out[0] = r.lower;
        out[1] = r.upper;
        return out;
    },
    py::arg("x_t"), py::arg("n_t"), py::arg("x_c"), py::arg("n_c"),
    py::arg("p_t_obs"), py::arg("p_c_obs"),
    py::arg("alpha") = 0.05,
    py::arg("pval_epsilon") = 1e-7,
    "Miettinen-Nurminen score confidence interval for a risk difference "
    "(inverts the constrained score test via bisection). Returns (lower, upper).");

    m.def("mn_pvalue", &mn_pvalue_cpp,
    py::arg("x_t"), py::arg("n_t"), py::arg("x_c"), py::arg("n_c"),
    py::arg("delta"), py::arg("p_t_obs"), py::arg("p_c_obs"),
    "Miettinen-Nurminen two-sided score p-value for testing risk difference = delta.");

    m.def("newcombe_independent_ci", [](double x1, double n1, double x2, double n2, double alpha) {
        NewcombeCIBounds r = newcombe_independent_ci_internal(x1, n1, x2, n2, alpha);
        py::tuple out(2);
        out[0] = r.lower;
        out[1] = r.upper;
        return out;
    },
    py::arg("x1"), py::arg("n1"), py::arg("x2"), py::arg("n2"),
    py::arg("alpha") = 0.05,
    "Newcombe hybrid score confidence interval for independent proportions "
    "(Method 10). Returns (lower, upper).");
}
