// Compiles directly against EDI/src's portable (EDI_CORE_ONLY-safe) scalar
// math headers -- fast_gamma_functions.h, fast_erfc.h,
// _helper_functions_core.h, ordinal_fixed_link_helpers.h -- via the CMake
// include path set in ../CMakeLists.txt. No file under python/ is a copy of
// anything in EDI/src; every function below is a thin vectorizing wrapper
// around an inline scalar function that already lives in one of those
// headers (used elsewhere by EDI's own model-fitting kernels), so there is
// nothing new to validate numerically here -- just the elementwise loop.

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "fast_gamma_functions.h"
#include "fast_erfc.h"
#include "_helper_functions_core.h"
#include "ordinal_fixed_link_helpers.h"

namespace py = pybind11;

namespace {

template <typename Func>
Eigen::VectorXd vectorize1(Func f, const Eigen::Ref<const Eigen::VectorXd>& x) {
    Eigen::VectorXd out(x.size());
    for (Eigen::Index i = 0; i < x.size(); ++i) out[i] = f(x[i]);
    return out;
}

} // namespace

void bind_fast_math(py::module_& m) {
    m.def("fast_pchisq_upper", &fast_pchisq_upper,
          py::arg("statistic"), py::arg("df"),
          "Upper-tail chi-squared p-value P(X > statistic) for X ~ chi-squared(df). "
          "Matches R's pchisq(statistic, df, lower.tail=FALSE) / scipy.stats.chi2.sf(statistic, df).");

    m.def("fast_pchisq_upper", [](const Eigen::Ref<const Eigen::VectorXd>& statistic,
                                   const Eigen::Ref<const Eigen::VectorXd>& df) {
        Eigen::VectorXd out(statistic.size());
        for (Eigen::Index i = 0; i < statistic.size(); ++i) out[i] = fast_pchisq_upper(statistic[i], df[i]);
        return out;
    },
    py::arg("statistic"), py::arg("df"),
    "Vectorized (elementwise) overload of fast_pchisq_upper.");

    m.def("fast_digamma", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return fast_digamma(v); }, x);
    }, py::arg("x"), "Vectorized digamma (elementwise). Matches scipy.special.digamma.");

    m.def("fast_trigamma", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return fast_trigamma(v); }, x);
    }, py::arg("x"), "Vectorized trigamma (elementwise). Matches scipy.special.polygamma(1, x).");

    m.def("fast_lgamma", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return fast_lgamma(v); }, x);
    }, py::arg("x"), "Vectorized log-gamma (elementwise). Matches scipy.special.gammaln.");

    m.def("fast_lbeta", [](const Eigen::Ref<const Eigen::VectorXd>& a,
                            const Eigen::Ref<const Eigen::VectorXd>& b) {
        Eigen::VectorXd out(a.size());
        for (Eigen::Index i = 0; i < a.size(); ++i) out[i] = fast_lbeta(a[i], b[i]);
        return out;
    }, py::arg("a"), py::arg("b"), "Vectorized (elementwise) log-beta. Matches scipy.special.betaln.");

    m.def("fast_dnbinom_mu", [](const Eigen::Ref<const Eigen::VectorXd>& x, double size, double mu, bool return_log) {
        Eigen::VectorXd out(x.size());
        for (Eigen::Index i = 0; i < x.size(); ++i) out[i] = fast_dnbinom_mu(x[i], size, mu, return_log);
        return out;
    }, py::arg("x"), py::arg("size"), py::arg("mu"), py::arg("return_log") = false,
    "Vectorized (elementwise, fixed size/mu) mean-parameterized negative-binomial "
    "density. Matches scipy.stats.nbinom.logpmf(x, size, size/(size+mu)) when "
    "return_log=True.");

    m.def("fast_qnorm", [](const Eigen::Ref<const Eigen::VectorXd>& p) {
        return vectorize1([](double v) { return fast_qnorm(v); }, p);
    }, py::arg("p"), "Vectorized standard normal quantile (elementwise). Matches scipy.stats.norm.ppf.");

    m.def("fast_log_pnorm", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return fast_log_pnorm(v); }, x);
    }, py::arg("x"), "Vectorized log standard normal CDF (elementwise). Matches scipy.stats.norm.logcdf.");

    m.def("fast_log_dnorm", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return fast_log_dnorm(v); }, x);
    }, py::arg("x"), "Vectorized log standard normal PDF (elementwise). Matches scipy.stats.norm.logpdf.");

    m.def("fast_erfc", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return fast_erfc(v); }, x);
    }, py::arg("x"), "Vectorized complementary error function (elementwise). Matches scipy.special.erfc.");

    m.def("pnorm_fast", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return pnorm_fast(v); }, x);
    }, py::arg("x"), "Vectorized standard normal CDF (elementwise). Matches scipy.stats.norm.cdf.");

    m.def("dnorm_fast", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return dnorm_fast(v); }, x);
    }, py::arg("x"), "Vectorized standard normal PDF (elementwise). Matches scipy.stats.norm.pdf.");

    m.def("fast_atan", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return edi_ordinal::fast_atan(v); }, x);
    }, py::arg("x"), "Vectorized arctangent (elementwise, Cephes-style minimax approximation). Matches numpy.arctan.");

    m.def("fast_log1pexp", [](const Eigen::Ref<const Eigen::VectorXd>& x) {
        return vectorize1([](double v) { return fast_log1pexp(v); }, x);
    }, py::arg("x"), "Vectorized log(1+exp(x)) (elementwise). Matches numpy.logaddexp(0, x).");
}
