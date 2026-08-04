"""Parity test for edi_kernels.fast_hurdle_negbin against an R-generated fixture.

This does NOT re-verify the statistical correctness of the hurdle negative
binomial fit itself -- that's the same fast_hurdle_negbin_internal object
code the R package's own test suite already exercises via
fast_hurdle_negbin_cpp. What this test covers, which nothing on the R side
can: the pybind11 binding layer itself (argument marshaling, defaults,
edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled path (this
extension links fetched Eigen + LBFGSpp, not the R build's RcppEigen/
RcppNumerical).

X_hurdle is the design matrix for the hurdle (zero/positive) part, fit as
an INDEPENDENT logistic regression -- structurally different from
statsmodels' HurdleCountModel, which uses a censored-count likelihood for
the same split (see python_bindings_package_spec.md's TODO-4 status note).
This kernel's own result shape reflects that split explicitly: "b"/
"converged" describe the truncated-NegBin count part, "hurdle_b"/
"hurdle_converged" describe the logistic hurdle part, "theta_hat" is the
count part's dispersion -- there is no single unified "params"/"neg_loglik"
pair the way most other kernels in this suite have.

The expected values below were computed once via:
    EDI:::fast_hurdle_negbin_cpp(X, y, X_hurdle)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(81). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_hurdle_negbin

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(81)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    X_hurdle = X.copy()
    mu = np.exp(0.6 + 0.3 * X[:, 1] - 0.2 * X[:, 2])
    r = 4.0
    p_nb = r / (r + mu)
    y_count = rng.negative_binomial(r, p_nb).astype(float)
    zero_p = 1.0 / (1.0 + np.exp(-(0.1 - 0.4 * X_hurdle[:, 1])))
    is_zero = rng.binomial(1, zero_p, n)
    y = np.where(is_zero == 1, 0.0, y_count)
    return X, y, X_hurdle


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_hurdle_negbin_cpp(X, y, X_hurdle)
R_B = np.array([0.444387907923267, 0.48344186998546, -0.257048212113382])
R_HURDLE_B = np.array([-0.559785934856885, 0.653844825842017, -0.142486317032633])
R_THETA_HAT = 2.62142098732817


def test_matches_r_fixture():
    X, y, X_hurdle = _synthetic_data()
    res = fast_hurdle_negbin(X, y, X_hurdle)

    assert res["converged"] is True
    assert res["hurdle_converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["hurdle_b"] == pytest.approx(R_HURDLE_B, abs=ATOL, rel=RTOL)
    assert res["theta_hat"] == pytest.approx(R_THETA_HAT, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, X_hurdle = _synthetic_data()
    res = fast_hurdle_negbin(X, y, X_hurdle)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert res["hurdle_b"].shape == (p,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["hurdle_converged"], bool)
