"""Parity test for edi_kernels.fast_probit_regression against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the probit fit
itself -- that's the same fast_probit_regression_internal object code the R
package's own test suite already exercises via fast_probit_regression_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, dict field conversion)
and the EDI_CORE_ONLY compiled path (this extension links fetched Eigen +
LBFGSpp, not the R build's RcppEigen/RcppNumerical).

The expected values below were computed once via:
    EDI:::fast_probit_regression_cpp(X, y, estimate_only = FALSE)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(104).
"""
import numpy as np
import pytest
from scipy.stats import norm

from edi_kernels import fast_probit_regression

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(104)
    n = 200
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n),
    ])
    beta_true = np.array([0.2, 0.8, -0.4])
    p = norm.cdf(X @ beta_true)
    y = rng.binomial(1, p).astype(float)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_probit_regression_cpp(X, y, estimate_only = FALSE)
R_B = np.array([0.138800520046567, 0.901697766675145, -0.37137322952509])
R_NEG_LOGLIK = 106.523599522208
R_ITERATIONS = 7


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_probit_regression(X, y, estimate_only=False)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)
    assert res["iterations"] == R_ITERATIONS


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_probit_regression(X, y, estimate_only=False)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert res["mu"].shape == (X.shape[0],)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_true_still_returns_b():
    X, y = _synthetic_data()
    res = fast_probit_regression(X, y, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
