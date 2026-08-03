"""Parity test for edi_kernels.fast_robust_regression against an R-generated
fixture.

This does NOT re-verify the statistical correctness of the robust (MM-
estimation) fit itself -- that's the same fast_robust_regression_internal
object code the R package's own test suite already exercises via
fast_robust_regression_cpp. What this test covers, which nothing on the R
side can: the pybind11 binding layer itself (argument marshaling, defaults,
dict field conversion) and the EDI_CORE_ONLY compiled path (this extension
links fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

The expected values below were computed once via:
    EDI:::fast_robust_regression_cpp(X, y)
in R (default method="MM"), on the exact same synthetic dataset (with a
handful of injected outliers) generated below with
numpy.random.default_rng(102).
"""
import numpy as np
import pytest

from edi_kernels import fast_robust_regression

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(102)
    n = 150
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n),
    ])
    beta_true = np.array([1.0, 0.5, -0.3])
    y = X @ beta_true + rng.normal(scale=0.5, size=n)
    outlier_idx = rng.choice(n, 5, replace=False)
    y[outlier_idx] += rng.normal(scale=8, size=5)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_robust_regression_cpp(X, y)
R_COEFFICIENTS = np.array(
    [0.95890410926024, 0.533261327370603, -0.309509629741258]
)
R_SCALE = 0.490313028232971
R_ITERATIONS = 11


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_robust_regression(X, y)

    assert res["converged"] is True
    assert res["coefficients"] == pytest.approx(R_COEFFICIENTS, abs=ATOL, rel=RTOL)
    assert res["scale"] == pytest.approx(R_SCALE, abs=ATOL, rel=RTOL)
    assert res["iterations"] == R_ITERATIONS


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_robust_regression(X, y)

    p = X.shape[1]
    assert res["coefficients"].shape == (p,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["scale"], float)
