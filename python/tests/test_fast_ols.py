"""Parity test for edi_kernels.fast_ols against an R-generated fixture.

This does NOT re-verify the statistical correctness of the OLS fit itself
-- that's the same fast_ols_internal object code the R package's own test
suite already exercises via fast_ols_cpp. What this test covers, which
nothing on the R side can: the pybind11 binding layer itself (argument
marshaling, defaults, dict field conversion) and the EDI_CORE_ONLY compiled
path (this extension links fetched Eigen, not the R build's RcppEigen).

The expected values below were computed once via:
    EDI:::fast_ols_cpp(X, y)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(101).
"""
import numpy as np
import pytest

from edi_kernels import fast_ols

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(101)
    n = 150
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n),
    ])
    beta_true = np.array([1.0, 0.5, -0.3])
    y = X @ beta_true + rng.normal(scale=0.7, size=n)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_ols_cpp(X, y)
R_B = np.array([0.984509430074946, 0.583033444887532, -0.257116081275132])


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_ols(X, y)

    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_ols(X, y)

    p = X.shape[1]
    assert res["b"].shape == (p,)
