"""Parity test for edi_kernels.fast_log_binomial_regression_with_var
against an R-generated fixture.

Secondary entry point into the same file as
test_fast_log_binomial_regression.py (fast_log_binomial_regression.cpp's
fit_constrained_binomial_with_var_cpp_impl) -- covered separately here
since it's its own bound Python function with its own arg list (adds `j`,
drops `weights`/`estimate_only`).

This does NOT re-verify the statistical correctness of the fit itself --
see test_fast_log_binomial_regression.py's docstring for what a parity
test in this suite does and doesn't cover.

The expected values below were computed once via:
    EDI:::fast_log_binomial_regression_with_var_cpp(X, y, tol = 1e-8)
in R (tol set explicitly to match the Python binding's default -- see
test_fast_log_binomial_regression.py's docstring for the R/Python tol
default mismatch this works around), on the exact same synthetic dataset
generated below with numpy.random.default_rng(111). Do not regenerate this
fixture casually -- if it needs updating, regenerate from R and update the
comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_log_binomial_regression_with_var

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(111)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    p_small = 0.08 + 0.05 * X[:, 1]
    y = rng.binomial(1, p_small, n).astype(float)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_log_binomial_regression_with_var_cpp(X, y, tol = 1e-8)
R_B = np.array([-2.34826908323323, 0.106716389582667, 0.271106901658374])
R_NEG_LL = 100.319039980906


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_log_binomial_regression_with_var(X, y, tol=1e-8)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["neg_ll"] == pytest.approx(R_NEG_LL, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_log_binomial_regression_with_var(X, y, tol=1e-8)

    assert res["b"].shape == (X.shape[1],)
    assert isinstance(res["converged"], bool)
