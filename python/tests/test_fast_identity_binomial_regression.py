"""Parity test for edi_kernels.fast_identity_binomial_regression against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the identity-link
(risk-difference) fit itself -- that's the same
fit_constrained_binomial_cpp_impl object code the R package's own test
suite already exercises via fast_identity_binomial_regression_cpp. What
this test covers, which nothing on the R side can: the pybind11 binding
layer itself (argument marshaling, defaults, dict field conversion, the
log/identity link dispatch in bind_constrained_binomial) and the
EDI_CORE_ONLY compiled path (this extension links fetched Eigen, not the R
build's RcppEigen).

The expected values below were computed once via:
    EDI:::fast_identity_binomial_regression_cpp(X, y)
in R (note: this R export has no estimate_only parameter -- it always
computes the full fit), on the exact same synthetic dataset generated below
with numpy.random.default_rng(106).

IMPORTANT default mismatch found while building this fixture: R's tol
default for this kernel is 1e-8 (see EDI/R/RcppExports.R), but the Python
binding's tol default is 1e-6 (python/cpp/bindings_binary.cpp,
bind_constrained_binomial) -- a real, worth-fixing drift, not intentional
(same issue as fast_log_binomial_regression). The call below passes
tol=1e-8 explicitly to match R's fixture; without that, results agree only
to ~1e-9-1e-8 territory at best, not reliably within this test's 1e-9.
"""
import numpy as np
import pytest

from edi_kernels import fast_identity_binomial_regression

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(106)
    n = 200
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n) * 0.2,
    ])
    p = 0.3 + 0.1 * X[:, 1] - 0.05 * X[:, 2]
    p = np.clip(p, 0.05, 0.95)
    y = rng.binomial(1, p).astype(float)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_identity_binomial_regression_cpp(X, y)
R_B = np.array([0.335591533171959, 0.0100650398927152, 0.245688016207728])
R_ITERATIONS = 6


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_identity_binomial_regression(X, y, estimate_only=False, tol=1e-8)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["iterations"] == R_ITERATIONS


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_identity_binomial_regression(X, y, estimate_only=False, tol=1e-8)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert res["mu_hat"].shape == (X.shape[0],)
    assert res["fisher_information"].shape == (p, p)
    assert isinstance(res["converged"], bool)
