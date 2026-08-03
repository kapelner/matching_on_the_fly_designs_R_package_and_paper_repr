"""Parity test for edi_kernels.fast_beta_regression against an R-generated
fixture.

This does NOT re-verify the statistical correctness of the beta-regression
fit itself -- that's the same fast_beta_regression_internal object code the
R package's own test suite already exercises via fast_beta_regression_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, dict field conversion)
and the EDI_CORE_ONLY compiled path (this extension links fetched Eigen +
LBFGSpp, not the R build's RcppEigen/RcppNumerical).

Note: res["neg_loglik"] is always NaN for this kernel, same root cause as
fast_neg_bin -- ModelResult's neg_ll field is never populated in
EDI/src/fast_beta_regression.cpp (R's own $neg_loglik is computed
separately and never flows through the ModelResult struct the Python
binding reads from). Not a binding bug; this test does not assert on it.
res["dispersion"] carries the real value (R's $phi).

The expected values below were computed once via:
    EDI:::fast_beta_regression_cpp(X, y, estimate_only = FALSE)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(109).
"""
import numpy as np
import pytest

from edi_kernels import fast_beta_regression

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(109)
    n = 200
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n) * 0.3,
    ])
    eta = 0.2 + 0.5 * X[:, 1] - 0.2 * X[:, 2]
    mu = 1 / (1 + np.exp(-eta))
    phi = 8.0
    y = rng.beta(mu * phi, (1 - mu) * phi)
    y = np.clip(y, 1e-6, 1 - 1e-6)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_beta_regression_cpp(X, y, estimate_only = FALSE)
R_B = np.array([0.156834383591728, 0.635923722518058, 0.191089027510867])
R_DISPERSION = 8.75339119583672  # R's $phi


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_beta_regression(X, y, estimate_only=False)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["dispersion"] == pytest.approx(R_DISPERSION, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_beta_regression(X, y, estimate_only=False)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["dispersion"], float)


def test_estimate_only_true_still_returns_b():
    X, y = _synthetic_data()
    res = fast_beta_regression(X, y, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
