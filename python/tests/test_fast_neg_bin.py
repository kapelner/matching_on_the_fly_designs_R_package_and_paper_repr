"""Parity test for edi_kernels.fast_neg_bin against an R-generated fixture.

This does NOT re-verify the statistical correctness of the negative-
binomial fit itself -- that's the same fast_neg_bin_internal object code
the R package's own test suite already exercises via fast_neg_bin_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, dict field conversion)
and the EDI_CORE_ONLY compiled path (this extension links fetched Eigen +
LBFGSpp, not the R build's RcppEigen/RcppNumerical).

IMPORTANT default mismatch found while building this fixture: R's
fast_neg_bin_cpp defaults to smart_cold_start = FALSE, but the Python
binding's smart_cold_start defaults to True (python/cpp/bindings_count.cpp).
Left as-is (not a bug -- both are valid starting strategies that converge
to the same optimum here, just via a different iteration count) but the R
fixture below was computed with smart_cold_start = TRUE explicitly so it's
directly comparable to Python's default call.

Also note: res["neg_loglik"] is always NaN for this kernel -- ModelResult's
neg_ll field is never populated in EDI/src/fast_negbin_regression.cpp (R's
own equivalent value lives in a separately-computed $logLik field that
never flows through the ModelResult struct the Python binding reads from).
Not a binding bug; this test does not assert on it. res["dispersion"]
carries the real value (R's $theta_hat).

The expected values below were computed once via:
    EDI:::fast_neg_bin_cpp(X, as.integer(y), estimate_only = FALSE, smart_cold_start = TRUE)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(108).
"""
import numpy as np
import pytest

from edi_kernels import fast_neg_bin

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(108)
    n = 200
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n),
    ])
    beta_true = np.array([0.5, 0.3, -0.2])
    mu = np.exp(X @ beta_true)
    theta = 3.0
    p_nb = theta / (theta + mu)
    y = rng.negative_binomial(theta, p_nb).astype(np.int32)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_neg_bin_cpp(X, as.integer(y), estimate_only = FALSE, smart_cold_start = TRUE)
R_B = np.array([0.362440207398258, 0.376465607321876, -0.0794374115758531])
R_DISPERSION = 6.44141227855442  # R's $theta_hat
R_ITERATIONS = 18


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_neg_bin(X, y, estimate_only=False)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["dispersion"] == pytest.approx(R_DISPERSION, abs=ATOL, rel=RTOL)
    assert res["iterations"] == R_ITERATIONS


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_neg_bin(X, y, estimate_only=False)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["dispersion"], float)


def test_estimate_only_true_still_returns_b():
    X, y = _synthetic_data()
    res = fast_neg_bin(X, y, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
