"""Parity test for edi_kernels.fast_cpoisson_combined against an R-generated fixture.

This does NOT re-verify the statistical correctness of the combined
conditional-Poisson (matched valid pairs) + Poisson (reservoir singletons)
fit itself -- that's the same fast_cpoisson_combined_internal object code
the R package's own test suite already exercises via
fast_cpoisson_combined_with_var_cpp. What this test covers, which nothing
on the R side can: the pybind11 binding layer itself (argument marshaling,
defaults, edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled
path (this extension links fetched Eigen + LBFGSpp, not the R build's
RcppEigen/RcppNumerical).

The expected values below were computed once via:
    EDI:::fast_cpoisson_combined_with_var_cpp(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r)
in R (the with_var R export is the one with a working R-facing entry point
for this kernel; it shares the same fast_cpoisson_combined_internal call as
the plain Python binding, just also computing the variance) on the exact
same synthetic dataset generated below with numpy.random.default_rng(102)
(see package_metadata/python_bindings_package_spec.md "Testing"). Do not
regenerate this fixture casually -- if it needs updating, regenerate from R
and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_cpoisson_combined

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(102)

    # Matched valid pairs.
    nd = 40
    p = 1
    X_diff_v = rng.normal(size=(nd, p)) * 0.5
    beta_true_c = np.array([0.4])
    n_k_v = rng.integers(4, 15, nd).astype(float)
    eta = X_diff_v @ beta_true_c
    prob = 1.0 / (1.0 + np.exp(-eta))
    yT_v = rng.binomial(n_k_v.astype(int), prob).astype(float)

    # Reservoir singletons.
    nR = 50
    X_r = rng.normal(size=(nR, p)) * 0.5
    w_r = rng.binomial(1, 0.5, nR).astype(float)
    mu_r = np.exp(0.3 + 0.4 * w_r + 0.4 * X_r[:, 0])
    y_r = rng.poisson(mu_r).astype(float)

    return yT_v, n_k_v, X_diff_v, y_r, w_r, X_r


# R reference, EDI 1.0.0, computed 2026-08-03 via:
#   EDI:::fast_cpoisson_combined_with_var_cpp(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r)
R_B = np.array([0.453730142216324, 0.182136088071295, 0.31087300526784])
R_NEG_LOGLIK = 360.571939250395


def test_matches_r_fixture():
    yT_v, n_k_v, X_diff_v, y_r, w_r, X_r = _synthetic_data()
    res = fast_cpoisson_combined(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    yT_v, n_k_v, X_diff_v, y_r, w_r, X_r = _synthetic_data()
    res = fast_cpoisson_combined(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r)

    assert res["b"].shape == (3,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_still_returns_b():
    yT_v, n_k_v, X_diff_v, y_r, w_r, X_r = _synthetic_data()
    res = fast_cpoisson_combined(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
