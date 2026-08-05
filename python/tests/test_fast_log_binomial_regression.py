"""Parity test for edi_kernels.fast_log_binomial_regression against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the log-binomial
(relative-risk) fit itself -- that's the same fit_constrained_binomial_cpp_impl
object code the R package's own test suite already exercises via
fast_log_binomial_regression_cpp. What this test covers, which nothing on
the R side can: the pybind11 binding layer itself (argument marshaling,
defaults, dict field conversion, the log/identity link dispatch in
bind_constrained_binomial) and the EDI_CORE_ONLY compiled path (this
extension links fetched Eigen, not the R build's RcppEigen).

The expected values below were computed once via:
    EDI:::fast_log_binomial_regression_cpp(X, y, estimate_only = FALSE)
in R, on the exact same synthetic dataset (small risk, for IRLS stability
under the log link) generated below with numpy.random.default_rng(105).

NOTE (corrected 2026-08-05): an earlier pass through this file claimed R's
tol default here was 1e-8 versus Python's 1e-6, and called it a real
mismatch worth fixing. That was wrong -- EDI/R/RcppExports.R and
EDI/src/fast_log_binomial_regression.cpp both show R's own default is
tol=1e-6, identical to the Python binding's (python/cpp/bindings_binary.cpp,
bind_constrained_binomial). No R/Python default drift exists for this
kernel's tol. (The likely source of the earlier confusion:
fast_neg_bin_cpp, a different kernel, has an unused eps_f=1e-8 parameter
sitting next to eps_g=1e-6 in its own signature.) Every call below still
passes tol=1e-8 explicitly -- not to match a differing R default, but
because it's a tighter convergence criterion than either side's 1e-6
default, chosen for fixture numerical stability at this test's
atol=rtol=1e-9. Both the R fixture and the Python call use the same
explicit tol=1e-8, so this is an apples-to-apples comparison either way.
"""
import numpy as np
import pytest

from edi_kernels import fast_log_binomial_regression

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(105)
    n = 200
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n) * 0.3,
    ])
    eta = -1.8 + 0.3 * X[:, 1] - 0.1 * X[:, 2]
    p = np.exp(eta)
    p = np.clip(p, 0.01, 0.5)
    y = rng.binomial(1, p).astype(float)
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_log_binomial_regression_cpp(X, y, estimate_only = FALSE)
R_B = np.array([-2.36178693647407, 0.549834908505586, -0.200845985175905])
R_ITERATIONS = 12


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_log_binomial_regression(X, y, estimate_only=False, tol=1e-8)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["iterations"] == R_ITERATIONS


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_log_binomial_regression(X, y, estimate_only=False, tol=1e-8)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert res["mu_hat"].shape == (X.shape[0],)
    assert res["fisher_information"].shape == (p, p)
    assert isinstance(res["converged"], bool)


def test_estimate_only_true_still_returns_b():
    X, y = _synthetic_data()
    res = fast_log_binomial_regression(X, y, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
