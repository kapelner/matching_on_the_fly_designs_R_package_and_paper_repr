"""Parity test for edi_kernels.fast_weibull_frailty against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Weibull AFT
frailty fit itself -- that's the same fast_weibull_frailty_internal object
code the R package's own test suite already exercises via
fast_weibull_frailty_cpp. What this test covers, which nothing on the R
side can: the pybind11 binding layer itself (argument marshaling, defaults,
edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled path (this
extension links fetched Eigen + LBFGSpp, not the R build's RcppEigen/
RcppNumerical).

X, y (event time), dead (event indicator) follow the same convention as
test_fast_weibull_regression.py; group_id adds the random-intercept
(frailty) grouping.

*** group_id must partition observations into groups of size EXACTLY 1 or
2 *** -- see test_fast_logistic_glmm.py's docstring for why (shared
GLMM/CLMM/LMM-family constraint).

The expected values below were computed once via:
    EDI:::fast_weibull_frailty_cpp(X, y, dead, group_id)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(51). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_weibull_frailty

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(51)
    n = 300
    X = np.column_stack([
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    group_id = (np.repeat(np.arange(n // 2), 2) + 1).astype(np.int32)
    beta_true = np.array([0.4, -0.3])
    b_re = np.repeat(rng.normal(scale=0.3, size=n // 2), 2)
    eta = X @ beta_true + b_re
    scale = np.exp(-eta)
    y = rng.exponential(scale)
    censor_time = rng.exponential(3.0, n)
    dead = (y <= censor_time).astype(float)
    y_obs = np.minimum(y, censor_time)
    return X, y_obs, dead, group_id


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_weibull_frailty_cpp(X, y, dead, group_id)
R_B = np.array([-0.316123234717711, 0.23856536066327])
R_LOG_SIGMA_EPS = 0.0230973368182207
R_LOG_SIGMA_U = -4.70307170180294
R_NEG_LOGLIK = 200.386107163062


def test_matches_r_fixture():
    X, y, dead, group_id = _synthetic_data()
    res = fast_weibull_frailty(X, y, dead, group_id)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["log_sigma_eps"] == pytest.approx(R_LOG_SIGMA_EPS, abs=ATOL, rel=RTOL)
    assert res["log_sigma_u"] == pytest.approx(R_LOG_SIGMA_U, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, dead, group_id = _synthetic_data()
    res = fast_weibull_frailty(X, y, dead, group_id)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)
