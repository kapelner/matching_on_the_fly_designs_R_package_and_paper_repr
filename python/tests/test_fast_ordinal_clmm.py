"""Parity test for edi_kernels.fast_ordinal_clmm against an R-generated fixture.

This does NOT re-verify the statistical correctness of the ordinal CLMM fit
itself -- that's the same fast_ordinal_clmm_internal object code the R
package's own test suite already exercises via fast_ordinal_clmm_cpp. What
this test covers, which nothing on the R side can: the pybind11 binding
layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

fast_ordinal_clmm generalizes fast_ordinal_glmm to a `link` argument
("logit", "probit", "cauchit", "cloglog"); this test exercises the default
"logit" link. X has NO intercept column; y is 1-indexed; j_T is the 0-based
treatment column index.

*** group_id must partition observations into groups of size EXACTLY 1 or
2 *** -- see test_fast_logistic_glmm.py's docstring for why (shared
GLMM/CLMM/LMM-family constraint).

The expected values below were computed once via:
    EDI:::fast_ordinal_clmm_cpp(X, y, group_id, K = 3L, j_T = 0L, link = "logit")
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(41). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_ordinal_clmm

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(41)
    n = 300
    X = np.column_stack([
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    group_id = (np.repeat(np.arange(n // 2), 2) + 1).astype(np.int32)
    alpha_true = np.array([-1.0, 0.5])
    beta_true = np.array([0.5, -0.3])
    b_re = np.repeat(rng.normal(scale=0.4, size=n // 2), 2)
    eta = X @ beta_true + b_re
    y = np.zeros(n)
    for i in range(n):
        cum = 1.0 / (1.0 + np.exp(-(alpha_true - eta[i])))
        u = rng.uniform()
        y[i] = 1 + np.sum(u > cum)
    return X, y.astype(np.int32), group_id


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_ordinal_clmm_cpp(X, y, group_id, K = 3L, j_T = 0L, link = "logit")
R_B = np.array([0.49344668004228, -0.455881644568393])
R_ALPHA = np.array([-1.03927476132197, 0.367853367499999])
R_LOG_SIGMA = -2.99621926508333
R_NEG_LOGLIK = 305.063716714842


def test_matches_r_fixture():
    X, y, group_id = _synthetic_data()
    res = fast_ordinal_clmm(X, y, group_id, K=3, j_T=0, link="logit")

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["alpha"] == pytest.approx(R_ALPHA, abs=ATOL, rel=RTOL)
    assert res["log_sigma"] == pytest.approx(R_LOG_SIGMA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, group_id = _synthetic_data()
    res = fast_ordinal_clmm(X, y, group_id, K=3, j_T=0, link="logit")

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert res["alpha"].shape == (2,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)
