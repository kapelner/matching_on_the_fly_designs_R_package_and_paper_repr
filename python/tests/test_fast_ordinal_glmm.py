"""Parity test for edi_kernels.fast_ordinal_glmm against an R-generated fixture.

This does NOT re-verify the statistical correctness of the ordinal GLMM fit
itself -- that's the same fast_ordinal_glmm_internal object code the R
package's own test suite already exercises via fast_ordinal_glmm_cpp. What
this test covers, which nothing on the R side can: the pybind11 binding
layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

X has NO intercept column (per the comment on fast_ordinal_glmm_cpp's X
parameter in EDI/src/fast_ordinal_glmm.cpp -- the K-1 alpha thresholds
serve that role); y is 1-indexed (values in {1,...,K}); j_T is the 0-based
treatment column index.

*** group_id must partition observations into groups of size EXACTLY 1 or
2 *** -- see test_fast_logistic_glmm.py's docstring for why (shared
GLMM/CLMM/LMM-family constraint).

The expected values below were computed once via:
    EDI:::fast_ordinal_glmm_cpp(X, y, group_id, K = 3L, j_T = 0L)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(31). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_ordinal_glmm

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(31)
    n = 300
    X = np.column_stack([
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    group_id = (np.repeat(np.arange(n // 2), 2) + 1).astype(np.int32)
    alpha_true = np.array([-1.0, 0.5])  # K=3
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
#   EDI:::fast_ordinal_glmm_cpp(X, y, group_id, K = 3L, j_T = 0L)
R_B = np.array([0.421439074385497, -0.648357579117052])
R_ALPHA = np.array([-0.95801084589414, 0.382523315421605])
R_LOG_SIGMA = -2.99535909312164
R_NEG_LOGLIK = 304.056594986312


def test_matches_r_fixture():
    X, y, group_id = _synthetic_data()
    res = fast_ordinal_glmm(X, y, group_id, K=3, j_T=0)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["alpha"] == pytest.approx(R_ALPHA, abs=ATOL, rel=RTOL)
    assert res["log_sigma"] == pytest.approx(R_LOG_SIGMA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, group_id = _synthetic_data()
    res = fast_ordinal_glmm(X, y, group_id, K=3, j_T=0)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert res["alpha"].shape == (2,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)
