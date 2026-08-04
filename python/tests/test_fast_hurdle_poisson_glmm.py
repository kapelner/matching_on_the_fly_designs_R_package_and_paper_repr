"""Parity test for edi_kernels.fast_hurdle_poisson_glmm against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the Hurdle-Poisson
GLMM fit itself -- that's the same fast_hurdle_poisson_glmm_internal object
code the R package's own test suite already exercises via
fast_hurdle_poisson_glmm_cpp. What this test covers, which nothing on the
R side can: the pybind11 binding layer itself (argument marshaling,
defaults, edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled
path (this extension links fetched Eigen + LBFGSpp, not the R build's
RcppEigen/RcppNumerical).

X includes an intercept column (same convention as fast_logistic_glmm);
j_T is the 0-based treatment column index.

*** group_id must partition observations into groups of size EXACTLY 1 or
2 (matched pairs / singleton reservoir subjects) *** -- see
test_fast_logistic_glmm.py's docstring for why (the whole GLMM/CLMM/LMM
family shares this constraint; larger groups risk silent memory
corruption).

The expected values below were computed once via:
    EDI:::fast_hurdle_poisson_glmm_cpp(X, y, group_id, j_T = 1L)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(21). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_hurdle_poisson_glmm

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(21)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    group_id = (np.repeat(np.arange(n // 2), 2) + 1).astype(np.int32)
    beta_true = np.array([0.5, 0.4, -0.2])
    b_re = np.repeat(rng.normal(scale=0.4, size=n // 2), 2)
    mu = np.exp(X @ beta_true + b_re)
    y = rng.poisson(mu).astype(float)
    y[rng.uniform(size=n) < 0.25] = 0  # extra zeros for a genuine hurdle structure
    return X, y, group_id


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_hurdle_poisson_glmm_cpp(X, y, group_id, j_T = 1L)
R_B = np.array([0.668790651794333, 0.397003368297639, -0.280651086672844])
R_LOG_SIGMA = -1.11401467283622
R_NEG_LOGLIK = 317.46109281751


def test_matches_r_fixture():
    X, y, group_id = _synthetic_data()
    res = fast_hurdle_poisson_glmm(X, y, group_id, j_T=1)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["log_sigma"] == pytest.approx(R_LOG_SIGMA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, group_id = _synthetic_data()
    res = fast_hurdle_poisson_glmm(X, y, group_id, j_T=1)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)
