"""Parity test for edi_kernels.fast_logistic_glmm against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Logistic GLMM
fit itself -- that's the same fast_logistic_glmm_internal object code the R
package's own test suite already exercises via fast_logistic_glmm_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

X INCLUDES an intercept column (unlike fast_ordinal_glmm/fast_ordinal_clmm,
which do not) -- see the comment on fast_logistic_glmm_cpp's X_r parameter
in EDI/src/fast_logistic_glmm.cpp. j_T is the 0-based treatment column
index.

*** group_id must partition observations into groups of size EXACTLY 1 or
2 (matched pairs / singleton reservoir subjects) *** -- this whole GLMM/
CLMM/LMM family uses an internal fixed-size-2 Eigen type (SM2/SV2) as an
analytical shortcut for the random-intercept covariance structure; a
larger group causes silent memory corruption (a real, confirmed,
thoroughly root-caused bug from earlier development on this package). This
fixture uses all-pairs data (group size exactly 2, no singletons).

The expected values below were computed once via:
    EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 1L)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(11). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_logistic_glmm

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(11)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    group_id = (np.repeat(np.arange(n // 2), 2) + 1).astype(np.int32)
    beta_true = np.array([0.3, 0.6, -0.4])
    b_re = np.repeat(rng.normal(scale=0.5, size=n // 2), 2)
    eta = X @ beta_true + b_re
    p = 1.0 / (1.0 + np.exp(-eta))
    y = rng.binomial(1, p, n).astype(float)
    return X, y, group_id


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_logistic_glmm_cpp(X, y, group_id, j_T = 1L)
R_B = np.array([0.460728584117953, 0.407526507792535, -0.509899263656389])
R_LOG_SIGMA = -2.9958677618216
R_NEG_LOGLIK = 184.903423932874


def test_matches_r_fixture():
    X, y, group_id = _synthetic_data()
    res = fast_logistic_glmm(X, y, group_id, j_T=1)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["log_sigma"] == pytest.approx(R_LOG_SIGMA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, group_id = _synthetic_data()
    res = fast_logistic_glmm(X, y, group_id, j_T=1)

    p = X.shape[1]
    assert res["b"].shape == (p,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_invalid_group_size_is_avoided_by_construction():
    """Documents the matched-pairs-only constraint this whole family shares
    (see module docstring) -- group sizes > 2 are invalid input, not
    exercised here on purpose (would risk memory corruption, not a clean
    assertion failure)."""
    _, _, group_id = _synthetic_data()
    _, counts = np.unique(group_id, return_counts=True)
    assert set(counts.tolist()) <= {1, 2}
