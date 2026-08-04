"""Parity test for edi_kernels.fast_zinb against an R-generated fixture.

This does NOT re-verify the statistical correctness of the zero-inflated
negative binomial fit itself -- that's the same fast_zinb_internal object
code the R package's own test suite already exercises via fast_zinb_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

Xc is the conditional (count-model) design matrix, Xz is the zero-inflation
design matrix (here the same matrix, both with an intercept column).
Returned 'params' packs [beta_cond, beta_zi, log_theta].

The expected values below were computed once via:
    EDI:::fast_zinb_cpp(Xc, Xz, y)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(61). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_zinb

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(61)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    Xz = X.copy()
    mu = np.exp(0.5 + 0.3 * X[:, 1] - 0.2 * X[:, 2])
    r = 5.0
    p_nb = r / (r + mu)
    y_count = rng.negative_binomial(r, p_nb).astype(float)
    infl_p = 1.0 / (1.0 + np.exp(-(0.2 - 0.5 * Xz[:, 1])))
    zero_infl = rng.binomial(1, infl_p, n)
    y = np.where(zero_infl == 1, 0.0, y_count)
    return X, Xz, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_zinb_cpp(Xc, Xz, y)
# params = [beta_cond (3), beta_zi (3), log_theta]
R_PARAMS = np.array([
    0.548385872088908, 0.452758008008324, -0.0943453756487485,
    0.209582564701352, -0.349470389920744, -0.111976267188071,
    2.93221168997351,
])
R_NEG_LOGLIK = 413.741399696402


def test_matches_r_fixture():
    Xc, Xz, y = _synthetic_data()
    res = fast_zinb(Xc, Xz, y)

    assert res["converged"] is True
    assert res["params"] == pytest.approx(R_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    Xc, Xz, y = _synthetic_data()
    res = fast_zinb(Xc, Xz, y)

    assert res["params"].shape == (7,)  # 3 + 3 + 1
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)
