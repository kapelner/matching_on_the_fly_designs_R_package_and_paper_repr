"""Parity test for edi_kernels.fast_ordinal_cauchit_regression against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the proportional-odds
(cauchit link) fit itself -- that's the same fast_ordinal_cauchit_regression_internal
object code the R package's own test suite already exercises via
fast_ordinal_cauchit_regression_cpp. What this test covers, which nothing
on the R side can: the pybind11 binding layer itself (argument marshaling,
defaults, edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled
path (this extension links fetched Eigen + LBFGSpp, not the R build's
RcppEigen/RcppNumerical).

X has NO intercept column; y is 1-indexed (values in {1,...,K}).

The expected values below were computed once via:
    EDI:::fast_ordinal_cauchit_regression_cpp(X, y)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(123) (same recipe as test_fast_ordinal_regression.py).
Do not regenerate this fixture casually -- if it needs updating, regenerate
from R and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_ordinal_cauchit_regression

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(123)
    n = 300
    X_full = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    alpha_true = np.array([-1.0, 0.2, 1.3])
    beta_true = np.array([0.6, -0.4])
    eta = X_full[:, 1:] @ beta_true
    y = np.zeros(n)
    for i in range(n):
        cum = 1.0 / (1.0 + np.exp(-(alpha_true - eta[i])))
        u = rng.uniform()
        y[i] = 1 + np.sum(u > cum)
    return X_full[:, 1:], y.astype(float)


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_ordinal_cauchit_regression_cpp(X, y)
R_B = np.array([0.762940225704674, -0.380031592847121])
R_ALPHA = np.array([-1.04420895747289, 0.316682491390824, 1.60973658205326])
R_NEG_LOGLIK = 400.860778031237
R_VCOV_DIAG = np.array([0.0338536734628008, 0.0163533910514878, 0.0449028286619556,
                         0.0445505718991333, 0.0118208258597687])


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_ordinal_cauchit_regression(X, y)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["alpha"] == pytest.approx(R_ALPHA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)
    assert np.diag(res["vcov"]) == pytest.approx(R_VCOV_DIAG, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_ordinal_cauchit_regression(X, y)

    p = X.shape[1]
    n_alpha = 3
    total = p + n_alpha
    assert res["b"].shape == (p,)
    assert res["alpha"].shape == (n_alpha,)
    assert res["vcov"].shape == (total, total)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_omits_vcov():
    X, y = _synthetic_data()
    res = fast_ordinal_cauchit_regression(X, y, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
    assert "vcov" not in res
