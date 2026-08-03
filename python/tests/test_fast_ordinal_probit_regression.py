"""Parity test for edi_kernels.fast_ordinal_probit_regression against an
R-generated fixture. See test_fast_ordinal_regression.py for the rationale
this test class shares (binding-layer + EDI_CORE_ONLY path coverage, not
re-verifying the fit's statistical correctness).

The expected values below were computed once via:
    EDI:::fast_ordinal_probit_regression_cpp(X, y)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(123). Do not regenerate this fixture casually.
"""
import numpy as np
import pytest

from edi_kernels import fast_ordinal_probit_regression

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


# R reference, EDI 1.0.0, computed 2026-08-03 via:
#   EDI:::fast_ordinal_probit_regression_cpp(X, y)
R_B = np.array([0.431103435204438, -0.233086633969437])
R_ALPHA = np.array([-0.603651825585411, 0.22385857570018, 1.03105594451497])
R_NEG_LOGLIK = 401.07098254633


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_ordinal_probit_regression(X, y)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["alpha"] == pytest.approx(R_ALPHA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_ordinal_probit_regression(X, y)

    assert res["b"].shape == (X.shape[1],)
    assert res["alpha"].shape == (3,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_matches_full_b():
    X, y = _synthetic_data()
    res = fast_ordinal_probit_regression(X, y, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
