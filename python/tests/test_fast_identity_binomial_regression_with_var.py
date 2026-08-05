"""Parity test for edi_kernels.fast_identity_binomial_regression_with_var
against an R-generated fixture.

Secondary entry point into the same file as
test_fast_identity_binomial_regression.py -- see that file's docstring and
test_fast_log_binomial_regression_with_var.py's docstring (its logit-link
sibling) for context on this pair of "_with_var" entry points.

The expected values below were computed once via:
    EDI:::fast_identity_binomial_regression_with_var_cpp(X, y, tol = 1e-8)
in R -- tol passed explicitly to match the Python call below (a
tighter-than-default convergence criterion for fixture stability, not a
workaround for an R/Python default mismatch -- both sides actually default
to tol=1e-6, see test_fast_log_binomial_regression.py's docstring), on the
exact same synthetic dataset generated below with
numpy.random.default_rng(111) (the y draw here uses a fresh continuation
of the same rng stream as test_fast_log_binomial_regression_with_var.py's
data -- keep both files' data-generation code in sync if either is
regenerated, since they intentionally share one rng stream). Do not
regenerate this fixture casually -- if it needs updating, regenerate from R
and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_identity_binomial_regression_with_var

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(111)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    # Draw (and discard) the same first y as the log-binomial sibling test
    # so this file's y draw lands on the identical second segment of the
    # shared rng stream used when the fixture was computed.
    p_small = 0.08 + 0.05 * X[:, 1]
    rng.binomial(1, p_small, n)

    p2 = np.clip(0.3 + 0.1 * X[:, 1] - 0.03 * X[:, 2], 0.02, 0.98)
    y2 = rng.binomial(1, p2, n).astype(float)
    return X, y2


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_identity_binomial_regression_with_var_cpp(X, y, tol = 1e-8)
R_B = np.array([0.249771559457616, 0.194646444716197, -0.0553245578510556])
R_NEG_LL = 184.294560022455


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_identity_binomial_regression_with_var(X, y, tol=1e-8)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["neg_ll"] == pytest.approx(R_NEG_LL, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_identity_binomial_regression_with_var(X, y, tol=1e-8)

    assert res["b"].shape == (X.shape[1],)
    assert isinstance(res["converged"], bool)
