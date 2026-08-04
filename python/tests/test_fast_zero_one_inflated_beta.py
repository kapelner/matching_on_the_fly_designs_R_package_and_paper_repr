"""Parity test for edi_kernels.fast_zero_one_inflated_beta against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the zero-one-
inflated beta fit itself -- that's the same
fast_zero_one_inflated_beta_internal object code the R package's own test
suite already exercises via fast_zero_one_inflated_beta_cpp. What this test
covers, which nothing on the R side can: the pybind11 binding layer itself
(argument marshaling, defaults, edi::to_py_dict field conversion) and the
EDI_CORE_ONLY compiled path (this extension links fetched Eigen + LBFGSpp,
not the R build's RcppEigen/RcppNumerical).

X_zero_one is the design matrix for the zero/one-inflation component; y is
in [0, 1] inclusive (may include exact 0s and 1s, unlike plain beta
regression). Returned 'params' packs [beta, log_phi, zero_one_b0,
zero_one_b1].

Note: fast_zero_one_inflated_beta_cpp's R-facing return list has no
"converged"/"iterations"/"gradient_norm" fields at all (verified via
`names(EDI:::fast_zero_one_inflated_beta_cpp(...))`) even though the
underlying LikelihoodFitResult the Python binding surfaces does carry
them -- the R wrapper just doesn't expose those fields. So this test only
cross-checks params/neg_loglik against R (the only fields both sides
expose) and asserts converged==True as a Python-side-only sanity check,
not a cross-language comparison.

The expected values below were computed once via:
    EDI:::fast_zero_one_inflated_beta_cpp(X, X_zero_one, y)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(101). Do not regenerate this fixture casually --
if it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_zero_one_inflated_beta

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(101)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    X_zero_one = X.copy()
    mu = 1.0 / (1.0 + np.exp(-(0.3 + 0.4 * X[:, 1] - 0.3 * X[:, 2])))
    phi = 8.0
    y_beta = rng.beta(mu * phi, (1 - mu) * phi)

    p0 = 1.0 / (1.0 + np.exp(-(-1.0 + 0.3 * X_zero_one[:, 1])))
    p1 = 1.0 / (1.0 + np.exp(-(-1.2 - 0.2 * X_zero_one[:, 1])))
    u = rng.uniform(size=n)
    y = y_beta.copy()
    y[u < p0] = 0.0
    y[(u >= p0) & (u < p0 + p1)] = 1.0
    return X, X_zero_one, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_zero_one_inflated_beta_cpp(X, X_zero_one, y)
# params = [beta (3), log_phi, zero_one_b0 (3), zero_one_b1 (3)]
R_PARAMS = np.array([
    0.352794130617712, 0.325391374378021, -0.282430219264307,
    1.88727372055993,
    -0.629974194384907, 0.00854523430697365, -0.00755452556415076,
    -0.60718765372455, 0.0122774565246871, -0.058961448979751,
])
R_NEG_LOGLIK = 262.24798255249


def test_matches_r_fixture():
    X, X_zero_one, y = _synthetic_data()
    res = fast_zero_one_inflated_beta(X, X_zero_one, y)

    assert res["params"] == pytest.approx(R_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, X_zero_one, y = _synthetic_data()
    res = fast_zero_one_inflated_beta(X, X_zero_one, y)

    assert res["params"].shape == (10,)  # 3 + 1 + 3 + 3
    assert isinstance(res["neg_loglik"], float)
    assert res["converged"] is True  # Python-side-only sanity check, see docstring
