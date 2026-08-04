"""Parity test for edi_kernels.fast_truncated_negbin_count against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the zero-truncated
negative binomial fit itself -- that's the same
fast_truncated_negbin_count_internal object code the R package's own test
suite already exercises via fast_truncated_negbin_count_cpp. What this test
covers, which nothing on the R side can: the pybind11 binding layer itself
(argument marshaling, defaults, edi::to_py_dict field conversion) and the
EDI_CORE_ONLY compiled path (this extension links fetched Eigen + LBFGSpp,
not the R build's RcppEigen/RcppNumerical).

y must be strictly positive (zero-truncated) -- unlike fast_neg_bin, y=0 is
invalid input here, not just an unusual observation.

The expected values below were computed once via:
    EDI:::fast_truncated_negbin_count_cpp(X, y)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(91). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_truncated_negbin_count

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(91)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    mu = np.exp(0.6 + 0.3 * X[:, 1] - 0.2 * X[:, 2])
    r = 4.0
    p_nb = r / (r + mu)
    y = rng.negative_binomial(r, p_nb).astype(float)
    y[y == 0] = 1.0  # zero-truncated: no zeros allowed
    return X, y


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_truncated_negbin_count_cpp(X, y)
# params = [beta (3), log_theta]
R_PARAMS = np.array([0.261203526147079, 0.401564681977834, -0.27401774478916, 0.798572995719343])
R_NEG_LL = 466.977810330508


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_truncated_negbin_count(X, y)

    assert res["converged"] is True
    assert res["params"] == pytest.approx(R_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_ll"] == pytest.approx(R_NEG_LL, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_truncated_negbin_count(X, y)

    p = X.shape[1]
    assert res["params"].shape == (p + 1,)  # beta + log_theta
    assert isinstance(res["converged"], bool)
