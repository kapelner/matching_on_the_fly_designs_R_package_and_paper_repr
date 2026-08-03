"""Parity test for edi_kernels.fast_dep_cens_transform_optim against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the dependent-
censoring transformation-model fit itself -- that's the same
fast_dep_cens_transform_optim_internal object code the R package's own
test suite already exercises via fast_dep_cens_transform_optim_cpp. What
this test covers, which nothing on the R side can: the pybind11 binding
layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

The expected values below were computed once via:
    EDI:::fast_dep_cens_transform_optim_cpp(X, y, dead)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(105) (see package_metadata/python_bindings_package_spec.md
"Testing"). Do not regenerate this fixture casually -- if it needs
updating, regenerate from R and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_dep_cens_transform_optim

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(105)
    n = 120
    X = np.column_stack([
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    y = rng.exponential(1.0 / np.exp(0.2 * X[:, 0] + 0.1 * X[:, 1]), n)
    dead = rng.binomial(1, 0.8, n).astype(float)
    return X, y, dead


# R reference, EDI 1.0.0, computed 2026-08-03 via:
#   EDI:::fast_dep_cens_transform_optim_cpp(X, y, dead)
R_PARAMS = np.array([
    -0.386563443389422, -0.205453049292381, 2.37533966048015,
    -0.102197528898895, 0.506734651472898, 0.774679274179001,
    -0.654186559210211,
])
R_NEG_LOGLIK = 167.102935134467


def test_matches_r_fixture():
    X, y, dead = _synthetic_data()
    res = fast_dep_cens_transform_optim(X, y, dead)

    assert res["converged"] is True
    assert res["params"] == pytest.approx(R_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, dead = _synthetic_data()
    res = fast_dep_cens_transform_optim(X, y, dead)

    assert res["params"].shape == (7,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_still_returns_params():
    X, y, dead = _synthetic_data()
    res = fast_dep_cens_transform_optim(X, y, dead, estimate_only=True)

    assert res["params"] == pytest.approx(R_PARAMS, abs=1e-6, rel=1e-6)
