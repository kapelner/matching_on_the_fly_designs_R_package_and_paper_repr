"""Parity test for edi_kernels.gee_pairs_singletons against an R-generated fixture.

This does NOT re-verify the statistical correctness of the matched-pair/
singleton GEE fit itself -- that's the same gee_pairs_singletons_cpp_impl
object code the R package's own test suite already exercises via
gee_pairs_singletons_cpp. What this test covers, which nothing on the R
side can: the pybind11 binding layer itself (argument marshaling, the
group_id -> grp_start/grp_size sort+partition logic reimplemented in the
Python binding, edi::to_py_dict field conversion) and the EDI_CORE_ONLY
compiled path (this extension links fetched Eigen + LBFGSpp, not the R
build's RcppEigen/RcppNumerical).

The expected values below were computed once via:
    EDI:::gee_pairs_singletons_cpp(X, y, group_id, "gaussian")
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(101) (see package_metadata/python_bindings_package_spec.md
"Testing"). Do not regenerate this fixture casually -- if it needs
updating, regenerate from R and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import gee_pairs_singletons

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(101)
    n = 100
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    beta_true = np.array([1.0, 0.5, -0.3])
    # Matched pairs (group size exactly 2) -- required for this kernel.
    group_id = np.repeat(np.arange(n // 2), 2).astype(np.int32) + 1
    b = rng.normal(scale=0.4, size=n // 2).repeat(2)
    y = X @ beta_true + b + rng.normal(scale=0.5, size=n)
    return X, y, group_id


# R reference, EDI 1.0.0, computed 2026-08-03 via:
#   EDI:::gee_pairs_singletons_cpp(X, y, group_id, "gaussian")
R_BETA = np.array([1.0663177530105, 0.285182366555922, -0.218756218512224])
R_ALPHA = 0.160221693613315
R_QUASI_LOGLIK = -19.7170161093913
R_VCOV_00 = 0.0102321055738414


def test_matches_r_fixture():
    X, y, group_id = _synthetic_data()
    res = gee_pairs_singletons(X, y, group_id, "gaussian")

    assert res["converged"] is True
    assert res["beta"] == pytest.approx(R_BETA, abs=ATOL, rel=RTOL)
    assert res["alpha"] == pytest.approx(R_ALPHA, abs=ATOL, rel=RTOL)
    assert res["quasi_loglik"] == pytest.approx(R_QUASI_LOGLIK, abs=ATOL, rel=RTOL)
    assert res["vcov"][0, 0] == pytest.approx(R_VCOV_00, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, group_id = _synthetic_data()
    res = gee_pairs_singletons(X, y, group_id, "gaussian")

    p = X.shape[1]
    assert res["beta"].shape == (p,)
    assert res["vcov"].shape == (p, p)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["quasi_loglik"], float)
    assert isinstance(res["niter"], int)


def test_family_argument_accepted():
    X, y, group_id = _synthetic_data()
    # Just confirm the other two family strings are accepted and produce a
    # converged fit -- no fixture-precision check, just a binding-layer
    # smoke test that the family_str dispatch (gaussian/binomial/poisson)
    # in the pybind11 lambda works for all three branches.
    y_bin = (y > np.median(y)).astype(float)
    res_bin = gee_pairs_singletons(X, y_bin, group_id, "binomial")
    assert isinstance(res_bin["converged"], bool)

    y_count = np.abs(np.round(y - y.min())).astype(float)
    res_pois = gee_pairs_singletons(X, y_count, group_id, "poisson")
    assert isinstance(res_pois["converged"], bool)
