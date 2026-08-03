"""Parity test for edi_kernels.fast_ordinal_regression against an R-generated fixture.

This does NOT re-verify the statistical correctness of the proportional-odds
fit itself -- that's the same fast_ordinal_regression_internal object code
the R package's own test suite already exercises via fast_ordinal_regression_cpp.
What this test covers, which nothing on the R side can: the pybind11 binding
layer (argument marshaling, defaults, edi::to_py_dict field conversion) and
the EDI_CORE_ONLY compiled path (this extension links fetched Eigen +
LBFGSpp, not the R build's RcppEigen/RcppNumerical).

X has NO intercept column -- the proportional-odds model's K-1 alpha
thresholds already serve that role; an intercept column would be collinear
with them.

The expected values below were computed once via:
    EDI:::fast_ordinal_regression_cpp(X, y)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(123) (see package_metadata/python_bindings_package_spec.md
"Testing"). Do not regenerate this fixture casually -- if it needs updating,
regenerate from R and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_ordinal_regression

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
#   EDI:::fast_ordinal_regression_cpp(X, y)
R_B = np.array([0.746048853778645, -0.403753252367445])
R_ALPHA = np.array([-1.00692105860466, 0.368033067942919, 1.70976960260878])
R_NEG_LOGLIK = 400.548667410759
R_VCOV_00 = 0.0260708285429491
R_VCOV_11 = 0.0225646931130649


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_ordinal_regression(X, y)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["alpha"] == pytest.approx(R_ALPHA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)
    assert res["vcov"][0, 0] == pytest.approx(R_VCOV_00, abs=ATOL, rel=RTOL)
    assert res["vcov"][1, 1] == pytest.approx(R_VCOV_11, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_ordinal_regression(X, y)

    p = X.shape[1]
    n_alpha = 3
    assert res["b"].shape == (p,)
    assert res["alpha"].shape == (n_alpha,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_omits_vcov():
    X, y = _synthetic_data()
    res = fast_ordinal_regression(X, y, estimate_only=True)

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
    assert "vcov" not in res or res.get("vcov") is None
