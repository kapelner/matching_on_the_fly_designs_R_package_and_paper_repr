"""Parity test for edi_kernels.fast_coxph_regression against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Cox PH fit itself
-- that's the same fast_coxph_regression_internal object code the R
package's own (much larger) test suite already exercises via
fast_coxph_regression_cpp. What this test covers, which nothing on the R
side can: the pybind11 binding layer itself (argument marshaling, defaults,
edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled path (this
extension links fetched Eigen + LBFGSpp, not the R build's RcppEigen/
RcppNumerical).

The expected values below were computed once via:
    EDI:::fast_coxph_regression_cpp(X, y, dead)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(123) (see package_metadata/python_bindings_package_spec.md
"Testing"). Do not regenerate this fixture casually -- if it needs updating,
regenerate from R and update the comment with the date/EDI version.

Note: fast_coxph_regression's Python binding is deliberately unstratified
with no cluster-robust vcov (see python/cpp/bindings_survival.cpp); X has no
intercept column, matching the Cox partial-likelihood convention (the
baseline hazard absorbs it) -- this matches the R side's own test
convention in EDI/tests/testthat/test-rcpp-fitting-real-data.R.
"""
import numpy as np
import pytest

from edi_kernels import fast_coxph_regression

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(123)
    n = 150
    X = np.column_stack([
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    beta_true = np.array([0.5, -0.3])
    eta = X @ beta_true
    scale = np.exp(-eta)
    y = rng.exponential(scale)
    censor_time = rng.exponential(3.0, n)
    dead = (y <= censor_time).astype(float)
    y_obs = np.minimum(y, censor_time)
    return X, y_obs, dead


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_coxph_regression_cpp(X, y, dead)
R_COEFFICIENTS = np.array([0.654770101904559, -0.319961675829669])
R_NEG_LL = 492.766701831856
R_VCOV_DIAG = np.array([0.0346329954487925, 0.0091551597979166])
R_ITERATIONS = 4


def test_matches_r_fixture():
    X, y, dead = _synthetic_data()
    res = fast_coxph_regression(X, y, dead)

    assert res["converged"] is True
    assert res["coefficients"] == pytest.approx(R_COEFFICIENTS, abs=ATOL, rel=RTOL)
    assert res["neg_ll"] == pytest.approx(R_NEG_LL, abs=ATOL, rel=RTOL)
    assert np.diag(res["vcov"]) == pytest.approx(R_VCOV_DIAG, abs=ATOL, rel=RTOL)
    assert res["iterations"] == R_ITERATIONS


def test_result_shape_and_types():
    X, y, dead = _synthetic_data()
    res = fast_coxph_regression(X, y, dead)

    p = X.shape[1]
    assert res["coefficients"].shape == (p,)
    assert res["vcov"].shape == (p, p)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_ll"], float)


def test_estimate_only_omits_vcov_fields():
    X, y, dead = _synthetic_data()
    res = fast_coxph_regression(X, y, dead, estimate_only=True)

    # estimate_only skips the vcov computation entirely (the key is absent,
    # not set to None) but still includes fisher_information/gradient_norm
    # -- see fast_coxph_regression_internal in EDI/src/fast_coxph_regression.cpp.
    assert "vcov" not in res
    assert res["fisher_information"].shape == (X.shape[1], X.shape[1])
    assert res["coefficients"] == pytest.approx(R_COEFFICIENTS, abs=1e-6, rel=1e-6)
