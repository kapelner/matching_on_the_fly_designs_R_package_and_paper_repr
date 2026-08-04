"""Parity test for edi_kernels.fast_weibull_regression against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Weibull AFT fit
itself -- that's the same fast_weibull_regression_internal object code the R
package's own test suite already exercises via fast_weibull_regression_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

Same synthetic-data recipe as test_fast_coxph_regression.py (same seed,
same X/y/dead construction) so both survival kernels are exercised on
directly comparable data.

The expected values below were computed once via:
    EDI:::fast_weibull_regression_cpp(X, y, dead)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(123). Do not regenerate this fixture casually --
if it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_weibull_regression

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
#   EDI:::fast_weibull_regression_cpp(X, y, dead)
# params = [beta_0, beta_1, log_sigma]
R_PARAMS = np.array([-0.580572508466287, 0.35362562665742, 0.115502572892534])
R_NEG_LOGLIK = 100.088153200844
R_VCOV_DIAG = np.array([0.0210927071459095, 0.0105260329864121, 0.00435393717820719])


def test_matches_r_fixture():
    X, y, dead = _synthetic_data()
    res = fast_weibull_regression(X, y, dead)

    assert res["converged"] is True
    assert res["params"] == pytest.approx(R_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)
    assert np.diag(res["vcov"]) == pytest.approx(R_VCOV_DIAG, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, dead = _synthetic_data()
    res = fast_weibull_regression(X, y, dead)

    p = X.shape[1]
    total = p + 1  # beta + log_sigma
    assert res["params"].shape == (total,)
    assert res["vcov"].shape == (total, total)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_still_fits_params():
    X, y, dead = _synthetic_data()
    res = fast_weibull_regression(X, y, dead, estimate_only=True)

    # estimate_only returns beta/log_sigma split into separate "b"/"log_sigma"
    # keys rather than the combined "params" vector the full fit returns --
    # and omits "vcov"/"neg_loglik" entirely rather than setting them to None.
    assert "vcov" not in res
    assert "params" not in res
    p = X.shape[1]
    assert res["b"] == pytest.approx(R_PARAMS[:p], abs=1e-6, rel=1e-6)
    assert res["log_sigma"] == pytest.approx(R_PARAMS[p], abs=1e-6, rel=1e-6)
