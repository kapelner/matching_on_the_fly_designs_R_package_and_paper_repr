"""Parity test for edi_kernels.fast_poisson_glmm against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Poisson GLMM fit
itself -- that's the same fast_poisson_glmm_internal object code the R
package's own (much larger) test suite already exercises via
fast_poisson_glmm_cpp. What this test covers, which nothing on the R side
can: the pybind11 binding layer itself (argument marshaling, defaults,
edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled path (this
extension links fetched Eigen + LBFGSpp, not the R build's RcppEigen/
RcppNumerical).

The expected values below were computed once via:
    EDI:::fast_poisson_glmm_cpp(X, y, group_id, j_T = 1L)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(42) (see package_metadata/python_bindings_package_spec.md
"Testing"). Do not regenerate this fixture casually -- if it needs updating,
regenerate from R and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_poisson_glmm

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(42)
    n = 200
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n),
    ])
    group_id = np.repeat(np.arange(n // 4), 4).astype(np.int32)
    beta_true = np.array([0.5, 0.3, -0.2])
    eta = X @ beta_true
    mu = np.exp(eta + rng.normal(scale=0.3, size=n // 4).repeat(4))
    y = rng.poisson(mu).astype(float)
    return X, y, group_id


# R reference, EDI 1.0.0, computed 2026-08-02 via:
#   EDI:::fast_poisson_glmm_cpp(X, y, group_id, j_T = 1L)
R_B = np.array([0.396523267557659, 0.325077586190454, -0.181464003915205])
R_LOG_SIGMA = -0.987593452504449
R_NEG_LOGLIK = 344.026116874617912
R_VCOV_00 = 0.010051111518828
R_VCOV_11 = 0.012777563541711


def test_matches_r_fixture():
    X, y, group_id = _synthetic_data()
    res = fast_poisson_glmm(X, y, group_id, j_T=1)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["log_sigma"] == pytest.approx(R_LOG_SIGMA, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)
    assert res["vcov"][0, 0] == pytest.approx(R_VCOV_00, abs=ATOL, rel=RTOL)
    assert res["vcov"][1, 1] == pytest.approx(R_VCOV_11, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, group_id = _synthetic_data()
    res = fast_poisson_glmm(X, y, group_id, j_T=1)

    p = X.shape[1]
    total = p + 1
    assert res["b"].shape == (p,)
    assert res["vcov"].shape == (total, total)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_omits_vcov_fields():
    X, y, group_id = _synthetic_data()
    res = fast_poisson_glmm(X, y, group_id, j_T=1, estimate_only=True)

    assert res["vcov"] is None
    assert res["fisher_information"] is None
    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
