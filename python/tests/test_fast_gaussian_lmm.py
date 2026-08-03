"""Parity test for edi_kernels.fast_gaussian_lmm against an R-generated fixture.

This does NOT re-verify the statistical correctness of the Gaussian LMM fit
itself -- that's the same fast_gaussian_lmm_internal object code the R
package's own (much larger) test suite already exercises via
fast_gaussian_lmm_cpp. What this test covers, which nothing on the R side
can: the pybind11 binding layer itself (argument marshaling, defaults,
edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled path (this
extension links fetched Eigen + LBFGSpp, not the R build's RcppEigen/
RcppNumerical).

fast_gaussian_lmm (like the rest of the GLMM/CLMM/LMM family) uses an
internal fixed-size-2 Eigen type as an analytical shortcut for the
random-intercept variance structure -- group_id MUST partition observations
into groups of size exactly 1 or 2 (matched pairs + singleton reservoir
subjects), never larger, or the fit corrupts memory.

The expected values below were computed once via:
    EDI:::fast_gaussian_lmm_cpp(X, y, group_id)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(101). Do not regenerate this fixture casually --
if it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_gaussian_lmm

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(101)
    n = 200
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n),
        rng.normal(size=n),
    ])
    group_id = np.repeat(np.arange(n // 2), 2).astype(np.int32) + 1
    beta_true = np.array([1.0, 0.5, -0.3])
    b_re = rng.normal(scale=0.4, size=n // 2).repeat(2)
    y = X @ beta_true + b_re + rng.normal(scale=0.5, size=n)
    return X, y, group_id


# R reference, EDI 1.0.0, computed 2026-08-03 via:
#   EDI:::fast_gaussian_lmm_cpp(X, y, group_id)
R_B = np.array([
    1.00773196170803, 0.445034732202075, -0.27302794853661,
    -0.736898694315722, -1.00781610999521,
])
R_NITER = 8
R_SSQ_B_T = 0.00623585036175557
R_NEG_LOGLIK = 174.999915889948
R_VCOV_00 = 0.0038835084661506
R_VCOV_11 = 0.00623585036175557


def test_matches_r_fixture():
    X, y, group_id = _synthetic_data()
    res = fast_gaussian_lmm(X, y, group_id)

    assert res["converged"] is True
    assert res["niter"] == R_NITER
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["ssq_b_T"] == pytest.approx(R_SSQ_B_T, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)
    assert res["vcov"][0, 0] == pytest.approx(R_VCOV_00, abs=ATOL, rel=RTOL)
    assert res["vcov"][1, 1] == pytest.approx(R_VCOV_11, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, group_id = _synthetic_data()
    res = fast_gaussian_lmm(X, y, group_id)

    p = X.shape[1]
    total = p + 2  # beta + log_sigma_e + log_sigma_b
    assert res["b"].shape == (total,)
    assert res["vcov"].shape == (total, total)
    assert res["fisher_information"].shape == (total, total)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_matches_full_fit_beta():
    X, y, group_id = _synthetic_data()
    res = fast_gaussian_lmm(X, y, group_id, estimate_only=True)

    # Unlike fast_poisson_glmm (which includes vcov/fisher_information as
    # None in its estimate_only branch), fast_gaussian_lmm's estimate_only
    # branch omits those keys entirely -- see fast_gaussian_lmm_internal in
    # EDI/src/fast_gaussian_lmm.cpp. Documenting the actual contract here
    # rather than assuming parity with the poisson_glmm kernel.
    assert "vcov" not in res
    assert "fisher_information" not in res
    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)


def test_invalid_group_size_raises_or_is_avoided_by_construction():
    """Documents the matched-pairs-only constraint this whole family shares
    (see module docstring) -- group sizes > 2 are invalid input, not
    something this kernel is expected to validate defensively. This test
    exists so the constraint is visible in the test file, not just in a
    comment; it does not call the kernel with bad data (that risks the
    exact memory corruption the docstring warns about)."""
    X, y, group_id = _synthetic_data()
    _, counts = np.unique(group_id, return_counts=True)
    assert set(counts.tolist()) <= {1, 2}
