"""Parity test for edi_kernels.fast_clogit_plus_glmm against an R-generated fixture.

This does NOT re-verify the statistical correctness of the combined
conditional-logit (discordant pairs) + random-intercept-GLMM (concordant
pairs) fit itself -- that's the same fast_clogit_plus_glmm_internal object
code the R package's own test suite already exercises via
fast_clogit_plus_glmm_cpp. What this test covers, which nothing on the R
side can: the pybind11 binding layer itself (argument marshaling, defaults,
edi::to_py_dict field conversion) and the EDI_CORE_ONLY compiled path (this
extension links fetched Eigen + LBFGSpp, not the R build's
RcppEigen/RcppNumerical).

The concordant-pair GLMM half uses the same fixed-size-2 (SM2/SV2) Eigen
analytical shortcut as the rest of the GLMM/CLMM/LMM family -- group_conc
MUST be matched pairs (group size exactly 2); a larger group silently
corrupts memory (see this session's fast_gaussian_lmm investigation).

The expected values below were computed once via:
    EDI:::fast_clogit_plus_glmm_cpp(X_disc, y_disc, X_conc, y_conc, group_conc, TRUE, TRUE)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(103) (see package_metadata/python_bindings_package_spec.md
"Testing"). Do not regenerate this fixture casually -- if it needs
updating, regenerate from R and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_clogit_plus_glmm

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(103)

    # Discordant pairs -> conditional logit part.
    n_disc = 80
    p = 2
    X_disc = rng.normal(size=(n_disc, p))
    beta_true = np.array([0.6, -0.4])
    eta_d = X_disc @ beta_true
    y_disc = rng.binomial(1, 1 / (1 + np.exp(-eta_d)), n_disc).astype(float)

    # Concordant pairs -> random-intercept GLMM part, matched pairs (size 2 ONLY).
    n_conc = 80
    X_conc = rng.normal(size=(n_conc, p))
    group_conc = np.repeat(np.arange(n_conc // 2), 2).astype(np.int32) + 1
    b_conc = rng.normal(scale=0.5, size=n_conc // 2).repeat(2)
    eta_c = X_conc @ beta_true + b_conc
    y_conc = rng.binomial(1, 1 / (1 + np.exp(-eta_c)), n_conc).astype(float)

    return X_disc, y_disc, X_conc, y_conc, group_conc


# R reference, EDI 1.0.0, computed 2026-08-03 via:
#   EDI:::fast_clogit_plus_glmm_cpp(X_disc, y_disc, X_conc, y_conc, group_conc, TRUE, TRUE)
R_B = np.array([0.757061395065093, 0.141544204860354, -0.353967177233617])
R_NEG_LOGLIK = 105.177227452469


def test_matches_r_fixture():
    X_disc, y_disc, X_conc, y_conc, group_conc = _synthetic_data()
    res = fast_clogit_plus_glmm(X_disc, y_disc, X_conc, y_conc, group_conc, True, True)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X_disc, y_disc, X_conc, y_conc, group_conc = _synthetic_data()
    res = fast_clogit_plus_glmm(X_disc, y_disc, X_conc, y_conc, group_conc, True, True)

    assert res["b"].shape == (3,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_still_returns_b():
    X_disc, y_disc, X_conc, y_conc, group_conc = _synthetic_data()
    res = fast_clogit_plus_glmm(
        X_disc, y_disc, X_conc, y_conc, group_conc, True, True, estimate_only=True
    )

    assert res["b"] == pytest.approx(R_B, abs=1e-6, rel=1e-6)
