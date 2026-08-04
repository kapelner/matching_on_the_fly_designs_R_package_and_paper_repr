"""Parity test for edi_kernels.fast_zero_augmented_poisson against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the zero-augmented
Poisson fit itself -- that's the same fast_zap_internal object code the R
package's own test suite already exercises via fast_zero_augmented_poisson_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

is_hurdle picks between two genuinely different models sharing one kernel:
False = zero-inflated Poisson, True = hurdle Poisson -- both variants are
covered here on the same X/y/Xzi data. Xzi is the zero-inflation/hurdle
design matrix. Returned 'params' packs [beta_cond, beta_zi].

The expected values below were computed once via:
    EDI:::fast_zero_augmented_poisson_cpp(X, y, Xzi, is_hurdle = FALSE/TRUE)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(71). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_zero_augmented_poisson

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(71)
    n = 300
    X = np.column_stack([
        np.ones(n),
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    Xzi = X.copy()
    mu = np.exp(0.6 + 0.3 * X[:, 1] - 0.2 * X[:, 2])
    y_count = rng.poisson(mu).astype(float)
    infl_p = 1.0 / (1.0 + np.exp(-(0.1 - 0.4 * Xzi[:, 1])))
    zero_infl = rng.binomial(1, infl_p, n)
    y = np.where(zero_infl == 1, 0.0, y_count)
    return X, y, Xzi


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_zero_augmented_poisson_cpp(X, y, Xzi, is_hurdle = FALSE/TRUE)
R_ZIP_PARAMS = np.array([
    0.454031351041028, 0.380676748540686, -0.348763233239726,
    -0.117788283549738, -0.325753378018191, -0.0712384642478158,
])
R_ZIP_NEG_LOGLIK = 207.744362282367

R_HURDLE_PARAMS = np.array([
    0.458817047452881, 0.371801420635466, -0.351153719656712,
    0.361886124223062, -0.519088193215444, 0.155301714950873,
])
R_HURDLE_NEG_LOGLIK = 207.81619218983


def test_matches_r_fixture_zero_inflated():
    X, y, Xzi = _synthetic_data()
    res = fast_zero_augmented_poisson(X, y, Xzi, is_hurdle=False)

    assert res["converged"] is True
    assert res["params"] == pytest.approx(R_ZIP_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_ZIP_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_matches_r_fixture_hurdle():
    X, y, Xzi = _synthetic_data()
    res = fast_zero_augmented_poisson(X, y, Xzi, is_hurdle=True)

    assert res["converged"] is True
    assert res["params"] == pytest.approx(R_HURDLE_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_HURDLE_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, Xzi = _synthetic_data()
    res = fast_zero_augmented_poisson(X, y, Xzi, is_hurdle=False)

    assert res["params"].shape == (6,)  # 3 + 3
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)
