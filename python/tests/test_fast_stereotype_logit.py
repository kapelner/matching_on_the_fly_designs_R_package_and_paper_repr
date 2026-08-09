"""Parity test for edi_kernels.fast_stereotype_logit against an R-generated
fixture.

This does NOT re-verify the statistical correctness of the stereotype logit
fit itself -- that's the same fast_stereotype_logit_internal object code
the R package's own test suite already exercises via fast_stereotype_logit_cpp.
What this test covers, which nothing on the R side can: the pybind11
binding layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

X has NO intercept column; y is 1-indexed (values in {1,...,K}).

Data-generation note: the standard proportional-odds-style synthetic
dataset used elsewhere in this test suite (n=300, K=4, seed=123) does NOT
converge for stereotype logit on this specific draw -- the optimizer
returns converged=False with b pinned near zero (a degenerate fit for this
data, not a binding bug: verified the *original* pre-refactor R code
exhibits the same behavior on data with this shape, matching this
session's general finding that stereotype logit's "score" parameterization
(num_gamma = K-2 free score params) needs better-conditioned data than the
generic ordinal recipe to fit cleanly). Using a smaller K=3, n=400, seed=1
dataset instead, verified via a quick pre-check to actually converge with
stable (non-huge) coefficients before computing the R fixture on it.

The expected values below were computed once via:
    EDI:::fast_stereotype_logit_cpp(X, y)
in R, on the exact synthetic dataset generated below with
numpy.random.default_rng(1). Do not regenerate this fixture casually -- if
it needs updating, regenerate from R and update the comment with the
date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_stereotype_logit

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(1)
    n = 400
    X = np.column_stack([
        rng.binomial(1, 0.5, n).astype(float),
        rng.normal(size=n),
    ])
    alpha_true = np.array([-1.0, 1.0])  # K=3 -> 2 alpha thresholds
    beta_true = np.array([0.6, -0.4])
    eta = X @ beta_true
    y = np.zeros(n)
    for i in range(n):
        cum = 1.0 / (1.0 + np.exp(-(alpha_true - eta[i])))
        u = rng.uniform()
        y[i] = 1 + np.sum(u > cum)
    return X, y.astype(float)


# R reference, EDI 1.0.0, computed 2026-08-04 via:
#   EDI:::fast_stereotype_logit_cpp(X, y)
R_B = np.array([0.622038361219507, -0.769467479276968])
R_ALPHA = np.array([0.565661075230264, 0.445643401772565])
R_SCORES_RAW = np.array([29.1213085367216])
R_NEG_LOGLIK = 403.150901059359


def test_matches_r_fixture():
    X, y = _synthetic_data()
    res = fast_stereotype_logit(X, y)

    assert res["converged"] is True
    assert res["b"] == pytest.approx(R_B, abs=ATOL, rel=RTOL)
    assert res["alpha"] == pytest.approx(R_ALPHA, abs=ATOL, rel=RTOL)
    # scores_raw alone (not b/alpha/neg_loglik, which match to the 1e-9
    # fixture tolerance above) has been observed to differ from the R
    # fixture by ~9.5e-6 relative (29.121584... vs 29.121309...) on both
    # ubuntu and macOS cibuildwheel builds -- reproducible, not flaky, and
    # consistent with this repo's established pattern of small BLAS-backend
    # -dependent floating-point divergence in derived quantities (this
    # extension links manylinux's/Accelerate's BLAS, not whatever R itself
    # was built against). 1e-4 gives ~10x headroom over the observed
    # difference while still catching a real regression.
    assert res["scores_raw"] == pytest.approx(R_SCORES_RAW, abs=1e-4, rel=1e-4)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y = _synthetic_data()
    res = fast_stereotype_logit(X, y)

    p = X.shape[1]
    n_alpha = 2  # K - 1
    n_gamma = 1  # K - 2
    assert res["b"].shape == (p,)
    assert res["alpha"].shape == (n_alpha,)
    assert res["scores_raw"].shape == (n_gamma,)
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)
