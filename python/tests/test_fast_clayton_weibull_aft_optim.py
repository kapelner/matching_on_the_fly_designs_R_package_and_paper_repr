"""Parity test for edi_kernels.fast_clayton_weibull_aft_optim against an
R-generated fixture.

This does NOT re-verify the statistical correctness of the Clayton-copula
Weibull AFT fit itself -- that's the same
fast_clayton_weibull_aft_optim_internal object code the R package's own
test suite already exercises via fast_clayton_weibull_aft_optim_cpp. What
this test covers, which nothing on the R side can: the pybind11 binding
layer itself (argument marshaling, defaults, edi::to_py_dict field
conversion) and the EDI_CORE_ONLY compiled path (this extension links
fetched Eigen + LBFGSpp, not the R build's RcppEigen/RcppNumerical).

pair_idx and singleton_rows are 0-based row indices into X/y/dead in BOTH
the Python binding and the R export (confirmed directly: the R fixture
below was computed with 0-based pair_idx/singleton_rows and matches the
Python side to float precision -- there is no off-by-one index-base
mismatch between the two languages for this kernel).

warm_start_params has no default on either side (length p + 2: beta, then
log_sigma, then log_theta -- see ClaytonWeibullLikelihood::operator() in
EDI/src/fast_survival_models_optim.cpp) -- a zero vector is used here.

The expected values below were computed once via:
    EDI:::fast_clayton_weibull_aft_optim_cpp(X, y, dead, pair_idx, singleton_rows, warm_start_params)
in R, on the exact same synthetic dataset generated below with
numpy.random.default_rng(104) (see package_metadata/python_bindings_package_spec.md
"Testing"). Do not regenerate this fixture casually -- if it needs
updating, regenerate from R and update the comment with the date/EDI version.
"""
import numpy as np
import pytest

from edi_kernels import fast_clayton_weibull_aft_optim

ATOL = 1e-9
RTOL = 1e-9


def _synthetic_data():
    rng = np.random.default_rng(104)

    n_pairs = 30
    n_single = 20
    n_tot = 2 * n_pairs + n_single
    p = 1
    X = rng.normal(size=(n_tot, p))
    beta_true = np.array([0.3])
    sigma_true = 0.7
    eps = rng.gumbel(size=n_tot) * sigma_true
    log_y = 0.5 + X[:, 0] * beta_true[0] + eps
    y = np.exp(log_y)
    dead = rng.binomial(1, 0.85, n_tot).astype(float)

    pair_idx = np.column_stack(
        [np.arange(0, 2 * n_pairs, 2), np.arange(1, 2 * n_pairs, 2)]
    ).astype(np.int32)
    singleton_rows = np.arange(2 * n_pairs, n_tot).astype(np.int32)
    warm_start_params = np.zeros(p + 2)

    return X, y, dead, pair_idx, singleton_rows, warm_start_params


# R reference, EDI 1.0.0, computed 2026-08-03 via:
#   EDI:::fast_clayton_weibull_aft_optim_cpp(X, y, dead, pair_idx, singleton_rows, warm_start_params)
# (pair_idx/singleton_rows passed 0-based, matching the Python/binding convention)
R_PARAMS = np.array([0.276881674797705, 0.5928472709728, 0.113338358330275])
R_NEG_LOGLIK = 207.431114196846


def test_matches_r_fixture():
    X, y, dead, pair_idx, singleton_rows, warm_start_params = _synthetic_data()
    res = fast_clayton_weibull_aft_optim(X, y, dead, pair_idx, singleton_rows, warm_start_params)

    assert res["converged"] is True
    assert res["params"] == pytest.approx(R_PARAMS, abs=ATOL, rel=RTOL)
    assert res["neg_loglik"] == pytest.approx(R_NEG_LOGLIK, abs=ATOL, rel=RTOL)


def test_result_shape_and_types():
    X, y, dead, pair_idx, singleton_rows, warm_start_params = _synthetic_data()
    res = fast_clayton_weibull_aft_optim(X, y, dead, pair_idx, singleton_rows, warm_start_params)

    assert res["params"].shape == (3,)  # beta (1) + log_sigma + log_theta
    assert isinstance(res["converged"], bool)
    assert isinstance(res["neg_loglik"], float)


def test_estimate_only_still_returns_params():
    X, y, dead, pair_idx, singleton_rows, warm_start_params = _synthetic_data()
    res = fast_clayton_weibull_aft_optim(
        X, y, dead, pair_idx, singleton_rows, warm_start_params, estimate_only=True
    )

    assert res["params"] == pytest.approx(R_PARAMS, abs=1e-6, rel=1e-6)
