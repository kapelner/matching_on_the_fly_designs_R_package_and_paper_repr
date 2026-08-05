"""Systematic omitted-argument coverage for every bound EDI kernel.

Acceptance Criterion 6 (see python_bindings_package_spec.md) requires that
every `Rcpp::Nullable<T>` -> `std::optional<T>` argument have "a working
`std::optional`-based Python equivalent with a passing omitted-argument
test". The per-kernel parity test files (test_fast_*.py) already call most
kernels with few or no optional kwargs incidentally, as part of matching an
R fixture -- but that's incidental coverage of *some* defaults, not a
deliberate, systematic check that *every* optional argument on *every*
bound kernel has a working default.

This file is that deliberate check. For each of the 37 model-fitting
kernels (all bound functions except fast_pchisq_upper, which isn't a
Nullable-bearing model-fitting kernel), it calls the function with *only*
the required (no-default) positional/keyword arguments -- i.e. every
`std::optional<T>` parameter is omitted, not just left at an
explicitly-passed default -- and checks the call succeeds and returns a
sane, finite result. If any binding's default construction were ever
broken (e.g. a `std::optional` that isn't actually optional, or a default
that produces NaN/crashes), this is the test that would catch it.

This does NOT re-verify statistical correctness against R -- that's what
the parity tests (test_fast_*.py's test_matches_r_fixture) are for. Data
generation here is deliberately reused from each kernel's own parity test
module (via a plain sibling import -- pytest's default rootdir-relative
sys.path insertion makes this work without a package __init__.py) so the
two test files can't silently drift apart on what "valid input" means for
a given kernel.
"""
import numpy as np

from edi_kernels import (
    fast_ols,
    fast_robust_regression,
    fast_logistic_regression,
    fast_probit_regression,
    fast_log_binomial_regression,
    fast_log_binomial_regression_with_var,
    fast_identity_binomial_regression,
    fast_identity_binomial_regression_with_var,
    fast_poisson_regression,
    fast_neg_bin,
    fast_zinb,
    fast_zero_augmented_poisson,
    fast_hurdle_negbin,
    fast_truncated_negbin_count,
    fast_cpoisson_combined,
    fast_beta_regression,
    fast_zero_one_inflated_beta,
    fast_adjacent_category_logit,
    fast_continuation_ratio_regression,
    fast_ordinal_regression,
    fast_ordinal_probit_regression,
    fast_ordinal_cauchit_regression,
    fast_ordinal_cloglog_regression,
    fast_ordinal_clmm,
    fast_ordinal_glmm,
    fast_stereotype_logit,
    gee_pairs_singletons,
    fast_coxph_regression,
    fast_weibull_regression,
    fast_weibull_frailty,
    fast_clayton_weibull_aft_optim,
    fast_dep_cens_transform_optim,
    fast_poisson_glmm,
    fast_gaussian_lmm,
    fast_logistic_glmm,
    fast_hurdle_poisson_glmm,
    fast_clogit_plus_glmm,
)

from test_fast_ols import _synthetic_data as _ols_data
from test_fast_robust_regression import _synthetic_data as _robust_data
from test_fast_logistic_regression import _synthetic_data as _logistic_data
from test_fast_probit_regression import _synthetic_data as _probit_data
from test_fast_log_binomial_regression import _synthetic_data as _log_binom_data
from test_fast_log_binomial_regression_with_var import _synthetic_data as _log_binom_var_data
from test_fast_identity_binomial_regression import _synthetic_data as _identity_binom_data
from test_fast_identity_binomial_regression_with_var import _synthetic_data as _identity_binom_var_data
from test_fast_poisson_regression import _synthetic_data as _poisson_data
from test_fast_neg_bin import _synthetic_data as _neg_bin_data
from test_fast_zinb import _synthetic_data as _zinb_data
from test_fast_zero_augmented_poisson import _synthetic_data as _zap_data
from test_fast_hurdle_negbin import _synthetic_data as _hurdle_negbin_data
from test_fast_truncated_negbin_count import _synthetic_data as _trunc_negbin_data
from test_fast_cpoisson_combined import _synthetic_data as _cpoisson_data
from test_fast_beta_regression import _synthetic_data as _beta_data
from test_fast_zero_one_inflated_beta import _synthetic_data as _zoib_data
from test_fast_adjacent_category_logit import _synthetic_data as _adj_cat_data
from test_fast_continuation_ratio_regression import _synthetic_data as _cont_ratio_data
from test_fast_ordinal_regression import _synthetic_data as _ordinal_data
from test_fast_ordinal_probit_regression import _synthetic_data as _ordinal_probit_data
from test_fast_ordinal_cauchit_regression import _synthetic_data as _ordinal_cauchit_data
from test_fast_ordinal_cloglog_regression import _synthetic_data as _ordinal_cloglog_data
from test_fast_ordinal_clmm import _synthetic_data as _ordinal_clmm_data
from test_fast_ordinal_glmm import _synthetic_data as _ordinal_glmm_data
from test_fast_stereotype_logit import _synthetic_data as _stereotype_data
from test_gee_pairs_singletons import _synthetic_data as _gee_data
from test_fast_coxph_regression import _synthetic_data as _coxph_data
from test_fast_weibull_regression import _synthetic_data as _weibull_data
from test_fast_weibull_frailty import _synthetic_data as _weibull_frailty_data
from test_fast_clayton_weibull_aft_optim import _synthetic_data as _clayton_data
from test_fast_dep_cens_transform_optim import _synthetic_data as _dep_cens_data
from test_fast_poisson_glmm import _synthetic_data as _poisson_glmm_data
from test_fast_gaussian_lmm import _synthetic_data as _gaussian_lmm_data
from test_fast_logistic_glmm import _synthetic_data as _logistic_glmm_data
from test_fast_hurdle_poisson_glmm import _synthetic_data as _hurdle_poisson_glmm_data
from test_fast_clogit_plus_glmm import _synthetic_data as _clogit_plus_glmm_data


def _assert_valid_result(res):
    """Every optional arg was omitted; the call must still succeed and
    produce a well-formed, finite result using the binding's own C++
    defaults (not sentinel NaNs, not a crash)."""
    assert isinstance(res, dict)
    assert len(res) > 0
    for key, value in res.items():
        if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.floating):
            assert np.all(np.isfinite(value)), f"non-finite values in '{key}'"


def test_fast_ols_omitted():
    X, y = _ols_data()
    _assert_valid_result(fast_ols(X, y))


def test_fast_robust_regression_omitted():
    X, y = _robust_data()
    _assert_valid_result(fast_robust_regression(X, y))


def test_fast_logistic_regression_omitted():
    X, y = _logistic_data()
    _assert_valid_result(fast_logistic_regression(X, y))


def test_fast_probit_regression_omitted():
    X, y = _probit_data()
    _assert_valid_result(fast_probit_regression(X, y))


def test_fast_log_binomial_regression_omitted():
    X, y = _log_binom_data()
    _assert_valid_result(fast_log_binomial_regression(X, y))


def test_fast_log_binomial_regression_with_var_omitted():
    X, y = _log_binom_var_data()
    _assert_valid_result(fast_log_binomial_regression_with_var(X, y))


def test_fast_identity_binomial_regression_omitted():
    X, y = _identity_binom_data()
    _assert_valid_result(fast_identity_binomial_regression(X, y))


def test_fast_identity_binomial_regression_with_var_omitted():
    X, y = _identity_binom_var_data()
    _assert_valid_result(fast_identity_binomial_regression_with_var(X, y))


def test_fast_poisson_regression_omitted():
    X, y = _poisson_data()
    _assert_valid_result(fast_poisson_regression(X, y))


def test_fast_neg_bin_omitted():
    X, y = _neg_bin_data()
    _assert_valid_result(fast_neg_bin(X, y))


def test_fast_zinb_omitted():
    Xc, Xz, y = _zinb_data()
    _assert_valid_result(fast_zinb(Xc, Xz, y))


def test_fast_zero_augmented_poisson_omitted():
    X, y, Xzi = _zap_data()
    # is_hurdle has no default in the binding -- it's a required arg, not
    # an optional one, so it's supplied here alongside X/y/Xzi; everything
    # after it in the signature is omitted.
    _assert_valid_result(fast_zero_augmented_poisson(X, y, Xzi, False))


def test_fast_hurdle_negbin_omitted():
    X, y, X_hurdle = _hurdle_negbin_data()
    _assert_valid_result(fast_hurdle_negbin(X, y, X_hurdle))


def test_fast_truncated_negbin_count_omitted():
    X, y = _trunc_negbin_data()
    _assert_valid_result(fast_truncated_negbin_count(X, y))


def test_fast_cpoisson_combined_omitted():
    yT_v, n_k_v, X_diff_v, y_r, w_r, X_r = _cpoisson_data()
    _assert_valid_result(fast_cpoisson_combined(yT_v, n_k_v, X_diff_v, y_r, w_r, X_r))


def test_fast_beta_regression_omitted():
    X, y = _beta_data()
    _assert_valid_result(fast_beta_regression(X, y))


def test_fast_zero_one_inflated_beta_omitted():
    X, X_zero_one, y = _zoib_data()
    _assert_valid_result(fast_zero_one_inflated_beta(X, X_zero_one, y))


def test_fast_adjacent_category_logit_omitted():
    X, y = _adj_cat_data()
    _assert_valid_result(fast_adjacent_category_logit(X, y))


def test_fast_continuation_ratio_regression_omitted():
    X, y = _cont_ratio_data()
    _assert_valid_result(fast_continuation_ratio_regression(X, y))


def test_fast_ordinal_regression_omitted():
    X, y = _ordinal_data()
    _assert_valid_result(fast_ordinal_regression(X, y))


def test_fast_ordinal_probit_regression_omitted():
    X, y = _ordinal_probit_data()
    _assert_valid_result(fast_ordinal_probit_regression(X, y))


def test_fast_ordinal_cauchit_regression_omitted():
    X, y = _ordinal_cauchit_data()
    _assert_valid_result(fast_ordinal_cauchit_regression(X, y))


def test_fast_ordinal_cloglog_regression_omitted():
    X, y = _ordinal_cloglog_data()
    _assert_valid_result(fast_ordinal_cloglog_regression(X, y))


def test_fast_ordinal_clmm_omitted():
    X, y, group_id = _ordinal_clmm_data()
    # K/j_T have no defaults (required); link="logit" IS defaulted and is
    # deliberately omitted here.
    _assert_valid_result(fast_ordinal_clmm(X, y, group_id, K=3, j_T=0))


def test_fast_ordinal_glmm_omitted():
    X, y, group_id = _ordinal_glmm_data()
    _assert_valid_result(fast_ordinal_glmm(X, y, group_id, K=3, j_T=0))


def test_fast_stereotype_logit_omitted():
    X, y = _stereotype_data()
    _assert_valid_result(fast_stereotype_logit(X, y))


def test_gee_pairs_singletons_omitted():
    X, y, group_id = _gee_data()
    # family has no default (required).
    _assert_valid_result(gee_pairs_singletons(X, y, group_id, "gaussian"))


def test_fast_coxph_regression_omitted():
    X, y_obs, dead = _coxph_data()
    _assert_valid_result(fast_coxph_regression(X, y_obs, dead))


def test_fast_weibull_regression_omitted():
    X, y_obs, dead = _weibull_data()
    _assert_valid_result(fast_weibull_regression(X, y_obs, dead))


def test_fast_weibull_frailty_omitted():
    X, y_obs, dead, group_id = _weibull_frailty_data()
    _assert_valid_result(fast_weibull_frailty(X, y_obs, dead, group_id))


def test_fast_clayton_weibull_aft_optim_omitted():
    X, y, dead, pair_idx, singleton_rows, warm_start_params = _clayton_data()
    # warm_start_params has no default in this binding (required) -- every
    # arg after it (estimate_only, maxit, reltol, fixed_idx/values,
    # optimization_alg, warm_start_fisher_info) is what's actually omitted.
    _assert_valid_result(
        fast_clayton_weibull_aft_optim(X, y, dead, pair_idx, singleton_rows, warm_start_params)
    )


def test_fast_dep_cens_transform_optim_omitted():
    X, y, dead = _dep_cens_data()
    _assert_valid_result(fast_dep_cens_transform_optim(X, y, dead))


def test_fast_poisson_glmm_omitted():
    X, y, group_id = _poisson_glmm_data()
    # j_T has no default (required).
    _assert_valid_result(fast_poisson_glmm(X, y, group_id, j_T=1))


def test_fast_gaussian_lmm_omitted():
    X, y, group_id = _gaussian_lmm_data()
    _assert_valid_result(fast_gaussian_lmm(X, y, group_id))


def test_fast_logistic_glmm_omitted():
    X, y, group_id = _logistic_glmm_data()
    _assert_valid_result(fast_logistic_glmm(X, y, group_id, j_T=1))


def test_fast_hurdle_poisson_glmm_omitted():
    X, y, group_id = _hurdle_poisson_glmm_data()
    _assert_valid_result(fast_hurdle_poisson_glmm(X, y, group_id, j_T=1))


def test_fast_clogit_plus_glmm_omitted():
    X_disc, y_disc, X_conc, y_conc, group_conc = _clogit_plus_glmm_data()
    # has_discordant/has_concordant have no defaults (required).
    _assert_valid_result(
        fast_clogit_plus_glmm(X_disc, y_disc, X_conc, y_conc, group_conc, True, True)
    )
