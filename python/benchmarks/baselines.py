"""Canonical Python baseline registry for the EDI model-fitting kernels.

This is the "families with a confirmed clean baseline" subset of the
Python Baseline Registry table in
``package_metadata/python_bindings_package_spec.md`` (TODO-5): OLS, robust
regression, logistic, probit, log-binomial (+ identity-binomial), Poisson,
NegBin, ZINB/ZAP/hurdle-NegBin, beta regression, Cox PH, Weibull AFT,
proportional-odds/cauchit/cloglog ordinal, and GEE.

Every entry's ``fit`` callable takes already-constructed NumPy arrays (no
pandas/DataFrame overhead in the timed region -- see the "Apples-to-apples
canonical timing" rule in the spec's Benchmark Harness section) and returns
the fitted result object, so a benchmark runner can time ``fit(...)``
directly. This module intentionally does NOT do dataset generation, EDI-side
timing, or report rendering -- that's ``run_benchmark_audit.py``'s job (see
Package Layout in the spec doc); this module is only the baseline registry.

Findings this registry encodes (see python_bindings_package_spec.md's TODO-4
status note, verified 2026-08-03 against statsmodels==0.14.6,
scikit-survival==0.28.0, lifelines==0.30.0):

- ``OrderedModel`` supports cauchit/cloglog links, but only via a
  ``scipy.stats`` distribution *instance* passed as ``distr=`` -- the
  ``distr=`` string shortcut only accepts ``'logit'``/``'probit'``.
  ``scipy.stats.cauchy`` is the cauchit link; ``scipy.stats.gumbel_l``'s CDF
  is exactly the cloglog inverse link (``1 - exp(-exp(x))``).
- EDI's ``fast_robust_regression`` defaults to ``method="MM"`` (a two-stage
  S-then-M estimator); ``statsmodels.robust.robust_linear_model.RLM``'s own
  default resolves to plain ``HuberT`` M-estimation. These are different
  estimator classes, not just different tuning constants -- call EDI with
  ``method="M"`` (its ``c=1.345`` default already matches ``HuberT``'s
  tuning constant) for a fair comparison against this baseline.
- ``HurdleCountModel``'s zero mechanism (a *censored* Poisson/NegBin
  likelihood for the zero/positive split) is structurally different from
  EDI's ``fast_hurdle_negbin`` (an independent logistic regression on the
  zero/positive indicator). Timing comparisons are still fair (same
  asymptotic complexity class); coefficient-parity checks are not, without
  accounting for the different parameterization.
- ``sksurv.linear_model.CoxPHSurvivalAnalysis.fit`` takes raw NumPy arrays
  (a 2D ``X`` and a structured ``y``); ``lifelines.WeibullAFTFitter.fit``
  requires a ``pandas.DataFrame`` -- that DataFrame construction is real,
  unavoidable overhead for the Weibull AFT bare-metal row, not padding to
  strip out.
- ``geepack::geeglm``'s ``corstr`` and ``statsmodels.genmod.
  generalized_estimating_equations.GEE``'s ``cov_struct`` both default to
  an independence working-correlation structure -- confirmed aligned, no
  mismatch to account for.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Optional

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

import statsmodels.api as sm
from statsmodels.discrete.discrete_model import NegativeBinomial
from statsmodels.discrete.count_model import (
    ZeroInflatedPoisson,
    ZeroInflatedNegativeBinomialP,
)
from statsmodels.discrete.truncated_model import HurdleCountModel
from statsmodels.othermod.betareg import BetaModel
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.robust.robust_linear_model import RLM
from statsmodels.genmod.generalized_estimating_equations import GEE

from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv

from lifelines import WeibullAFTFitter


@dataclass(frozen=True)
class Baseline:
    """One canonical Python baseline for an EDI model-fitting kernel."""

    package: str
    function: str
    fit: Callable[..., Any]
    notes: str = ""


# ── Continuous ───────────────────────────────────────────────────────────
def _ols_lstsq(X: np.ndarray, y: np.ndarray):
    """Bare-metal row: no model-object overhead, matches numpy.linalg.lstsq."""
    return np.linalg.lstsq(X, y, rcond=None)


def _ols_statsmodels(X: np.ndarray, y: np.ndarray):
    return sm.OLS(y, X).fit()


def _robust_regression(X: np.ndarray, y: np.ndarray, *, c: float = 1.345):
    # method="M" (plain HuberT) to match statsmodels RLM's own default --
    # see module docstring. EDI's own default is method="MM"; call EDI with
    # method="M" explicitly when pairing against this baseline.
    return RLM(y, X, M=sm.robust.norms.HuberT(t=c)).fit()


# ── Binary / proportion ──────────────────────────────────────────────────
def _glm_binomial(link: str):
    families = {
        "logit": sm.families.Binomial(),
        "probit": sm.families.Binomial(link=sm.families.links.Probit()),
        "log": sm.families.Binomial(link=sm.families.links.Log()),
        "identity": sm.families.Binomial(link=sm.families.links.Identity()),
    }
    fam = families[link]

    def fit(X: np.ndarray, y: np.ndarray, **kw):
        try:
            return sm.GLM(y, X, family=fam).fit(**kw)
        except Exception:
            # Non-canonical links (log, identity) can overshoot into an
            # infeasible probability region under pure IRLS from a naive
            # start. Fall back to a bounded quasi-Newton optimizer from a
            # safe start, matching R's glm.fit(..., start=) discipline for
            # these links.
            start = kw.pop("start_params", None)
            if start is None:
                p = X.shape[1]
                start = np.zeros(p)
                start[0] = -2.0 if link == "log" else 0.5
            return sm.GLM(y, X, family=fam).fit(
                start_params=start, method="lbfgs", maxiter=300, **kw
            )

    return fit


def _beta_regression(X: np.ndarray, y: np.ndarray, **kw):
    return BetaModel(y, X).fit(disp=0, **kw)


# ── Count ─────────────────────────────────────────────────────────────────
def _glm_poisson(X: np.ndarray, y: np.ndarray, **kw):
    return sm.GLM(y, X, family=sm.families.Poisson()).fit(**kw)


def _negbin(X: np.ndarray, y: np.ndarray, **kw):
    return NegativeBinomial(y, X).fit(disp=0, **kw)


def _zinb(X_cond: np.ndarray, X_infl: np.ndarray, y: np.ndarray, **kw):
    return ZeroInflatedNegativeBinomialP(y, X_cond, exog_infl=X_infl).fit(
        disp=0, **kw
    )


def _zap(X_cond: np.ndarray, X_infl: np.ndarray, y: np.ndarray, *, is_hurdle: bool = False, **kw):
    """EDI's fast_zero_augmented_poisson: is_hurdle picks zero-inflated vs
    hurdle Poisson, matching the EDI-side `is_hurdle` flag exactly."""
    if is_hurdle:
        return HurdleCountModel(y, X_cond, dist="poisson", zerodist="poisson").fit(
            disp=0, **kw
        )
    return ZeroInflatedPoisson(y, X_cond, exog_infl=X_infl).fit(disp=0, **kw)


def _hurdle_negbin(X_cond: np.ndarray, X_hurdle: np.ndarray, y: np.ndarray, **kw):
    # See module docstring: structurally different zero mechanism than
    # EDI's independent-logistic hurdle. X_hurdle is accepted for interface
    # symmetry with the EDI-side call but unused -- HurdleCountModel has no
    # separate zero-part design matrix, it reuses X_cond for both parts.
    del X_hurdle
    return HurdleCountModel(y, X_cond, dist="negbin", zerodist="negbin").fit(
        disp=0, **kw
    )


# ── Survival ─────────────────────────────────────────────────────────────
def _coxph(X: np.ndarray, time: np.ndarray, event: np.ndarray, **kw):
    y_struct = Surv.from_arrays(event=event.astype(bool), time=time)
    return CoxPHSurvivalAnalysis(**kw).fit(X, y_struct)


def _weibull_aft(X: np.ndarray, time: np.ndarray, event: np.ndarray, **kw):
    """Builds the required pandas.DataFrame from raw arrays -- this
    construction is itself real overhead relative to EDI's array-native
    call, not something to strip out of the timed region (see module
    docstring)."""
    p = X.shape[1]
    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(p)])
    df["__time"] = time
    df["__event"] = event
    return WeibullAFTFitter(**kw).fit(df, duration_col="__time", event_col="__event")


# ── Ordinal ──────────────────────────────────────────────────────────────
def _ordered(distr):
    def fit(X: np.ndarray, y: np.ndarray, **kw):
        return OrderedModel(y, X, distr=distr).fit(method="bfgs", disp=0, **kw)

    return fit


# ── GEE ──────────────────────────────────────────────────────────────────
def _gee(X: np.ndarray, y: np.ndarray, group_id: np.ndarray, *, family: str = "gaussian", **kw):
    families = {
        "gaussian": sm.families.Gaussian(),
        "binomial": sm.families.Binomial(),
        "poisson": sm.families.Poisson(),
    }
    return GEE(y, X, groups=group_id, family=families[family], cov_struct=None).fit(**kw)


# ── Registry ─────────────────────────────────────────────────────────────
# Keyed by the bound edi_kernels function name (see python/src/edi_kernels/
# __init__.py) so a benchmark runner can do BASELINES[kernel.__name__] or
# iterate this dict directly. Families with no confirmed clean baseline
# (GLMM/CLMM/LMM, adjacent-category/continuation-ratio/stereotype ordinal,
# zero-one-inflated beta, Weibull frailty, KK-combined estimators) are
# intentionally absent -- see "Baseline Gaps" in the spec doc, not silently
# missing rows.
BASELINES: dict[str, Baseline] = {
    "fast_ols": Baseline(
        "numpy", "numpy.linalg.lstsq", _ols_lstsq,
        "Bare-metal row. See fast_ols (statsmodels) for the OLS-model-object row.",
    ),
    "fast_ols (statsmodels)": Baseline("statsmodels", "OLS", _ols_statsmodels),
    "fast_robust_regression": Baseline(
        "statsmodels", "RLM(M=HuberT)", _robust_regression,
        "EDI defaults to method='MM' (two-stage S-then-M); call EDI with "
        "method='M' for a fair comparison against this baseline "
        "(single-stage HuberT M-estimation).",
    ),
    "fast_logistic_regression": Baseline(
        "statsmodels", "GLM(family=Binomial())", _glm_binomial("logit")
    ),
    "fast_probit_regression": Baseline(
        "statsmodels", "GLM(family=Binomial(link=Probit()))", _glm_binomial("probit")
    ),
    "fast_log_binomial_regression": Baseline(
        "statsmodels", "GLM(family=Binomial(link=Log()))", _glm_binomial("log"),
        "Log-link Binomial GLM can diverge under raw IRLS from a naive "
        "start; fit() falls back to lbfgs from a safe start if IRLS raises.",
    ),
    "fast_identity_binomial_regression": Baseline(
        "statsmodels", "GLM(family=Binomial(link=Identity()))", _glm_binomial("identity"),
        "Same IRLS-divergence fallback as the log-binomial row.",
    ),
    "fast_poisson_regression": Baseline(
        "statsmodels", "GLM(family=Poisson())", _glm_poisson
    ),
    "fast_neg_bin": Baseline(
        "statsmodels", "NegativeBinomial", _negbin,
        "Jointly estimates dispersion via MLE, matching glm.nb's semantics "
        "more closely than GLM(family=NegativeBinomial(alpha=fixed)).",
    ),
    "fast_zinb": Baseline(
        "statsmodels", "ZeroInflatedNegativeBinomialP", _zinb,
        "fit(X_cond, X_infl, y) -- pass EDI's Xc/Xz as X_cond/X_infl.",
    ),
    "fast_zero_augmented_poisson": Baseline(
        "statsmodels",
        "ZeroInflatedPoisson / HurdleCountModel(dist='poisson')",
        _zap,
        "fit(X_cond, X_infl, y, is_hurdle=...) -- is_hurdle mirrors EDI's "
        "own is_hurdle flag exactly (False: zero-inflated, True: hurdle).",
    ),
    "fast_hurdle_negbin": Baseline(
        "statsmodels", "HurdleCountModel(dist='negbin', zerodist='negbin')",
        _hurdle_negbin,
        "Structurally different hurdle mechanism than EDI's independent-"
        "logistic hurdle -- fair for timing, not for coefficient parity. "
        "See module docstring.",
    ),
    "fast_beta_regression": Baseline("statsmodels", "BetaModel", _beta_regression),
    "fast_coxph_regression": Baseline(
        "scikit-survival", "CoxPHSurvivalAnalysis", _coxph,
        "Preferred over lifelines.CoxPHFitter -- takes raw NumPy arrays, "
        "no DataFrame construction overhead.",
    ),
    "fast_weibull_regression": Baseline(
        "lifelines", "WeibullAFTFitter", _weibull_aft,
        "Requires a pandas.DataFrame; fit() builds it from raw arrays "
        "inside the timed region -- that construction is real overhead, "
        "not padding to strip out.",
    ),
    "fast_ordinal_regression": Baseline(
        "statsmodels", "OrderedModel(distr='logit')", _ordered("logit")
    ),
    "fast_ordinal_probit_regression": Baseline(
        "statsmodels", "OrderedModel(distr='probit')", _ordered("probit")
    ),
    "fast_ordinal_cauchit_regression": Baseline(
        "statsmodels", "OrderedModel(distr=scipy.stats.cauchy)",
        _ordered(sp_stats.cauchy),
        "Not the distr='...' string shortcut (only 'logit'/'probit' work "
        "as strings) -- pass a scipy.stats.rv_continuous instance instead.",
    ),
    "fast_ordinal_cloglog_regression": Baseline(
        "statsmodels", "OrderedModel(distr=scipy.stats.gumbel_l)",
        _ordered(sp_stats.gumbel_l),
        "scipy.stats.gumbel_l's CDF is exactly the cloglog inverse link "
        "(1 - exp(-exp(x))).",
    ),
    "fast_gee": Baseline(
        "statsmodels", "GEE(cov_struct=Independence())", _gee,
        "fit(X, y, group_id, family=...). Both sides default to an "
        "independence working-correlation structure -- confirmed aligned.",
    ),
}
