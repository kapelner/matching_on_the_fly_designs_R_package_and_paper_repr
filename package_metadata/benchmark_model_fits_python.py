#!/usr/bin/env python3
"""Python analog of benchmark/benchmark_model_fits.R.

Generates package_metadata/benchmark_model_fits_python.html: the same
Class/Response/EDI-Time/Canonical-Pkg/Canonical-Func/Canonical-Time/Speedup/
Timing-Pval table shape as benchmark_model_fits_R.html, same three-color row
coding (green = EDI faster + significant, grey = NA timing comparison, blue =
no canonical Python implementation exists), same underlying dataset spec
(N=1000 for most families, N=500 for survival, 4 continuous covariates + a
balanced binary treatment).

EDI has no Python bindings yet (see python_bindings_package_spec.md) — the
`python/` scaffold in this repo currently only exposes the `fast_math` scalar
utility kernels via pybind11, not the 33 model-fitting kernels this table
covers. So every row's "EDI Time" is NA today by construction, not because of
a failed fit: `time_edi()` below is a single, clearly-marked no-op that
returns NaN. Once a kernel is bound (`edi_kernels.fast_<name>`), wiring it in
is a one-line change to that row's `edi_kernel` callable — see the
`EDI_KERNELS` registry below, which already names the target function per row
so that future edit is mechanical, not a research task.

Canonical Python baselines were researched against
python_bindings_package_spec.md's "Python Baseline Registry" (2026-07-28/29)
and verified directly against the installed package versions before wiring
(see the venv this was developed/run against: numpy, scipy, pandas,
statsmodels, scikit-survival, lifelines — see requirements noted at the
bottom of this file). Families with no clean, actively-maintained Python
canonical equivalent are marked as Baseline Gaps (blue rows, NA canonical
timing too) rather than forcing a mismatched comparison — same discipline the
R report and the python_bindings_package_spec.md doc both already apply.
"""
import time
import warnings
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats as sstats

from scipy import special as sspecial

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
from statsmodels.regression.quantile_regression import QuantReg

from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.util import Surv

from lifelines import WeibullAFTFitter, KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import restricted_mean_survival_time

warnings.filterwarnings("ignore")  # matches R harness's options(warn = -1)
np.seterr(all="ignore")

rng = np.random.default_rng(42)

# ── Global config (mirrors benchmark_model_fits.R) ──────────────────────────
N_GLM = 1000
N_SURV = 500
P = 4                       # 4 continuous covariates + intercept + treatment
B_TIME = 30                 # cold timing replicates per side
TARGET_BATCH_MS = 200.0
MIN_RESOLVED_BATCH_MS = 10.0
MAX_INNER_REPS = 100_000

EDI_PY_AVAILABLE = False    # flip once `edi_kernels` (see spec doc) is installed


def time_edi(edi_kernel, fit_closure=None):
    """Bare-metal EDI timing. Always NA today (see module docstring).

    `fit_closure`, when EDI_PY_AVAILABLE flips True, should be a zero-arg
    callable performing exactly one EDI fit (same convention as the
    canonical-side closures below) built from `edi_kernels.<edi_kernel>`.
    """
    if not EDI_PY_AVAILABLE or fit_closure is None:
        return float("nan"), np.array([])
    return collect_timing_ms(fit_closure)


# ── Timing harness (mirrors collect_timing_ms() in benchmark_model_fits.R) ──
def collect_timing_ms(fn, times=B_TIME, target_batch_ms=TARGET_BATCH_MS,
                       max_inner_reps=MAX_INNER_REPS):
    fn()  # warm-up / validation call, discarded (matches R's convention)

    inner_reps = 1
    batch_ms = 0.0
    while inner_reps < max_inner_reps:
        t0 = time.perf_counter()
        for _ in range(inner_reps):
            fn()
        batch_ms = (time.perf_counter() - t0) * 1000.0
        if np.isfinite(batch_ms) and batch_ms >= min(target_batch_ms, MIN_RESOLVED_BATCH_MS):
            break
        inner_reps = min(max_inner_reps, inner_reps * 2)

    if np.isfinite(batch_ms) and batch_ms > 0:
        inner_reps = max(1, min(max_inner_reps, int(round(inner_reps * target_batch_ms / batch_ms))))

    vals = np.empty(times)
    for i in range(times):
        t0 = time.perf_counter()
        for _ in range(inner_reps):
            fn()
        vals[i] = (time.perf_counter() - t0) * 1000.0 / inner_reps

    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        return float("nan"), np.array([])
    return float(np.median(vals)), vals


def timing_ttest_pval(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return float("nan")
    try:
        return float(sstats.ttest_ind(x, y, equal_var=False).pvalue)
    except Exception:
        return float("nan")


# ── Data generators (mirrors generate_data()/dataset spec in the R harness) ─
def _design(n, p, w):
    X = rng.normal(size=(n, p))
    beta = rng.normal(scale=0.5, size=p)
    eta = X @ beta + 0.5 * w
    Xd = np.column_stack([np.ones(n), w, X])   # intercept + treatment + covariates
    Xo = np.column_stack([w, X])               # no-intercept design (ordinal/survival)
    return X, Xd, Xo, eta


def _balanced_treatment(n):
    n_t = n // 2
    w = np.array([1.0] * n_t + [0.0] * (n - n_t))
    rng.shuffle(w)
    return w


def data_binary(n=N_GLM, p=P, link="logit"):
    w = _balanced_treatment(n)
    X, Xd, Xo, eta = _design(n, p, w)
    if link == "log-binomial":
        prob = np.minimum(0.5, np.exp(eta - 2))
    elif link == "identity-binomial":
        prob = np.clip(0.5 * np.tanh(eta) + 0.5, 0.02, 0.98)
    else:
        prob = 1 / (1 + np.exp(-eta))
    y = rng.binomial(1, prob).astype(float)
    return dict(X=X, Xd=Xd, Xo=Xo, w=w, y=y)


def data_count(n=N_GLM, p=P):
    w = _balanced_treatment(n)
    X, Xd, Xo, eta = _design(n, p, w)
    y = rng.poisson(np.exp(eta * 0.3)).astype(float)
    return dict(X=X, Xd=Xd, Xo=Xo, w=w, y=y)


def data_continuous(n=N_GLM, p=P):
    w = _balanced_treatment(n)
    X, Xd, Xo, eta = _design(n, p, w)
    y = eta + rng.normal(scale=0.5, size=n)
    return dict(X=X, Xd=Xd, Xo=Xo, w=w, y=y)


def data_beta(n=N_GLM, p=P):
    w = _balanced_treatment(n)
    X, Xd, Xo, eta = _design(n, p, w)
    mu = 1 / (1 + np.exp(-eta * 0.3))
    phi = 10.0
    y = rng.beta(mu * phi, (1 - mu) * phi)
    return dict(X=X, Xd=Xd, Xo=Xo, w=w, y=np.clip(y, 1e-6, 1 - 1e-6))


def data_ordinal(n=N_GLM, p=P):
    w = _balanced_treatment(n)
    X, Xd, Xo, eta = _design(n, p, w)
    p1 = 1 / (1 + np.exp(eta - 1))
    p2 = 1 / (1 + np.exp(eta + 1)) - p1
    p3 = 1 - p1 - p2
    probs = np.clip(np.column_stack([p1, p2, p3]), 1e-6, None)
    probs = probs / probs.sum(1, keepdims=True)
    y = np.array([1 + rng.choice(3, p=probs[i]) for i in range(n)])
    return dict(X=X, Xd=Xd, Xo=Xo, w=w, y=y)


def data_survival(n=N_SURV, p=P):
    w = _balanced_treatment(n)
    X, Xd, Xo, eta = _design(n, p, w)
    y = rng.exponential(1 / np.exp(eta * 0.3))
    dead = rng.binomial(1, 0.8, n)
    return dict(X=X, Xd=Xd, Xo=Xo, w=w, y=y, dead=dead)


def data_survival_strat(n=N_SURV, p=P):
    # Low-cardinality strata grid (2x3), matching the R harness's
    # make_true_stratified_survival_data(): forces a genuinely stratified fit
    # rather than the unstratified fallback, with strata columns excluded
    # from the linear predictor (nested inside the strata label).
    w = _balanced_treatment(n)
    grid = np.array([(a, b) for a in range(2) for b in range(3)])
    idx = rng.integers(0, 6, n)
    s1 = grid[idx, 0].astype(float)
    s2 = grid[idx, 1].astype(float)
    X = rng.normal(size=(n, p))
    beta = np.array([0.45, -0.35, 0.20, -0.15])[:p]
    eta = 0.5 * w + X @ beta[: X.shape[1]]
    y = rng.exponential(1 / np.exp(0.3 * s1 - 0.2 * s2 + eta * 0.2))
    dead = rng.binomial(1, 0.8, n)
    df = pd.DataFrame({
        "y": y, "dead": dead, "treatment": w,
        **{f"x{i}": X[:, i] for i in range(X.shape[1])},
        "strata": pd.Series(s1).astype(int).astype(str) + "_" + pd.Series(s2).astype(int).astype(str),
    })
    return dict(df=df)


# ── Canonical fit closures, one per model family (grouped where the R
#    harness itself reuses one canonical call across several EDI classes,
#    e.g. quasi/modified/robust Poisson all reduce to the same glm.fit) ─────
def _glm(y, X, family):
    return lambda: sm.GLM(y, X, family=family).fit()


def build_glm_binomial(link="logit"):
    d = data_binary(link=link)
    fam = {
        "logit": sm.families.Binomial(),
        "probit": sm.families.Binomial(link=sm.families.links.Probit()),
        "log-binomial": sm.families.Binomial(link=sm.families.links.Log()),
        "identity-binomial": sm.families.Binomial(link=sm.families.links.Identity()),
    }[link]
    start = None
    if link == "log-binomial":
        start = np.r_[-2.0, np.zeros(d["Xd"].shape[1] - 1)]
    elif link == "identity-binomial":
        start = np.r_[0.5, np.zeros(d["Xd"].shape[1] - 1)]
    y, Xd = d["y"], d["Xd"]
    if start is None:
        return lambda: sm.GLM(y, Xd, family=fam).fit()

    # Pure-IRLS log/identity-link Binomial GLMs can overshoot into an
    # infeasible probability region mid-iteration and blow up (a known
    # instability of these non-canonical links, not specific to this
    # dataset). Determine once, outside the timed region, whether IRLS
    # converges from this start; if not, fall back to a bounded
    # quasi-Newton optimizer (matches how R's own glm.fit(family=
    # binomial(link="log"/"identity"), start=...) call needs a careful
    # start to avoid the same failure mode).
    try:
        sm.GLM(y, Xd, family=fam).fit(start_params=start, maxiter=300)
        method = "IRLS"
    except Exception:
        method = "lbfgs"
    return lambda: sm.GLM(y, Xd, family=fam).fit(start_params=start, method=method, maxiter=300)


def build_glm_poisson():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return lambda: sm.GLM(y, Xd, family=sm.families.Poisson()).fit()


def build_ols():
    d = data_continuous()
    y, Xd = d["y"], d["Xd"]
    return lambda: np.linalg.lstsq(Xd, y, rcond=None)


def build_lpm():
    d = data_binary()
    y, Xd = d["y"], d["Xd"]
    return lambda: np.linalg.lstsq(Xd, y, rcond=None)


def build_negbin():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return lambda: NegativeBinomial(y, Xd).fit(disp=0)


def build_beta_regr():
    d = data_beta()
    y, Xd = d["y"], d["Xd"]
    return lambda: BetaModel(y, Xd).fit(disp=0)


def build_fractional_logit():
    d = data_beta()
    y, Xd = d["y"], d["Xd"]
    return lambda: sm.GLM(y, Xd, family=sm.families.Binomial()).fit()


def build_ordered(distr):
    d = data_ordinal()
    y, Xo = d["y"], d["Xo"]
    return lambda: OrderedModel(y, Xo, distr=distr).fit(method="bfgs", disp=0)


def build_hurdle(dist):
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return lambda: HurdleCountModel(y, Xd, dist=dist, zerodist="poisson").fit(disp=0)


def build_zip():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return lambda: ZeroInflatedPoisson(y, Xd, exog_infl=Xd).fit(disp=0)


def build_zinb():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return lambda: ZeroInflatedNegativeBinomialP(y, Xd, exog_infl=Xd).fit(disp=0)


def build_weibull_aft():
    d = data_survival()
    df = pd.DataFrame({"y": d["y"], "dead": d["dead"], "treatment": d["w"],
                        **{f"x{i}": d["X"][:, i] for i in range(d["X"].shape[1])}})
    return lambda: WeibullAFTFitter().fit(df, duration_col="y", event_col="dead")


def build_coxph():
    d = data_survival()
    Xo, y, dead = d["Xo"], d["y"], d["dead"]
    ysurv = Surv.from_arrays(dead.astype(bool), y)
    return lambda: CoxPHSurvivalAnalysis().fit(Xo, ysurv)


def build_coxph_strat():
    d = data_survival_strat()
    df = d["df"]
    covs = [c for c in df.columns if c not in ("y", "dead", "strata")]
    return lambda: CoxPHFitter().fit(df[["y", "dead", "strata"] + covs], duration_col="y", event_col="dead", strata=["strata"])


def build_logrank():
    d = data_survival()
    y, dead, w = d["y"], d["dead"].astype(bool), d["w"]
    yt, yc = y[w == 1], y[w == 0]
    dt, dc = dead[w == 1], dead[w == 0]
    return lambda: logrank_test(yt, yc, dt, dc)


def build_km_median_diff():
    d = data_survival()
    y, dead, w = d["y"], d["dead"], d["w"]
    yt, yc = y[w == 1], y[w == 0]
    dt, dc = dead[w == 1], dead[w == 0]

    def fn():
        kmt = KaplanMeierFitter().fit(yt, dt)
        kmc = KaplanMeierFitter().fit(yc, dc)
        return kmt.median_survival_time_ - kmc.median_survival_time_
    return fn


def build_rmst_diff():
    d = data_survival()
    y, dead, w = d["y"], d["dead"], d["w"]
    yt, yc = y[w == 1], y[w == 0]
    dt, dc = dead[w == 1], dead[w == 0]
    tmax = min(yt.max(), yc.max())

    def fn():
        kmt = KaplanMeierFitter().fit(yt, dt)
        kmc = KaplanMeierFitter().fit(yc, dc)
        return restricted_mean_survival_time(kmt, t=tmax) - restricted_mean_survival_time(kmc, t=tmax)
    return fn


def build_hl_wilcox():
    d = data_continuous(n=500)  # scale = 0.5 in the R harness (O(n^2) pairwise diff)
    y, w = d["y"], d["w"]
    yt, yc = y[w == 1], y[w == 0]
    return lambda: np.median(np.subtract.outer(yt, yc))


def build_robust_regr():
    d = data_continuous()
    y, Xd = d["y"], d["Xd"]
    return lambda: RLM(y, Xd).fit()


def build_quantreg():
    d = data_continuous()
    y, Xd = d["y"], d["Xd"]
    return lambda: QuantReg(y, Xd).fit(q=0.5)


def build_gcomp_binary(effect):
    d = data_binary()
    y, Xd = d["y"], d["Xd"]

    def fn():
        mod = sm.GLM(y, Xd, family=sm.families.Binomial()).fit()
        X1 = Xd.copy(); X1[:, 1] = 1
        X0 = Xd.copy(); X0[:, 1] = 0
        r1 = mod.predict(X1).mean()
        r0 = mod.predict(X0).mean()
        return (r1 - r0) if effect == "RD" else (r1 / r0)
    return fn


def build_gcomp_fractional():
    d = data_beta()
    y, Xd = d["y"], d["Xd"]

    def fn():
        mod = sm.GLM(y, Xd, family=sm.families.Binomial()).fit()
        X1 = Xd.copy(); X1[:, 1] = 1
        X0 = Xd.copy(); X0[:, 1] = 0
        return mod.predict(X1).mean() - mod.predict(X0).mean()
    return fn


def build_gcomp_ordinal():
    d = data_ordinal()
    y, Xo = d["y"], d["Xo"]

    def fn():
        mod = OrderedModel(y, Xo, distr="logit").fit(method="bfgs", disp=0)
        X1 = Xo.copy(); X1[:, 0] = 1
        X0 = Xo.copy(); X0[:, 0] = 0
        probs1 = mod.model.predict(mod.params, exog=X1)
        probs0 = mod.model.predict(mod.params, exog=X0)
        score = np.arange(1, probs1.shape[1] + 1)
        return (probs1 @ score).mean() - (probs0 @ score).mean()
    return fn


# ── Wald (point + SE + p-value) canonical fit closures ──────────────────────
# Full-inference analog of the closures above: each returns (est, se, pval)
# for the treatment coefficient instead of just a fitted model object, still
# timed the same way (one call = one full fit + SE + Wald p-value). Most
# statsmodels-based families share one generic wrapper since statsmodels
# always names the (first, non-inflation) treatment column consistently
# ("x1" after an auto-detected "const", or plain "x1" when there is no
# intercept in the design) regardless of family — verified directly against
# GLM/NegativeBinomial/BetaModel/OrderedModel/Hurdle/ZIP/ZINB/RLM/QuantReg
# exog_names before relying on it here.
def _treatment_idx(exog_names):
    idx = [i for i, nm in enumerate(exog_names) if nm == "x1"]
    return idx[-1] if idx else 1


def _sm_wald_wrapper(fit_fn):
    def fn():
        res = fit_fn()
        i = _treatment_idx(res.model.exog_names)
        return float(res.params[i]), float(res.bse[i]), float(res.pvalues[i])
    return fn


def build_glm_binomial_wald(link="logit"):
    d = data_binary(link=link)
    fam = {
        "logit": sm.families.Binomial(),
        "probit": sm.families.Binomial(link=sm.families.links.Probit()),
        "log-binomial": sm.families.Binomial(link=sm.families.links.Log()),
        "identity-binomial": sm.families.Binomial(link=sm.families.links.Identity()),
    }[link]
    start = None
    if link == "log-binomial":
        start = np.r_[-2.0, np.zeros(d["Xd"].shape[1] - 1)]
    elif link == "identity-binomial":
        start = np.r_[0.5, np.zeros(d["Xd"].shape[1] - 1)]
    y, Xd = d["y"], d["Xd"]
    if start is None:
        return _sm_wald_wrapper(lambda: sm.GLM(y, Xd, family=fam).fit())
    try:
        sm.GLM(y, Xd, family=fam).fit(start_params=start, maxiter=300)
        method = "IRLS"
    except Exception:
        method = "lbfgs"
    return _sm_wald_wrapper(lambda: sm.GLM(y, Xd, family=fam).fit(start_params=start, method=method, maxiter=300))


def build_glm_poisson_wald():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: sm.GLM(y, Xd, family=sm.families.Poisson()).fit())


def build_ols_wald():
    d = data_continuous()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: sm.OLS(y, Xd).fit())


def build_lpm_wald():
    d = data_binary()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: sm.OLS(y, Xd).fit())


def build_negbin_wald():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: NegativeBinomial(y, Xd).fit(disp=0))


def build_beta_regr_wald():
    d = data_beta()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: BetaModel(y, Xd).fit(disp=0))


def build_fractional_logit_wald():
    d = data_beta()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: sm.GLM(y, Xd, family=sm.families.Binomial()).fit())


def build_ordered_wald(distr):
    d = data_ordinal()
    y, Xo = d["y"], d["Xo"]
    return _sm_wald_wrapper(lambda: OrderedModel(y, Xo, distr=distr).fit(method="bfgs", disp=0))


def build_hurdle_wald(dist):
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: HurdleCountModel(y, Xd, dist=dist, zerodist="poisson").fit(disp=0))


def build_zip_wald():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: ZeroInflatedPoisson(y, Xd, exog_infl=Xd).fit(disp=0))


def build_zinb_wald():
    d = data_count()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: ZeroInflatedNegativeBinomialP(y, Xd, exog_infl=Xd).fit(disp=0))


def build_robust_regr_wald():
    d = data_continuous()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: RLM(y, Xd).fit())


def build_quantreg_wald():
    d = data_continuous()
    y, Xd = d["y"], d["Xd"]
    return _sm_wald_wrapper(lambda: QuantReg(y, Xd).fit(q=0.5))


def build_weibull_aft_wald():
    d = data_survival()
    df = pd.DataFrame({"y": d["y"], "dead": d["dead"], "treatment": d["w"],
                        **{f"x{i}": d["X"][:, i] for i in range(d["X"].shape[1])}})

    def fn():
        aft = WeibullAFTFitter().fit(df, duration_col="y", event_col="dead")
        row = aft.summary.loc[("lambda_", "treatment")]
        return float(row["coef"]), float(row["se(coef)"]), float(row["p"])
    return fn


def build_coxph_wald():
    # sksurv's CoxPHSurvivalAnalysis doesn't expose SE/p-values directly;
    # lifelines' CoxPHFitter does (via .summary), so the Wald-table row
    # switches packages for the unstratified Cox row specifically — same
    # "documented fallback" pattern benchmark_model_fits.R itself uses when
    # a canonical package's fast bare-metal entry point doesn't expose the
    # variance the full-inference table needs.
    d = data_survival()
    df = pd.DataFrame({"y": d["y"], "dead": d["dead"], "treatment": d["w"],
                        **{f"x{i}": d["X"][:, i] for i in range(d["X"].shape[1])}})

    def fn():
        cph = CoxPHFitter().fit(df, duration_col="y", event_col="dead")
        row = cph.summary.loc["treatment"]
        return float(row["coef"]), float(row["se(coef)"]), float(row["p"])
    return fn


def _delta_method_se(theta, cov, estimand):
    eps = 1e-6
    grad = np.zeros_like(theta)
    for j in range(len(theta)):
        tp = theta.copy(); tp[j] += eps
        tm = theta.copy(); tm[j] -= eps
        grad[j] = (estimand(tp) - estimand(tm)) / (2 * eps)
    return float(np.sqrt(max(grad @ cov @ grad, 0.0)))


def build_gcomp_binary_wald(effect):
    d = data_binary()
    y, Xd = d["y"], d["Xd"]

    def estimand(theta):
        X1 = Xd.copy(); X1[:, 1] = 1
        X0 = Xd.copy(); X0[:, 1] = 0
        p1 = (1 / (1 + np.exp(-X1 @ theta))).mean()
        p0 = (1 / (1 + np.exp(-X0 @ theta))).mean()
        return (p1 - p0) if effect == "RD" else (p1 / p0)

    def fn():
        mod = sm.GLM(y, Xd, family=sm.families.Binomial()).fit()
        theta = np.asarray(mod.params)
        cov = np.asarray(mod.cov_params())
        est = estimand(theta)
        se = _delta_method_se(theta, cov, estimand)
        pval = float(2 * sstats.norm.sf(abs(est / se))) if se > 0 else float("nan")
        return est, se, pval
    return fn


def build_gcomp_fractional_wald():
    d = data_beta()
    y, Xd = d["y"], d["Xd"]

    def estimand(theta):
        X1 = Xd.copy(); X1[:, 1] = 1
        X0 = Xd.copy(); X0[:, 1] = 0
        p1 = (1 / (1 + np.exp(-X1 @ theta))).mean()
        p0 = (1 / (1 + np.exp(-X0 @ theta))).mean()
        return p1 - p0

    def fn():
        mod = sm.GLM(y, Xd, family=sm.families.Binomial()).fit()
        theta = np.asarray(mod.params)
        cov = np.asarray(mod.cov_params())
        est = estimand(theta)
        se = _delta_method_se(theta, cov, estimand)
        pval = float(2 * sstats.norm.sf(abs(est / se))) if se > 0 else float("nan")
        return est, se, pval
    return fn


def build_gcomp_ordinal_wald():
    d = data_ordinal()
    y, Xo = d["y"], d["Xo"]

    def fn():
        mod = OrderedModel(y, Xo, distr="logit").fit(method="bfgs", disp=0)
        theta = np.asarray(mod.params)
        cov = np.asarray(mod.cov_params())

        def estimand(th):
            X1 = Xo.copy(); X1[:, 0] = 1
            X0 = Xo.copy(); X0[:, 0] = 0
            probs1 = mod.model.predict(th, exog=X1)
            probs0 = mod.model.predict(th, exog=X0)
            score = np.arange(1, probs1.shape[1] + 1)
            return (probs1 @ score).mean() - (probs0 @ score).mean()

        est = estimand(theta)
        se = _delta_method_se(theta, cov, estimand)
        pval = float(2 * sstats.norm.sf(abs(est / se))) if se > 0 else float("nan")
        return est, se, pval
    return fn


# Wald-table-only build functions — for EDI classes that have no row in the
# point-estimate table at all (R's own wald_specs list includes several
# nonparametric-test classes that only ever appear in the Wald table, never
# the Results table, e.g. exact/rank tests that are inherently "one shot,
# already a full inference" computations with no separate point-estimate
# form). Verified against python_bindings_package_spec.md-style research:
# statsmodels.stats.proportion.test_proportions_2indep(method='score',
# correction=True) is documented by statsmodels itself as the Miettinen-
# Nurminen 1985 correction; confint_proportions_2indep(method='newcomb') is
# statsmodels' literal Newcombe method; lifelines' logrank_test(weightings=
# 'wilcoxon') is the Gehan (generalized Wilcoxon) weighting, matching R's
# survival::survdiff(rho=1).
def build_ttest_pooled():
    d = data_continuous()
    y, w = d["y"], d["w"]
    yt, yc = y[w == 1], y[w == 0]
    return lambda: sstats.ttest_ind(yt, yc, equal_var=True)


def build_mannwhitney():
    d = data_continuous(n=500)  # matches the point-estimate Wilcox row's scale
    y, w = d["y"], d["w"]
    yt, yc = y[w == 1], y[w == 0]
    return lambda: sstats.mannwhitneyu(yt, yc)


def build_lin_wald():
    d = data_continuous()
    y, w, X = d["y"], d["w"], d["X"]
    Xc = X - X.mean(0)
    Xint = np.column_stack([np.ones(len(y)), w, Xc, Xc * w[:, None]])
    return _sm_wald_wrapper(lambda: sm.OLS(y, Xint).fit())


def build_fisher_exact():
    d = data_binary()
    w, y = d["w"], d["y"]
    tab = np.array([[np.sum((w == 1) & (y == 1)), np.sum((w == 1) & (y == 0))],
                     [np.sum((w == 0) & (y == 1)), np.sum((w == 0) & (y == 0))]])
    return lambda: sstats.fisher_exact(tab)


def build_binom_diff(method):
    d = data_binary()
    w, y = d["w"], d["y"]
    x1, n1 = int(np.sum(y[w == 1])), int(np.sum(w == 1))
    x2, n2 = int(np.sum(y[w == 0])), int(np.sum(w == 0))
    kwargs = dict(method=method)
    if method == "score":
        kwargs["correction"] = True
    from statsmodels.stats.proportion import confint_proportions_2indep
    return lambda: confint_proportions_2indep(x1, n1, x2, n2, **kwargs)


def build_ridit():
    d = data_continuous(n=1600)  # matches the R harness's n_override for this row
    y, w = d["y"], d["w"]
    yy = np.round(y - y.min()).astype(int) + 1

    def fn():
        vals, counts = np.unique(yy, return_counts=True)
        cum = np.cumsum(counts)
        prev = np.r_[0, cum[:-1]]
        ridit_map = dict(zip(vals, (prev + 0.5 * counts) / len(yy)))
        r = np.array([ridit_map[v] for v in yy])
        return r[w == 1].mean() - r[w == 0].mean()
    return fn


def build_gehan_wilcox():
    d = data_survival()
    y, dead, w = d["y"], d["dead"], d["w"]
    yt, yc = y[w == 1], y[w == 0]
    dt, dc = dead[w == 1], dead[w == 0]
    return lambda: logrank_test(yt, yc, dt, dc, weightings="wilcoxon")


# WALD_SPECS mirrors R's wald_specs class set exactly (41 classes) — not
# simply "the same 40 as the point-estimate table minus a few exclusions".
# R's own Wald table includes several nonparametric-test classes with no
# point-estimate-table row at all (mean-diff pooled-var t-test, Wilcoxon
# rank-sum, Lin's covariate-adjusted estimator, Fisher exact, Miettinen-
# Nurminen/Newcombe risk-difference CIs, Jonckheere-Terpstra, Ridit, Gehan-
# Wilcoxon, KM median difference) and *excludes* several ordinal-link
# Baseline Gap classes and a few incidence/count rows that the
# point-estimate table does include (identity-binomial, modified Poisson,
# fractional logit, ordered probit, cauchit, cloglog) — those simply have no
# row here, matching R exactly rather than a "40 minus 3" guess.
WALD_SPECS = [
    ("continuous", "InferenceAllSimpleMeanDiffPooledVar", "scipy", "stats.ttest_ind(pooled)", "out of python-bindings scope (nonparametric-test kernel)", build_ttest_pooled, False),
    ("continuous", "InferenceAllSimpleWilcox", "scipy", "stats.mannwhitneyu", "out of python-bindings scope (nonparametric-test kernel, EDI:::wilcox_hl_point_estimate_cpp)", build_mannwhitney, False),
    ("continuous", "InferenceContinLin", "statsmodels", "OLS(interaction)+summary", "(none bound yet — Lin's estimator has no distinct EDI kernel in this table)", build_lin_wald, False),
    ("incidence", "InferenceIncidLogRegr", "statsmodels", "GLM(Binomial)+summary", "fast_logistic_regression_with_var", lambda: build_glm_binomial_wald("logit"), False),
    ("continuous", "InferenceContinOLS", "statsmodels", "OLS+summary", "fast_ols_with_var", build_ols_wald, False),
    ("count", "InferenceCountPoisson", "statsmodels", "GLM(Poisson)+summary", "fast_poisson_regression_with_var", build_glm_poisson_wald, False),
    ("survival", "InferenceSurvivalCoxPHRegr", "lifelines", "CoxPHFitter+summary", "fast_coxph_regression_with_var", build_coxph_wald, False),
    ("count", "InferenceCountNegBin", "statsmodels", "NegativeBinomial+summary", "fast_neg_bin_with_var", build_negbin_wald, False),
    ("proportion", "InferencePropBetaRegr", "statsmodels", "BetaModel+summary", "fast_beta_regression_with_var", build_beta_regr_wald, False),
    ("ordinal", "InferenceOrdinalPropOddsRegr", "statsmodels", "OrderedModel(logit)+summary", "fast_ordinal_regression_with_var", lambda: build_ordered_wald("logit"), False),
    ("count", "InferenceCountHurdlePoisson", "statsmodels", "HurdleCountModel(poisson)+summary", "fast_zero_augmented_poisson_with_var(is_hurdle=True)", lambda: build_hurdle_wald("poisson"), False),
    ("count", "InferenceCountZeroInflatedPoisson", "statsmodels", "ZeroInflatedPoisson+summary", "fast_zero_augmented_poisson_with_var(is_hurdle=False)", build_zip_wald, False),
    ("count", "InferenceCountZeroInflatedNegBin", "statsmodels", "ZeroInflatedNegativeBinomialP+summary", "fast_zinb_with_var", build_zinb_wald, False),
    ("count", "InferenceCountHurdleNegBin", "statsmodels", "HurdleCountModel(negbin)+summary", "fast_hurdle_negbin_with_var", lambda: build_hurdle_wald("negbin"), False),
    ("count", "InferenceCountQuasiPoisson", "statsmodels", "GLM(Poisson)+summary", "fast_quasipoisson_regression_with_var", build_glm_poisson_wald, False),
    ("survival", "InferenceSurvivalWeibullRegr", "lifelines", "WeibullAFTFitter+summary", "fast_weibull_regression (SE unexposed at bare-metal level)", build_weibull_aft_wald, False),
    ("continuous", "InferenceContinRobustRegr", "statsmodels", "RLM+summary", "fast_robust_regression_with_var", build_robust_regr_wald, False),
    ("continuous", "InferenceContinQuantileRegr", "statsmodels", "QuantReg+summary", "(none — EDI's own R class delegates to quantreg::rq; no distinct EDI kernel)", build_quantreg_wald, False),
    ("incidence", "InferenceIncidExactFisher", "scipy", "stats.fisher_exact", "out of python-bindings scope (nonparametric-test kernel)", build_fisher_exact, False),
    ("incidence", "InferenceIncidLogBinomial", "statsmodels", "GLM(Binomial, log link)+summary", "fast_log_binomial_regression_with_var", lambda: build_glm_binomial_wald("log-binomial"), False),
    ("incidence", "InferenceIncidProbitRegr", "statsmodels", "GLM(Binomial, probit link)+summary", "fast_probit_regression_with_var", lambda: build_glm_binomial_wald("probit"), False),
    ("incidence", "InferenceIncidMiettinenNurminenRiskDiff", "statsmodels", "stats.proportion.confint_proportions_2indep(score)", "out of python-bindings scope (nonparametric-test kernel, EDI:::mn_pvalue_cpp)", lambda: build_binom_diff("score"), False),
    ("incidence", "InferenceIncidNewcombeRiskDiff", "statsmodels", "stats.proportion.confint_proportions_2indep(newcomb)", "out of python-bindings scope (nonparametric-test kernel, EDI:::newcombe_independent_ci_cpp)", lambda: build_binom_diff("newcomb"), False),
    ("ordinal", "InferenceOrdinalAdjCatLogitRegr", None, None, "fast_adjacent_category_logit_with_var", None, True),
    ("ordinal", "InferenceOrdinalContRatioRegr", None, None, "fast_continuation_ratio_regression_with_var", None, True),
    ("survival", "InferenceSurvivalLogRank", "lifelines", "statistics.logrank_test", "out of python-bindings scope (nonparametric-test kernel); identical call as the point-estimate row — logrank_test always returns a full test", build_logrank, False),
    ("survival", "InferenceSurvivalStratCoxPHRegr", "lifelines", "CoxPHFitter(strata=)", "fast_coxph_regression_with_var (prebuilt/stratified path); identical call as the point-estimate row — lifelines computes the variance unconditionally", build_coxph_strat, False),
    ("incidence", "InferenceIncidRiskDiff", "statsmodels", "OLS+summary (LPM)", "fast_ols_with_var", build_lpm_wald, False),
    ("count", "InferenceCountRobustPoisson", "statsmodels", "GLM(Poisson)+summary", "fast_poisson_regression_with_var", build_glm_poisson_wald, False),
    ("proportion", "InferencePropGCompMeanDiff", "statsmodels", "GLM(Binomial)+gcomp+delta-method SE", "fast_logistic_regression + gcomp utility (out of python-bindings scope)", build_gcomp_fractional_wald, False),
    ("incidence", "InferenceIncidGCompRiskRatio", "statsmodels", "GLM(Binomial)+gcomp(RR)+delta-method SE", "fast_logistic_regression + gcomp utility (out of python-bindings scope)", lambda: build_gcomp_binary_wald("RR"), False),
    ("incidence", "InferenceIncidGCompRiskDiff", "statsmodels", "GLM(Binomial)+gcomp(RD)+delta-method SE", "fast_logistic_regression + gcomp utility (out of python-bindings scope)", lambda: build_gcomp_binary_wald("RD"), False),
    ("ordinal", "InferenceOrdinalGCompMeanDiff", "statsmodels", "OrderedModel(logit)+gcomp+delta-method SE", "fast_ordinal_regression + gcomp utility (out of python-bindings scope)", build_gcomp_ordinal_wald, False),
    ("ordinal", "InferenceOrdinalJonckheereTerpstraTest", None, None, "out of python-bindings scope (nonparametric-test kernel, EDI:::fast_jonckheere_terpstra_cpp)", None, True),
    ("ordinal", "InferenceOrdinalRidit", "numpy", "manual ridit computation", "out of python-bindings scope (nonparametric-test kernel, EDI:::fast_ridit_analysis_cpp)", build_ridit, False),
    ("survival", "InferenceSurvivalGehanWilcox", "lifelines", "statistics.logrank_test(weightings='wilcoxon')", "out of python-bindings scope (nonparametric-test kernel)", build_gehan_wilcox, False),
    ("survival", "InferenceSurvivalKMDiff", "lifelines", "KaplanMeierFitter(median)+CI", "out of python-bindings scope (nonparametric-test kernel, EDI:::get_survival_stat_diff); identical call as the point-estimate row — lifelines computes the CI unconditionally", build_km_median_diff, False),
    ("incidence", "InferenceIncidKKCondLogitPlusGLMMOneLik", None, None, "fast_clogit_plus_glmm", None, True),
    ("count", "InferenceCountKKCondPoissonOneLik", None, None, "fast_cpoisson_combined_with_var", None, True),
    ("survival", "InferenceSurvivalKKWeibullFrailtyOneLik", None, None, "fast_weibull_frailty", None, True),
    ("proportion", "InferencePropZeroOneInflatedBetaRegr", None, None, "fast_zero_one_inflated_beta", None, True),
]

# R's own wald_specs list has no row for InferenceSurvivalRestrictedMeanDiff
# either (lifelines' restricted_mean_survival_time is a point-value-only
# function; a variance would need a separate, differently-shaped bootstrap
# computation) — this mirrors that omission exactly rather than being a
# unilateral Python-side choice.
WALD_EXCLUDED = {
    "InferenceSurvivalRestrictedMeanDiff": "Not in R's own Wald table either: lifelines' restricted_mean_survival_time (like R's survival:::survmean) returns only a point value; a variance requires a separate bootstrap-based computation, a different computational profile than every other Wald row here.",
}


# ── Model registry — one row per EDI Inference class, same 40-class subset
#    benchmark_model_fits.R's point-estimate table covers, in the same
#    "targeted subset" spirit. `edi_kernel` documents the future
#    edi_kernels.* binding this row will time once it exists (or, for
#    kernels permanently out of the model-fitting binding scope per
#    python_bindings_package_spec.md — nonparametric tests, gcomp utilities
#    — a note saying so instead of a function name). ─────────────────────
MODEL_SPECS = [
    # response, cls, pkg, func, edi_kernel, build, no_canonical
    ("incidence", "InferenceIncidLogRegr", "statsmodels", "GLM(Binomial)", "fast_logistic_regression", lambda: build_glm_binomial("logit"), False),
    ("continuous", "InferenceContinOLS", "numpy", "linalg.lstsq", "fast_ols", build_ols, False),
    ("count", "InferenceCountPoisson", "statsmodels", "GLM(Poisson)", "fast_poisson_regression", build_glm_poisson, False),
    ("survival", "InferenceSurvivalCoxPHRegr", "scikit-survival", "CoxPHSurvivalAnalysis", "fast_coxph_regression", build_coxph, False),
    ("count", "InferenceCountNegBin", "statsmodels", "NegativeBinomial", "fast_neg_bin", build_negbin, False),
    ("proportion", "InferencePropBetaRegr", "statsmodels", "BetaModel", "fast_beta_regression", build_beta_regr, False),
    ("ordinal", "InferenceOrdinalPropOddsRegr", "statsmodels", "OrderedModel(logit)", "fast_ordinal_regression", lambda: build_ordered("logit"), False),
    ("count", "InferenceCountHurdlePoisson", "statsmodels", "HurdleCountModel(poisson)", "fast_zero_augmented_poisson(is_hurdle=True)", lambda: build_hurdle("poisson"), False),
    ("count", "InferenceCountZeroInflatedPoisson", "statsmodels", "ZeroInflatedPoisson", "fast_zero_augmented_poisson(is_hurdle=False)", build_zip, False),
    ("count", "InferenceCountZeroInflatedNegBin", "statsmodels", "ZeroInflatedNegativeBinomialP", "fast_zinb", build_zinb, False),
    ("count", "InferenceCountHurdleNegBin", "statsmodels", "HurdleCountModel(negbin)", "fast_hurdle_negbin", lambda: build_hurdle("negbin"), False),
    ("count", "InferenceCountQuasiPoisson", "statsmodels", "GLM(Poisson)", "fast_poisson_regression", build_glm_poisson, False),
    ("survival", "InferenceSurvivalWeibullRegr", "lifelines", "WeibullAFTFitter", "fast_weibull_regression", build_weibull_aft, False),
    ("continuous", "InferenceContinRobustRegr", "statsmodels", "RLM", "fast_robust_regression", build_robust_regr, False),
    ("continuous", "InferenceContinQuantileRegr", "statsmodels", "QuantReg", "(none — EDI's own R class delegates to quantreg::rq.fit; no distinct EDI kernel)", build_quantreg, False),
    ("proportion", "InferencePropFractionalLogit", "statsmodels", "GLM(Binomial, fractional y)", "fast_logistic_regression", build_fractional_logit, False),
    ("incidence", "InferenceIncidLogBinomial", "statsmodels", "GLM(Binomial, log link)", "fast_log_binomial_regression", lambda: build_glm_binomial("log-binomial"), False),
    ("incidence", "InferenceIncidProbitRegr", "statsmodels", "GLM(Binomial, probit link)", "fast_probit_regression", lambda: build_glm_binomial("probit"), False),
    ("incidence", "InferenceIncidBinomialIdentityRiskDiff", "statsmodels", "GLM(Binomial, identity link)", "fast_identity_binomial_regression", lambda: build_glm_binomial("identity-binomial"), False),
    ("ordinal", "InferenceOrdinalAdjCatLogitRegr", None, None, "fast_adjacent_category_logit", None, True),
    ("ordinal", "InferenceOrdinalContRatioRegr", None, None, "fast_continuation_ratio_regression", None, True),
    ("ordinal", "InferenceOrdinalOrderedProbitRegr", "statsmodels", "OrderedModel(probit)", "fast_ordinal_probit_regression", lambda: build_ordered("probit"), False),
    ("ordinal", "InferenceOrdinalCloglogRegr", None, None, "fast_ordinal_cloglog_regression", None, True),
    ("ordinal", "InferenceOrdinalCauchitRegr", None, None, "fast_ordinal_cauchit_regression", None, True),
    ("survival", "InferenceSurvivalLogRank", "lifelines", "statistics.logrank_test", "out of python-bindings scope (nonparametric-test kernel, EDI:::fast_logrank_stats_cpp)", build_logrank, False),
    ("continuous", "InferenceAllSimpleWilcox", "numpy", "median(HL pairwise diff)", "out of python-bindings scope (nonparametric-test kernel, EDI:::wilcox_hl_point_estimate_cpp)", build_hl_wilcox, False),
    ("incidence", "InferenceIncidModifiedPoisson", "statsmodels", "GLM(Poisson)", "fast_poisson_regression", build_glm_poisson, False),
    ("survival", "InferenceSurvivalStratCoxPHRegr", "lifelines", "CoxPHFitter(strata=)", "fast_coxph_regression (prebuilt/stratified path)", build_coxph_strat, False),
    ("incidence", "InferenceIncidRiskDiff", "numpy", "linalg.lstsq (LPM)", "fast_ols", build_lpm, False),
    ("survival", "InferenceSurvivalKMDiff", "lifelines", "KaplanMeierFitter(median)", "out of python-bindings scope (nonparametric-test kernel, EDI:::get_survival_stat_diff)", build_km_median_diff, False),
    ("survival", "InferenceSurvivalRestrictedMeanDiff", "lifelines", "utils.restricted_mean_survival_time", "out of python-bindings scope (nonparametric-test kernel, EDI:::get_survival_stat_diff)", build_rmst_diff, False),
    ("count", "InferenceCountRobustPoisson", "statsmodels", "GLM(Poisson)", "fast_poisson_regression", build_glm_poisson, False),
    ("proportion", "InferencePropGCompMeanDiff", "statsmodels", "GLM(Binomial)+gcomp", "fast_logistic_regression + gcomp utility (out of python-bindings scope)", build_gcomp_fractional, False),
    ("incidence", "InferenceIncidGCompRiskRatio", "statsmodels", "GLM(Binomial)+gcomp(RR)", "fast_logistic_regression + gcomp utility (out of python-bindings scope)", lambda: build_gcomp_binary("RR"), False),
    ("incidence", "InferenceIncidGCompRiskDiff", "statsmodels", "GLM(Binomial)+gcomp(RD)", "fast_logistic_regression + gcomp utility (out of python-bindings scope)", lambda: build_gcomp_binary("RD"), False),
    ("ordinal", "InferenceOrdinalGCompMeanDiff", "statsmodels", "OrderedModel(logit)+gcomp", "fast_ordinal_regression + gcomp utility (out of python-bindings scope)", build_gcomp_ordinal, False),
    ("incidence", "InferenceIncidKKCondLogitPlusGLMMOneLik", None, None, "fast_clogit_plus_glmm", None, True),
    ("count", "InferenceCountKKCondPoissonOneLik", None, None, "fast_cpoisson_combined", None, True),
    ("survival", "InferenceSurvivalKKWeibullFrailtyOneLik", None, None, "fast_weibull_frailty", None, True),
    ("proportion", "InferencePropZeroOneInflatedBetaRegr", None, None, "fast_zero_one_inflated_beta", None, True),
]

NO_CANONICAL_NOTE = {
    "InferenceOrdinalAdjCatLogitRegr": "No identified Python package implements the adjacent-category logit link (R uses VGAM::vglm(acat())).",
    "InferenceOrdinalContRatioRegr": "No identified Python package implements the continuation-ratio link (R uses VGAM::vglm(cratio())).",
    "InferenceOrdinalCloglogRegr": "statsmodels.miscmodels.ordinal_model.OrderedModel's distr= argument only officially supports 'probit'/'logit' strings; cloglog is not a documented/tested option, so this is treated as a gap rather than an unverified custom-distribution hack.",
    "InferenceOrdinalCauchitRegr": "Same as cloglog: OrderedModel's distr= only documents 'probit'/'logit'.",
    "InferenceIncidKKCondLogitPlusGLMMOneLik": "KK combined (matched-pair + reservoir) joint-likelihood estimator: no canonical analog in either R or Python (see python_bindings_package_spec.md Baseline Gaps).",
    "InferenceCountKKCondPoissonOneLik": "KK combined (matched-pair + reservoir) joint-likelihood estimator: no canonical analog in either R or Python.",
    "InferenceSurvivalKKWeibullFrailtyOneLik": "Weibull AFT with shared log-normal frailty: no clean Python package (R side notes only a partial/PH-parameterized frailtypack analog).",
    "InferencePropZeroOneInflatedBetaRegr": "Zero-one-inflated beta regression: no canonical package in either R or Python.",
    "InferenceOrdinalJonckheereTerpstraTest": "No identified Python package implements the Jonckheere-Terpstra trend test (checked scipy, statsmodels, scikit-posthocs; R itself needs the specialty clinfun package, not base R either — no equally-specialized Python package was found).",
}


# ── Utility / math-kernel canonical baselines ───────────────────────────────
# Mirrors the R harness's "Utility / Math Kernel Performance" section: EDI's
# internal fast_* scalar math kernels (used inside the NegBin/Beta/ZINB/
# Hurdle likelihoods and probit cold-start heuristics above) benchmarked
# against scipy's/numpy's vectorized equivalents, over a fixed-length
# vector, mirroring benchmark/fast_math_utils_bench.cpp's methodology.
# EDI's python/ fast_math pybind11 stub currently only binds
# fast_pchisq_upper (see python/cpp/bindings_fast_math.cpp) — none of the 8
# functions below — so this table's EDI column is NA for the same "not
# bound yet" reason as the two tables above, not a different one.
N_UTIL = 5000


def build_util_digamma():
    x = rng.uniform(0.5, 50, N_UTIL)
    return lambda: sspecial.digamma(x)


def build_util_trigamma():
    x = rng.uniform(0.5, 50, N_UTIL)
    return lambda: sspecial.polygamma(1, x)


def build_util_lgamma():
    x = rng.uniform(0.5, 50, N_UTIL)
    return lambda: sspecial.gammaln(x)


def build_util_lbeta():
    a = rng.uniform(1, 20, N_UTIL)
    b = rng.uniform(1, 20, N_UTIL)
    return lambda: sspecial.betaln(a, b)


def build_util_qnorm():
    p = rng.uniform(0.001, 0.999, N_UTIL)
    return lambda: sstats.norm.ppf(p)


def build_util_log_pnorm():
    x = rng.normal(0, 3, N_UTIL)
    return lambda: sstats.norm.logcdf(x)


def build_util_log_dnorm():
    x = rng.normal(0, 3, N_UTIL)
    return lambda: sstats.norm.logpdf(x)


def build_util_dnbinom_mu():
    x = rng.poisson(5, N_UTIL)
    size, mu = 2.0, 5.0
    p_nb = size / (size + mu)
    return lambda: sstats.nbinom.logpmf(x, size, p_nb)


UTILITY_SPECS = [
    ("utility", "fast_digamma", "scipy", "special.digamma", "fast_math.fast_digamma (not yet bound — only fast_pchisq_upper is)", build_util_digamma, False),
    ("utility", "fast_trigamma", "scipy", "special.polygamma(1,.)", "fast_math.fast_trigamma (not yet bound)", build_util_trigamma, False),
    ("utility", "fast_lgamma", "scipy", "special.gammaln", "fast_math.fast_lgamma (not yet bound)", build_util_lgamma, False),
    ("utility", "fast_lbeta", "scipy", "special.betaln", "fast_math.fast_lbeta (not yet bound)", build_util_lbeta, False),
    ("utility", "fast_qnorm", "scipy", "stats.norm.ppf", "fast_math.fast_qnorm (not yet bound)", build_util_qnorm, False),
    ("utility", "fast_log_pnorm", "scipy", "stats.norm.logcdf", "fast_math.fast_log_pnorm (not yet bound)", build_util_log_pnorm, False),
    ("utility", "fast_log_dnorm", "scipy", "stats.norm.logpdf", "fast_math.fast_log_dnorm (not yet bound)", build_util_log_dnorm, False),
    ("utility", "fast_dnbinom_mu", "scipy", "stats.nbinom.logpmf", "fast_math.fast_dnbinom_mu (not yet bound)", build_util_dnbinom_mu, False),
]


# ── Formatting helpers (mirror format_ms/format_pval/row_bg_color in R) ────
def format_ms(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if x < 0.01:
        return f"{x:.3g}"
    return f"{x:.2f}"


def format_pval(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{x:.3g}"


def format_pval_stars(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    if x < 0.001:
        return "***"
    if x < 0.01:
        return "**"
    if x < 0.05:
        return "*"
    return ""


def row_bg_color(speedup, pval, no_canonical=False):
    if no_canonical:
        return "#cfe2ff"
    if speedup is None or not np.isfinite(speedup) or pval is None or not np.isfinite(pval):
        return "#eceff1"
    if pval < 0.05 and speedup > 1:
        return "#d9fdd3"
    return ""


# ── Runner ───────────────────────────────────────────────────────────────
def run_one(response, cls, pkg, func, edi_kernel, build, no_canonical):
    print(f"Benchmarking {cls}...")
    if no_canonical:
        return dict(cls=cls, response=response, edi_ms=float("nan"), pkg="None",
                    func="", canonical_ms=float("nan"),
                    speedup=float("nan"), pval=float("nan"), no_canonical=True,
                    edi_kernel=edi_kernel)

    edi_ms, edi_samples = time_edi(edi_kernel, None)

    try:
        fn = build()
        can_ms, can_samples = collect_timing_ms(fn)
    except Exception as e:
        print(f"  Canonical Error ({cls}): {e!r}")
        can_ms, can_samples = float("nan"), np.array([])

    speedup = (can_ms / edi_ms) if (np.isfinite(can_ms) and np.isfinite(edi_ms) and edi_ms > 0) else float("nan")
    pval = timing_ttest_pval(edi_samples, can_samples)
    return dict(cls=cls, response=response, edi_ms=edi_ms, pkg=pkg, func=func,
                canonical_ms=can_ms, speedup=speedup, pval=pval, no_canonical=False,
                edi_kernel=edi_kernel)


def run_all(specs, label):
    print(f"\n=== {label} ({len(specs)} rows) ===")
    rows = [run_one(*spec) for spec in specs]
    rows.sort(key=lambda r: (r["response"], r["cls"]))
    return rows


def render_table_html(rows):
    out = []
    for r in rows:
        bg = row_bg_color(r["speedup"], r["pval"], r["no_canonical"])
        style = f' style="background-color: {bg};"' if bg else ""
        speed_str = "NA" if not np.isfinite(r["speedup"]) else f"{r['speedup']:.2f}x"
        out.append(
            f'<tr{style}><td>{r["cls"]}</td><td>{r["response"]}</td>'
            f'<td>{format_ms(r["edi_ms"])}</td><td>{r["pkg"]}</td><td>{r["func"]}</td>'
            f'<td>{format_ms(r["canonical_ms"])}</td><td>{speed_str}</td>'
            f'<td>{format_pval(r["pval"])}</td><td>{format_pval_stars(r["pval"])}</td></tr>'
        )
    return "\n".join(out)


def render_gap_list_html(rows, notes):
    return "\n".join(
        f"<li><code>{r['cls']}</code> (future kernel <code>{r['edi_kernel']}</code>): "
        f"{notes.get(r['cls'], '')}</li>"
        for r in rows if r["no_canonical"]
    )


LEGEND_HTML = """<p class="legend">
<span style="background:#d9fdd3"></span>EDI faster, significant &nbsp;
<span style="background:#eceff1"></span>NA timing comparison &nbsp;
<span style="background:#cfe2ff"></span>no canonical Python implementation
</p>"""

TABLE_HEAD_HTML = ('<table>\n<thead><tr><th>Class</th><th>Response</th><th>EDI Time (ms)</th>'
                    '<th>Canonical Pkg</th><th>Canonical Func</th><th>Canonical Time (ms)</th>'
                    '<th>Speedup</th><th>Timing Pval</th><th></th></tr></thead>\n<tbody>\n')


def counts(rows):
    n_gap = sum(1 for r in rows if r["no_canonical"])
    n_ok = sum(1 for r in rows if not r["no_canonical"] and np.isfinite(r["canonical_ms"]))
    n_fail = len(rows) - n_gap - n_ok
    return n_ok, n_gap, n_fail


def main():
    model_rows = run_all(MODEL_SPECS, "Point-estimate")
    wald_rows = run_all(WALD_SPECS, "Wald (point + SE + p-value)")
    util_rows = run_all(UTILITY_SPECS, "Utility / math kernels")

    m_ok, m_gap, m_fail = counts(model_rows)
    w_ok, w_gap, w_fail = counts(wald_rows)
    u_ok, u_gap, u_fail = counts(util_rows)

    model_gap_html = render_gap_list_html(model_rows, NO_CANONICAL_NOTE)
    wald_gap_html = render_gap_list_html(wald_rows, NO_CANONICAL_NOTE)
    wald_excluded_html = "\n".join(f"<li><code>{cls}</code>: {note}</li>" for cls, note in WALD_EXCLUDED.items())

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    versions = []
    for modname in ("numpy", "scipy", "pandas", "statsmodels", "sksurv", "lifelines"):
        try:
            mod = __import__(modname)
            versions.append(f"{modname} {getattr(mod, '__version__', '?')}")
        except Exception:
            versions.append(f"{modname} (not installed)")

    html = f"""<!doctype html>
<html><head><meta charset="utf-8">
<title>EDI Python Baseline Benchmarks</title>
<style>
:root {{ color-scheme: light dark; }}
body {{ font-family: -apple-system, Segoe UI, Helvetica, Arial, sans-serif; max-width: 1200px; margin: 2rem auto; padding: 0 1rem; line-height: 1.5; }}
@media (prefers-color-scheme: dark) {{
  body {{ background: #0d1117; color: #e6edf3; }}
  a {{ color: #2f81f7; }}
  table, th, td {{ border-color: #30363d !important; }}
  code {{ background: #161b22; }}
}}
@media (prefers-color-scheme: light) {{
  body {{ background: #ffffff; color: #24292f; }}
  code {{ background: #f6f8fa; }}
}}
table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
th, td {{ border: 1px solid #d0d7de; padding: 4px 8px; text-align: left; font-size: 0.9rem; }}
th {{ background: rgba(127,127,127,0.15); }}
tr[style] td {{ color: #111; }}
code {{ padding: 1px 4px; border-radius: 4px; }}
.legend span {{ display: inline-block; width: 1em; height: 1em; margin-right: 0.4em; vertical-align: middle; border: 1px solid #8888; }}
nav a {{ margin-right: 1rem; }}
</style>
</head><body>
<h1>EDI Python Baseline Benchmarks</h1>
<p><em>Generated: {generated}</em></p>
<nav><a href="#results">Point-estimate</a><a href="#wald">Wald (full inference)</a><a href="#utility">Utility / math kernels</a></nav>

<p>Python analog of <a href="benchmark_model_fits_R.html">benchmark_model_fits_R.html</a> —
same three tables (point-estimate, Wald/full-inference, utility math kernels), same table shape,
same three-color row coding. Produced by <code>package_metadata/benchmark_model_fits_python.py</code>.</p>

<h2>Status: EDI Python bindings do not exist yet</h2>
<p>EDI's C++ model-fitting kernels have no Python bindings yet, and the <code>python/</code>
scaffold's nascent <code>fast_math</code> pybind11 stub only binds one utility function
(<code>fast_pchisq_upper</code>, not used by any row below) — see
<code>package_metadata/python_bindings_package_spec.md</code>. Every row's
<strong>EDI Time</strong> is therefore <code>NA</code> in all three tables today by
construction, not because of a failed fit: <code>time_edi()</code> in this
script is a single, clearly-marked no-op returning <code>NaN</code>. Once a
kernel is bound, wiring it in is a one-line change per row (the
<code>edi_kernel</code> field in each spec list already names the
target <code>edi_kernels.*</code> function) — flip
<code>EDI_PY_AVAILABLE = True</code> and pass the real fit closure. Until
then, these reports document what the <strong>canonical Python baselines cost
today</strong>, and which model families/functions have <strong>no canonical Python
baseline at all</strong> (regardless of EDI bindings).</p>

<h2>Benchmark Dataset Specification</h2>
<ul>
<li><strong>Sample size (N):</strong> 1,000 subjects for most families; 500 for survival families and the Wilcoxon Hodges-Lehmann row (O(n&sup2;) pairwise-difference computation, matching the R harness's <code>scale = 0.5</code>); 5,000 elements for the utility-function vectors.</li>
<li><strong>Predictors (p):</strong> intercept + a balanced binary treatment + 4 continuous covariates ~ Normal(0, 1); covariate coefficients ~ Normal(0, 0.5); treatment coefficient fixed at 0.5.</li>
<li><strong>Same family-generation formulas as <code>benchmark_model_fits.R</code></strong> (logistic/log-binomial/identity-binomial links, Poisson mean model, Beta(mu*phi, (1-mu)*phi) proportions, exponential survival times with ~20% censoring, 3-level ordinal construction).</li>
<li><strong>Stratified Cox exception:</strong> a low-cardinality (2x3) strata grid is injected before outcome generation, matching the R harness, so the row exercises a genuinely stratified fit.</li>
</ul>

<h2>Methodology</h2>
<ul>
<li><strong>Bare-metal canonical timing:</strong> each canonical row constructs the model object and calls its lowest-level <code>.fit()</code>/equivalent directly on pre-built NumPy arrays (or, where the package requires it — <code>lifelines</code>, <code>statsmodels.OrderedModel</code> internals — a pre-built <code>pandas.DataFrame</code>), inside the timed region only; data generation happens once, outside the timed closure.</li>
<li><strong>Point-estimate vs. Wald tables:</strong> the point-estimate table times a bare fit (matching each row's fastest available canonical entry point — <code>lstsq</code> instead of <code>OLS().fit()</code> where that's a real, distinct fast path); the Wald table times a full fit that also produces the treatment coefficient's standard error and two-sided p-value (<code>.bse</code>/<code>.pvalues</code> off the same fitted result for most statsmodels families; a package switch to <code>lifelines.CoxPHFitter</code> for the unstratified-Cox row specifically, since <code>scikit-survival</code>'s bare-metal <code>CoxPHSurvivalAnalysis</code> doesn't expose a variance; a finite-difference delta-method SE for the four G-computation rows, off <code>.cov_params()</code>). The Wald table's row set matches R's <code>wald_specs</code> exactly (41 classes) rather than reusing the point-estimate table's 40 — it adds several nonparametric-test-only classes with no point-estimate-table row at all (pooled-variance t-test, Wilcoxon rank-sum, Lin's estimator, Fisher exact, Miettinen-Nurminen/Newcombe risk-difference CIs, Jonckheere-Terpstra, Ridit, Gehan-Wilcoxon, KM median difference) and omits several point-estimate-table rows that R's own Wald table never covers (identity-binomial, modified Poisson, fractional logit, ordered probit, cauchit, cloglog). One row, restricted-mean-survival-time difference, is excluded from the Wald table for the same reason R's own Wald table excludes it — see below.</li>
<li><strong>Utility-function timing:</strong> each row calls the scipy/numpy vectorized function once over a fixed-length input vector, mirroring <code>benchmark/fast_math_utils_bench.cpp</code>'s apples-to-apples vectorized-vs-vectorized discipline (a scalar-loop binding would unfairly penalize a C++ side that hasn't been written yet).</li>
<li><strong>Averaging:</strong> medians over {B_TIME} cold timing samples via an adaptive-batch <code>time.perf_counter</code> harness (mirrors the R harness's adaptive <code>system.time</code> split, target {int(TARGET_BATCH_MS)}ms/batch).</li>
<li><strong>Significance:</strong> Welch's two-sample t-test (<code>scipy.stats.ttest_ind(..., equal_var=False)</code>) between the EDI and canonical timing replicate distributions — currently always <code>NA</code> since the EDI side has no samples yet in any of the three tables.</li>
<li><strong>Row highlighting:</strong> light green = <code>Speedup &gt; 1</code> and <code>Timing Pval &lt; 0.05</code>; light grey = <code>NA</code> timing comparison (EDI not bound yet, or a fit failed); light blue = no canonical Python implementation exists at all for this model family/function.</li>
<li><strong>Package versions this report was run against:</strong> {", ".join(versions)}.</li>
</ul>

<h2 id="results">Point-Estimate Results ({m_ok} of {len(model_rows)} rows timed, {m_gap} Baseline Gaps)</h2>
{LEGEND_HTML}
<p>Families with no clean, actively-maintained Python canonical equivalent — an absent/NA comparison here is more honest than a mismatched substitute baseline (same discipline <code>python_bindings_package_spec.md</code> and the R report apply):</p>
<ul>
{model_gap_html}
</ul>
{TABLE_HEAD_HTML}
{render_table_html(model_rows)}
</tbody>
</table>
<p><small>{m_ok} of {len(model_rows)} rows have a working canonical Python timing; {m_fail} canonical fit(s) failed on this run (see console log); {m_gap} are documented Baseline Gaps.</small></p>

<h2 id="wald">Wald Test Performance / Full Inference ({w_ok} of {len(wald_rows)} rows timed, {w_gap} Baseline Gaps)</h2>
{LEGEND_HTML}
<p>Same Baseline Gap families as the point-estimate table above, plus three rows excluded outright (no variance computation, not a Baseline Gap in the usual sense):</p>
<ul>
{wald_excluded_html}
</ul>
<details><summary>Baseline Gap families (same as point-estimate table)</summary><ul>
{wald_gap_html}
</ul></details>
{TABLE_HEAD_HTML}
{render_table_html(wald_rows)}
</tbody>
</table>
<p><small>{w_ok} of {len(wald_rows)} rows have a working canonical Python timing; {w_fail} canonical fit(s) failed on this run (see console log); {w_gap} are documented Baseline Gaps; {len(WALD_EXCLUDED)} rows from the point-estimate table are excluded entirely (see above).</small></p>

<h2 id="utility">Utility / Math Kernel Performance ({u_ok} of {len(util_rows)} functions timed)</h2>
{LEGEND_HTML}
<p>EDI's internal <code>fast_*</code> scalar math kernels (<code>fast_digamma</code>, <code>fast_trigamma</code>, <code>fast_lgamma</code>, <code>fast_lbeta</code>, <code>fast_qnorm</code>, <code>fast_log_pnorm</code>, <code>fast_log_dnorm</code>, <code>fast_dnbinom_mu</code>) vs. their scipy/numpy vectorized equivalents, over a length-{N_UTIL} vector. None of these 8 are bound in the <code>python/</code> <code>fast_math</code> stub yet (only <code>fast_pchisq_upper</code> is), so every row is grey today.</p>
{TABLE_HEAD_HTML}
{render_table_html(util_rows)}
</tbody>
</table>
<p><small>{u_ok} of {len(util_rows)} functions have a working canonical timing; {u_fail} failed on this run.</small></p>

<p><small>See <code>package_metadata/python_bindings_package_spec.md</code> for the full kernel-binding plan.</small></p>
</body></html>
"""

    out_path = "package_metadata/benchmark_model_fits_python.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"\nWrote {out_path}: point-estimate {m_ok}/{len(model_rows)} ok ({m_gap} gaps, {m_fail} failed), "
          f"Wald {w_ok}/{len(wald_rows)} ok ({w_gap} gaps, {w_fail} failed), "
          f"utility {u_ok}/{len(util_rows)} ok ({u_fail} failed).")


if __name__ == "__main__":
    main()
