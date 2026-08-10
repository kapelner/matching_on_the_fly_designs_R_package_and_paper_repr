# Sequential Inference for `DesignSeqOneByOne*` Designs

Generated: 2026-08-09

## Scope

This is a feasibility and design-space exploration, not an implementation
spec. It asks: can EDI support sequential hypothesis testing / sequential
confidence-interval construction (i.e. valid interim analyses performed
repeatedly *while* a `DesignSeqOneByOne*` experiment is still enrolling
subjects), and if so, how would that be wired into the existing
`Design` / `Inference` object pair?

Related background:

- [cold_starts.md](cold_starts.md)
- [warm_starts.html](warm_starts.html)
- [extending-edi-r6.md](extending-edi-r6.md)

## 1. What "sequential" already means in EDI vs. what is being asked for

EDI already uses the word "sequential" for the *design* side:
`DesignSeqOneByOne` (`EDI/R/design_seq_one_by_one_abstract.R:12`) and its 13
subclasses (Bernoulli, iBCRD, Efron, Atkinson, Urn, RandomBlockSize, SPBR,
PocockSimon, KK14, KK21, KK21stepwise, plus `DesignCustomSequential`) admit
subjects one at a time via `add_one_subject()` /
`add_one_subject_to_experiment_and_assign()`
(`EDI/R/design_seq_one_by_one_abstract.R:88,151`), and responses arrive
one at a time via `Design$add_one_subject_response()`
(`EDI/R/design_abstract.R:145`).

But "sequential" here refers only to *randomization* — the treatment
assignment mechanism adapts on the fly to the accruing covariate/outcome
history (matching, minimization, biased coins, urns). Nothing in the
package currently performs **sequential inference** in the classical sense:
repeatedly testing/estimating the treatment effect *as data accrue*, with
error-rate or coverage guarantees that hold *simultaneously across looks*.
Today, every `Inference*` class assumes the trial is over
(`des_obj$assert_all_responses_recorded()` at
`EDI/R/inference_all_abstract.R:56`) and produces exactly one p-value/CI
from a frozen dataset. That's a fixed-sample, single-look analysis, full
stop.

This document is about closing that gap: reusing this same
matching-on-the-fly machinery to do valid *interim* analyses, i.e. "we have
enrolled and observed t < n subjects — what can we say about the treatment
effect right now, and can we do that again at t+1, t+2, ... without
inflating type-I error or breaking coverage?"

## 2. A load-bearing fact: `Inference$new()` already half-works mid-trial

`Inference$initialize()` gates construction on
`des_obj$assert_all_responses_recorded()`
(`EDI/R/inference_all_abstract.R:56`), which calls
`assert_all_subjects_arrived()` (`EDI/R/design_abstract.R:250`):

```r
assert_all_subjects_arrived = function(){
  if (should_run_asserts()) {
    if (private$fixed_sample & private$t < private$n){
      stop("This experiment is incomplete as all n subjects haven't arrived yet.")
    }
  }
}
```

and then checks `sum(!is.na(private$y)) != length(private$w)`
(`EDI/R/design_abstract.R:258`) — i.e. "does every *currently assigned*
subject have a recorded response?"

Two consequences that are easy to miss:

- If the design was constructed **with a fixed `n`** (`fixed_sample = TRUE`),
  `Inference$new(des_obj)` is hard-blocked until all `n` subjects have
  arrived. No interim analysis is possible without bypassing this check.
- If the design was constructed **with `n = NULL`** (open-ended accrual —
  fully supported today; see `Design$initialize`,
  `EDI/R/design_abstract.R:82`), then `assert_all_subjects_arrived()` is a
  no-op, and `assert_all_responses_recorded()` only checks that responses
  are in for the subjects enrolled *so far*. In that mode,
  **`InferenceAllSimpleMeanDiff$new(des_obj)` (or any other concrete
  `Inference*` class) already runs successfully at any intermediate `t`**,
  today, with zero code changes.

That second bullet is the practical foothold for everything below: an
open-accrual `DesignSeqOneByOne*` object plus repeated
`Inference*$new(des_obj)` calls at increasing `t` already gives you a raw
per-look estimate/SE/p-value/CI. What's missing is not data plumbing — it's
the sequential-testing *layer* on top (spending the right amount of alpha
per look, or an anytime-valid statistic), plus a place to keep the
look-by-look bookkeeping.

Caveat specific to this package's matching designs: `DesignSeqOneByOneKK14`
(and KK21/Atkinson/SPBR) use `self$get_n()` internally as a scale parameter
for the assignment rule itself, not just for a completeness check — e.g. the
Mahalanobis-distance match/reservoir cutoff
`T_cutoff_sq = rank_prev * (n - 1) / (n - rank_prev) * F_crit`
(`EDI/R/design_seq_one_by_one_KK14.R:60-63`). `get_n()` returns `private$t`
when `fixed_sample = FALSE` (`EDI/R/design_abstract.R:353`), so running a KK
design "open-ended" to dodge the completeness assertion silently changes
the matching threshold at every step (it keeps chasing the current `t`
instead of a fixed target). This is a statistical side effect, not just a
bookkeeping inconvenience, and needs to be either (a) documented as a known
limitation of interim analysis on `n = NULL` KK designs, or (b) fixed by
giving these designs a separate "target/planned n for the matching rule"
argument, distinct from "hard stop enrollment at n" — see §6.

## 3. Where the state lives today

All of the state needed to reconstruct "the trial as of time t" already
lives in `Design`'s private environment (`EDI/R/design_abstract.R`), and is
shared by reference with any `DesignSeqOneByOne*` subclass:

| field | what it is | growth |
|---|---|---|
| `private$t` | number of subjects enrolled so far | scalar counter |
| `private$Xraw` | raw covariates, one row per subject | `rbindlist`-appended in `add_one_subject()` |
| `private$w` | treatment assignments (0/1) | appended in `add_one_subject_to_experiment_and_assign()` |
| `private$y`, `private$y_original` | responses | appended in `add_one_subject_response()` / `add_all_subject_responses()` |
| `private$dead` | censoring indicator | same |
| `private$X` | imputed/model-matrix covariates | recomputed lazily once `t > ncol(Xraw) + 2` |
| `private$m` (KK subclasses only) | match/reservoir index per subject | updated inside `assign_wt()` |
| `private$all_subject_data_cache` | **memoized per-`t`** slices of `X`/rank/etc. | keyed by `as.character(private$t)` (`EDI/R/design_abstract.R:591`) |

Note `all_subject_data_cache` is already keyed by `t` — but it is a
*performance* memo (avoids recomputing QR/rank/scaling on the growing `X`
matrix), not an analysis ledger. It does not, and should not, store
p-values, test statistics, or "have we looked here already" flags. It's
worth naming explicitly so nobody mistakes it for sequential-analysis
bookkeeping.

`Inference$initialize()` (`EDI/R/inference_all_abstract.R:34`) takes a
**snapshot**, not a live view: `private$y = private$des_obj_priv_int$y`,
`private$w = private$des_obj_priv_int$w`, `private$dead = ...$dead` are
ordinary R value copies (copy-on-modify). `private$des_obj = des_obj` and
`private$des_obj_priv_int = des_obj$.__enclos_env__$private` *are* live R6
environment references, but the `Inference` object's own working data
(`private$y`, `private$w`, `private$X`, and everything derived from them
into `private$cached_values`) is frozen at construction time. If the
`des_obj` continues to accrue subjects after an `Inference` object is
built, that `Inference` object does **not** see the new subjects — it just
keeps answering questions about the `t` it was built at.

## 4. Statistical approaches to sequential testing / CI construction

These are not mutually exclusive; different `Inference*` families in this
package are better suited to different approaches.

1. **Group-sequential testing with alpha-spending** (Pocock, O'Brien-Fleming,
   Lan-DeMets spending functions). Pre-specify (or spend flexibly at)
   `K` looks; each look consumes a slice of total alpha as a function of
   information fraction `t / n_planned`. Requires a fixed target `n` to
   define the information fraction — in tension with the `n = NULL`
   loophole in §2, so this approach pushes toward fixing §6 (allow interim
   analysis on fixed-`n` designs) rather than routing around it.

2. **Bayesian sequential monitoring** (posterior probability boundaries,
   e.g. stop for efficacy when `P(effect > 0 | data_t) > 0.99`). This is
   arguably the most natural fit for the package's existing
   `InferenceBayesianBootstrap` family
   (`EDI/R/inference_all_abstract_bayesian_bootstrap.R:7`) since a posterior
   (or posterior-like bootstrap) distribution can be recomputed at each `t`
   with no formal alpha-spending correction — the well-known cost is that
   frequentist operating characteristics (type-I error under repeated
   looks) must be checked by simulation rather than assumed, e.g. via
   `SimulationFramework` (`EDI/R/simulations_framework.R:237`).

3. **Anytime-valid inference / confidence sequences** (test martingales,
   e-processes, mixture SPRT; Howard et al. 2021, Waudby-Smith & Ramdas
   2020ish "confidence sequences"). These give a sequence of CIs
   `CI_t` such that `P(mu in CI_t for all t) >= 1 - alpha` — valid under
   *arbitrary*, even data-dependent, stopping/peeking, with no
   pre-specified look schedule or spending function needed. This is the
   most robust option and best matches "peek whenever you want," but needs
   a martingale construction tailored to each estimand (mean difference,
   risk difference, etc.) — more implementation work than 1–2.

4. **Sequential/conditional randomization tests**, exploiting infrastructure
   this package already has. The rand-CI/rand-test family inverts
   `des_obj$draw_ws_according_to_design(r)`
   (`EDI/R/design_abstract.R:346`) to build a reference distribution of the
   test statistic under the *actual* sequential assignment rule (Bernoulli,
   KK match/reservoir, Atkinson, etc.) conditional on the observed
   covariate history. Because `draw_ws_raw()` for every
   `DesignSeqOneByOne*` subclass replays the assignment rule forward using
   only information available through `t` (see
   `EDI/R/design_seq_one_by_one_KK14.R:106-112`,
   `EDI/R/design_seq_one_by_one_abstract.R:186-195`), the *same* machinery
   already used for a single terminal randomization test can, in principle,
   be re-run at every interim `t` to get a valid conditional permutation
   p-value **at that t**. What it does *not* give you for free is joint
   (across-`t`) type-I error control — that still needs a spending function
   (approach 1) or a martingale wrapper (approach 3) layered around the
   per-look conditional p-value.

5. **Stochastic curtailment / conditional power** for operational
   futility/efficacy stopping. Not a formal error-rate-controlling method
   on its own (typically paired with 1), but cheap to compute from any
   `Inference*` object's current estimate + SE, and useful as a
   `SequentialMonitor` diagnostic even before a full spending-function
   implementation exists.

The central statistical hazard to flag prominently, regardless of which
approach is chosen: **naively calling `Inference*$new(des_obj)` at several
values of `t` and reporting "significant" the first time `p < 0.05` is
invalid** (classic repeated-significance-testing / "peeking" inflation).
The `assert_all_responses_recorded()` loophole in §2 makes it *technically
easy* to do this today, which is exactly why it must not be exposed to
users without a sequential-testing wrapper on top — see §7.

## 5. How this would weave into `Inference*` objects

Recommendation: **`Inference*` objects stay single-look, ephemeral, and
naive.** Do not turn any `Inference*` class into a long-lived, mutable,
"keeps state across looks" object. Instead, add a new orchestration layer
above `Design`/`Inference`, analogous in spirit to how `SimulationFramework`
(`EDI/R/simulations_framework.R:237`) already sits *above* `Design`/
`Inference` pairs to run many replications without either of those classes
knowing about replication.

Proposed new class, e.g. `SequentialMonitor` (name TBD) — but see §7:
**this object must hold no state that isn't cheaply reconstructible**, so
that the only thing a user is ever required to persist between R sessions
is `des_obj` itself, the object they are already accruing subjects into.

```r
monitor = SequentialMonitor$new(
  des_obj,                       # the live, still-enrolling Design object
  inference_class = InferenceAllSimpleMeanDiff,
  inference_args = list(),       # extra args forwarded to Inference$new()
  spending_function = obf_spending(),   # or "bayesian", "confidence_sequence", ...
  planned_n = 200,                # information-fraction denominator (see S6)
  alpha = 0.05
)
# ^ cheap, deterministic, holds only *policy*, not history -- see S7.
# Safe to reconstruct with the same arguments in every new R session; it
# reads/writes its ledger through des_obj, not through its own fields.

# called by the trial operator whenever a look is desired (subject-count- or
# calendar-triggered -- the monitor does not decide *when* to look)
result_t = monitor$look()
# -> constructs a throwaway Inference*, computes estimate/test/CI,
#    updates alpha spent / martingale value, appends one row to des_obj's
#    ledger (S7), discards the Inference object, returns a decision:
#    continue / stop-efficacy / stop-futility
```

Each call to `look()`:

1. Constructs a **fresh** `Inference*` object from the design's *current*
   state (bypassing/relaxing `assert_all_subjects_arrived()` for fixed-`n`
   designs — see §6).
2. Computes the per-look statistic via the normal `Inference*` public API
   (`compute_estimate()`, `compute_asymp_two_sided_pval()`, etc. —
   `EDI/R/inference_all_abstract.R:160,142`).
3. Feeds that statistic into the sequential-testing layer (spending
   function boundary check, or e-process/martingale update), replaying
   prior looks read back from `des_obj$get_analysis_events()` (§7) so a
   freshly-reconstructed monitor picks up exactly where the last session
   left off.
4. Appends one entry to `des_obj`'s ledger via
   `des_obj$record_analysis_event(...)` (see §7) and returns a decision.
5. **Drops the `Inference*` object.** It is not stored, not reused, not
   mutated. The next `look()` builds a brand-new one from the (by-then
   further-accrued) `des_obj`.

## 6. Should `Inference*` objects be deleted after each look? Yes.

Treat every interim `Inference*` instance as write-once/throwaway. Reasons:

- **They cache stale-by-construction data.** `private$cached_values` and
  `private$y`/`private$w`/`private$X` are frozen snapshots (§3). There is
  no dependency-tracking or invalidation mechanism that would let an old
  `Inference` object notice the backing `des_obj` grew from `t` to `t+1`.
  Reusing one across looks would silently return `t`-stale results.
- **This matches existing usage elsewhere in the codebase.** Bootstrap,
  jackknife, and permutation infrastructure already follow a
  construct-use-discard pattern rather than mutating a shared `Inference`
  object in place (`Inference$duplicate()`,
  `EDI/R/inference_all_abstract.R:191`, exists precisely to hand out cheap
  disposable copies for exactly this kind of resampling use).
- **No explicit cleanup is required.** R6 objects are ordinary environments;
  an abandoned interim `Inference` object is garbage-collected once
  dereferenced. `private$finalize` (`EDI/R/inference_all_abstract.R:349`)
  is a no-op placeholder for cluster-ownership bookkeeping, not a general
  destructor — nothing needs to be added there for this feature.

The one place this recommendation should be revisited: if a chosen
sequential method needs warm-starting an optimizer's fit across looks for
performance (e.g. likelihood-based `Inference*` classes re-fitting from
scratch at every interim look on a growing `n` could get expensive late in
a long trial), it may be worth extracting *just* the warm-start state
(`private$fit_warm_start*`, `EDI/R/inference_all_abstract.R:437-501`) and
threading it through `SequentialMonitor` to the next look's fresh
`Inference*$new()`, rather than keeping the whole object alive. That's a
performance optimization, not a reason to keep the object itself around.

## 7. Where does "we've done inference at t" get recorded?

Currently: nowhere. `all_subject_data_cache` (§3) is a numerical-recompute
cache, not an audit trail — it doesn't store p-values, boundaries, or
"already looked at this t" flags, and it's invalidated/irrelevant to
monitoring logic.

**Revised recommendation (superseding an earlier draft of this section):
put the ledger on `Design`/`DesignSeqOneByOne`, not on `SequentialMonitor`.**

The original tradeoff considered here was framed as "layering purity"
(keep `Design` ignorant of `Inference`/hypothesis testing) vs. "convenience"
(keep the ledger next to `des_obj`). That framing missed the deciding
factor: **persistence across R sessions.** EDI has no built-in save/load API
today — the only `saveRDS`/`readRDS` calls in the package are internal to
`SimulationFramework`'s replication cache
(`EDI/R/simulations_framework.R:2429,2463`) — so between-session
persistence of a `DesignSeqOneByOne*` object is already just base-R
`saveRDS(des_obj, ...)` / `readRDS(...)`, done directly by the user on the
object they're actively calling `add_one_subject_response()` on. If the
ledger instead lives inside a separate `SequentialMonitor` object:

- The user must remember to *also* persist `monitor`, not just `des_obj`.
- Worse, R6/environment reference sharing only survives serialization
  *within a single `saveRDS()` call*. If `des_obj` and `monitor` are ever
  saved in separate calls (the natural thing to do, since `des_obj` is
  updated far more often — every subject — than a monitor's config), then
  on reload `monitor$des_obj` is a stale, disconnected copy from whenever
  `monitor` was last saved, not the live `des_obj` the user just loaded.
  This silently produces a monitor whose ledger and decisions are wrong
  relative to the design it now points at, with no error raised.

Putting the ledger on `Design` sidesteps this entirely: it's just more
private state on the one object the user was always going to persist, and
it comes along for free through `duplicate()`/`clone()`
(`EDI/R/design_abstract.R:441`) the same way `private$Xraw`/`private$y`/
`private$w` do. Concretely, add to the shared `DesignSeqOneByOne` (or
`Design`) private state:

```r
private = list(
  analysis_log = list(),   # list of event records, each tagged with t
  ...
)
```

with two new generic public methods:

```r
record_analysis_event = function(payload){   # payload: any list, e.g.
  private$analysis_log[[length(private$analysis_log) + 1L]] =           # list(t=.., estimate=.., se=..,
    c(list(t = private$t, timestamp = Sys.time()), payload)             #      boundary=.., decision=..)
  invisible(self)
},
get_analysis_events = function(){ private$analysis_log }
```

`Design` does not need to know what a p-value or a spending function is —
it just appends and returns opaque, monitor-supplied payloads keyed by
`t`. That keeps the layering objection mostly intact (Design still knows
nothing about *statistics*, only that "analysis events happen and get
logged"), while giving up the harder, correctness-relevant claim from the
original draft (that this state should live off of `Design` altogether).

`SequentialMonitor` then becomes the thin object described in §5: it holds
only *policy* (`inference_class`, `spending_function`, `planned_n`,
`alpha`) — nothing that would be lost or go stale if the object itself is
never serialized and is instead freshly reconstructed with the same
arguments at the start of every R session. Its `look()` method reads prior
looks back via `des_obj$get_analysis_events()` to resume a spending
function or martingale exactly where it left off, and writes new ones via
`des_obj$record_analysis_event()`. This adds one new field and two new
methods to the `DesignSeqOneByOne` abstract class (one place, inherited by
all 13 concrete subclasses) — not the 13-file change the original draft of
this section implied.

## 8. Phased plan

- [ ] TODO-1: **Phase 0 (no core changes, proof of concept):** demonstrate
  group-sequential O'Brien-Fleming monitoring on top of
  `InferenceAllSimpleMeanDiff` for a continuous response, using an
  `n = NULL` (open-accrual) `DesignSeqOneByOneBernoulli` design and manual
  per-look `Inference*$new(des_obj)` calls, validated by simulation via
  `SimulationFramework`. This exercises §2's existing loophole directly and
  requires zero package changes — good for validating the statistical
  approach before investing in plumbing.
- [ ] TODO-2: **Phase 1:** add `record_analysis_event()` / `get_analysis_events()` to
  `DesignSeqOneByOne` (§7), then build the thin, statelessly-reconstructible
  `SequentialMonitor` class (§5) on top: spending-function plug-in point,
  `look()` API, reading/writing its ledger through `des_obj`.
- [ ] TODO-3: **Phase 2:** add an explicit, narrowly-scoped bypass for
  `assert_all_subjects_arrived()` (e.g. an `interim = TRUE` argument
  threaded through `Inference$initialize()` and
  `assert_all_responses_recorded()`) so fixed-`n` designs support interim
  analysis the same way `n = NULL` designs already accidentally do, without
  weakening the completeness check for ordinary terminal analyses. This
  also removes the need for the `n = NULL` KK-matching-threshold caveat in
  §2, since fixed-`n` designs would no longer need to fake open-ended
  accrual to unlock interim looks.
- [ ] TODO-4: **Phase 3:** audit each concrete `DesignSeqOneByOne*` assignment rule
  (KK14, KK21, KK21stepwise, Atkinson, SPBR, PocockSimon, RandomBlockSize)
  for whether its `draw_ws_raw()`/`assign_wt()` genuinely conditions only
  on history through `t` with no forward-looking dependence on final `n`
  beyond the matching-threshold caveat already identified — a prerequisite
  for approach 4 (sequential randomization tests) to be valid at arbitrary
  interim `t`.
- [ ] TODO-5: **Phase 4:** extend beyond group-sequential/randomization approaches to
  anytime-valid confidence sequences (approach 3) for the highest-value
  `Inference*` families (mean difference, risk difference), since these
  remove the requirement to pre-commit to a look schedule at all.

## 9. Non-goals

- Not proposing to make any `Inference*` class itself sequential/stateful.
- Not proposing to change what `Design$add_one_subject_response()` /
  `add_all_subject_responses()` do — accrual mechanics are unaffected.
- Not claiming the `n = NULL` loophole in §2 is a *supported* way to do
  interim analysis today — it's a description of current (undocumented,
  accidental) behavior that motivates and de-risks Phase 0, not a
  recommendation to build on it long-term given the KK-threshold caveat.
- Not selecting a single "winning" sequential-testing approach in §4 — that
  choice should be driven by which `Inference*` families and response types
  are prioritized, and is out of scope for this feasibility document.
