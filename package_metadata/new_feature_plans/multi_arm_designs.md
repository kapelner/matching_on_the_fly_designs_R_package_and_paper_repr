# Multi-Arm (K > 2) Designs

Generated: 2026-08-09

## Scope

Every `Design*`/`Inference*` class in EDI today assumes exactly two arms:
treatment and control. This is a feasibility and design-space exploration of
what it would take to support `K > 2` arms — e.g. a three-dose trial, a
factorial arm with more than two cells, a multi-treatment comparative
effectiveness trial — and it follows a phasing the user proposed directly:

- **Phase 1**: generalize `Design*` assignment mechanics to `K` arms for
  every design family *except* the KK matching-on-the-fly family, and give
  `Inference*` a way to consume a `K`-arm design by having the user pick two
  arms and subset down to them, so every existing `Inference*` class works
  completely unmodified.
- **Phase 2**: generalize the KK matching-on-the-fly designs
  (`DesignSeqOneByOneKK14`/`KK21`/`KK21stepwise`, `DesignFixedBinaryMatch`,
  `DesignFixedMatchingGreedyPairSwitching`) to `K` arms — deferred because
  their core mechanism is pairwise, not because it is conceptually
  uninteresting.
- **Phase 3**: native multi-arm `Inference`, either as (a) an orchestration
  layer that runs multiple existing two-arm `Inference*` objects over
  Phase 1's subsets and combines them, or (b) genuinely new `Inference*`
  classes whose estimand is a `K`-1-dimensional contrast instead of a
  scalar.

This document does not select a final design for any phase. It maps exactly
where the two-arm assumption is load-bearing today (§1), which concrete
`Design*` subclasses generalize cheaply vs. not at all (§2–§4), what the
`Inference` contract's scalar shape costs a native multi-arm implementation
(§5–§6), and proposes a concrete phased plan with file-level scope (§8).

Related background: [extending-edi-r6.md](extending-edi-r6.md) (the
external extension contract, whose scalar `estimate`/`se`/`df` shape is
directly relevant to §6); [sequential_inference.md](sequential_inference.md)
(a structurally similar feasibility document — its §5 "thin orchestration
layer above Design/Inference" pattern is reused in §6a of this document).

## 1. Where "two arms" is load-bearing today

The two-arm assumption is not one check, it's an encoding convention
repeated at every layer:

- **`overwrite_all_subject_assignments(w)`** (`EDI/R/design_abstract.R:234-242`)
  hard-asserts `w` contains only `{-1, +1}`, and `DesignFixed`'s
  `assign_w_to_all_subjects()` applies the same `(w+1)/2` storage transform
  when a precomputed `{-1,+1}` vector is supplied directly
  (`EDI/R/design_fixed_abstract.R:51-53`):
  ```r
  assertIntegerish(w, lower = -1, upper = 1, any.missing = FALSE, len = private$t)
  if (any(!(w %in% c(-1L, 1L)))) {
    stop("overwrite_all_subject_assignments: w must contain only -1 (control) or +1 (treated).")
  }
  private$w = (as.numeric(w) + 1L) / 2L
  ```
  Internally `private$w` is always stored as `{0, 1}`; the public API only
  ever speaks `{-1, +1}`.
- **`Design$get_w()`** (`design_abstract.R:339-341`) and
  **`draw_ws_according_to_design()`** (`design_abstract.R:346-349`) both
  apply `2L * x - 1L` to go from the internal `{0,1}` storage to the public
  `{-1,+1}` convention. This transform is correct only when `private$w` is
  binary; §3 shows a live example of it silently producing garbage when it
  isn't.
- **`assert_even_allocation()`** (`design_abstract.R:279-285`) and several
  concrete designs' constructors hard-stop unless `prob_T == 0.5` exactly —
  `DesignFixedGreedy` (`design_fixed_greedy.R:45-47`), `DesignFixedAOptimal`
  (`design_fixed_a_optimal.R:42-44`), `DesignFixedDOptimal`
  (`design_fixed_d_optimal.R:39-41`), `DesignFixedBinaryMatch`
  (`design_fixed_binary_match.R:52-54`). `prob_T` itself is a scalar
  "probability of treatment" — a shape that is inherently binary, not just a
  restrictive default. A `K`-arm design needs a `K`-vector of allocation
  probabilities (or counts) everywhere `prob_T` appears.
- **`Inference$initialize()`** (`inference_all_abstract.R:34-107`) copies
  `private$w = private$des_obj_priv_int$w` verbatim — a snapshot of the
  `Design`'s internal `{0,1}` vector, with no recoding step and no check
  that it actually is binary.
- **Every concrete `Inference*` class's estimator** consumes that `w` via
  the literal pattern `private$w == 1` / `private$w == 0` to split the
  sample into two groups — e.g. `InferenceAllSimpleMeanDiff$compute_estimate`
  (`inference_all_mean_diff.R:81-82`), its bootstrap-weighted variant
  (`inference_all_mean_diff.R:109-110`), its fast bootstrap re-draw loop
  (`inference_all_mean_diff.R:169,186`), and the closed-form SE decomposition
  in `private$shared()` (`inference_all_mean_diff.R:254`, which literally
  builds `X = cbind(1, private$w)` — a single scalar treatment column).
  Likelihood-based classes do the same thing one level up: e.g.
  `ordinal_cond_clogit_shared_multi` builds
  `X_full = cbind(treatment = expanded$w, X_stack)`
  (`inference_ordinal_KK_cond_logit_abstract.R:94`) — again, exactly one
  treatment column.
- **The extension contract itself is scalar-shaped.**
  `InferenceCustomAsymp$fit()` is documented to return "`estimate`: required
  numeric scalar treatment-effect estimate" and "`se`: optional numeric
  scalar standard error" (`extending-edi-r6.md:23-24`). Every one of the
  ~100 concrete `Inference*` classes in `EDI/R/inference_*.R` implements
  this same scalar contract; nothing in the class hierarchy currently
  represents a vector-valued estimand.
- **C++ kernels are binary-shaped too.** `apply_treatment_and_noise_cpp`
  (used by `SimulationFramework`, not by `Inference` directly) applies the
  simulated treatment effect via `w[i] == 1` at six call sites in
  `EDI/src/simulation_dgp.cpp:50,58,68,78,89,106`. Randomization-inference
  kernels like `compute_simple_mean_diff_parallel_cpp` (called from
  `inference_all_mean_diff.R:193`) take a `{0,1}` `w_mat` directly. None of
  these have a `K`-arm counterpart today.

The upshot: "two-arm" isn't a single guard clause to relax. It's a
convention — `{0,1}` internal storage, `{-1,+1}` public API, a scalar
`prob_T`, and a scalar estimand contract — repeated independently at the
`Design` layer, the `Inference` layer, and the C++ kernel layer.

## 2. `DesignFixedFactorial`: the closest existing precedent (and a seam that has since been fixed)

`DesignFixedFactorial` (`EDI/R/design_fixed_factorial.R`) is the one design
in the package built around a `factors` list (e.g. `list(A = 2, B = 2)`)
and `expand.grid()` over factor levels — the right shape, in general, for
`K`-arm treatment when `K` is a combination count. It remains the closest
existing precedent for whatever Phase 1's `K`-arm mechanics end up looking
like.

**Update (2026-08-09): the encoding bug this section originally documented
has been fixed, and the design has been restricted to two arms in the
meantime.** The original version of this section flagged, as live evidence
for §1's claim that the `{0,1}`/`{-1,+1}` convention is genuinely
load-bearing, that `DesignFixedFactorial` stored `private$w` as a raw
combination *index* (`1..num_combinations`) and overrode `get_w()` to
return it unmodified — bypassing the `2L * private$w - 1L` transform every
other design relies on. An empirical repro at the time confirmed this broke
real `Inference*` objects, not just a hypothetical: pairing a 2-level
`DesignFixedFactorial` with `InferenceAllSimpleMeanDiff` (which splits on
`private$w == 1` / `private$w == 0`) silently returned `NA_real_`, since
`w` only ever took values `{1,2}`, never `0`; pairing it with
`InferenceContinOLS` (which regresses on raw `private$w` as a numeric
column) returned a *finite* estimate with the correct magnitude but the
*wrong sign*, since OLS on `{1,2}` estimates "effect of level 1 → level 2,"
the negative of "effect of level 1 (treated) vs. level 0 (control)" under
the rest of the package's convention. Neither failure mode raised an error
or a warning.

That gap was also invisible to `SimulationFramework`'s own compatibility
gate: `.build_valid_combos_for_current_cell()`
(`simulations_framework.R:2731-2797`) checks whether a design/inference
combo is usable by constructing a *throwaway* design with
`skip_assignment = TRUE`, whose validation-only path hardcodes
`priv$w = rep(c(0L, 1L), length.out = n)` directly into the private field
(`simulations_framework.R:3999`) — bypassing `DesignFixedFactorial$get_w()`
entirely. So the gate always saw a well-formed binary `w` and never caught
that the design's *real* `get_w()` (used later, during actual replication)
didn't match. Since `DesignFixedFactorial` is in `SimulationFramework`'s
default design list (`simulations_framework.R:4120`, auto-injected as a
single `factors = list(treatment = 2L)` two-level factor whenever a
design's constructor takes a `factors` argument and the caller didn't
supply one, `simulations_framework.R:2986-2987,3981`) and
`InferenceAllSimpleMeanDiff` is in the default inference list
(`:4137`), this meant an out-of-the-box `SimulationFramework$new(...)$run()`
call was, by default, silently producing `NA`/wrong-signed rows for that
combination.

**The fix applied:** `design_fixed_factorial.R` no longer overrides
`get_w()`, `draw_ws_according_to_design()`, or `assign_w_to_all_subjects()`
— it inherits all three from `Design`/`DesignFixed` unmodified, and its
`draw_ws_raw()` now emits `{0,1}` (not `{1,2}`) so `private$w` follows
the same internal convention as every other design in §1's inventory.
`get_w_factorial()` was updated to convert back to 1-indexed factor levels
(`private$w + 1L`) for its own data-frame output. As a **stopgap** — since
generalizing `Inference*` consumption to genuinely non-binary `w` is
exactly what Phases 1/3 of this document are about, and hadn't been built
yet — the constructor now hard-asserts `prod(unlist(factors)) == 2`
(gated by `should_run_asserts()`, mirroring how `DesignFixedGreedy`/
`DesignFixedAOptimal`/`DesignFixedDOptimal`/`DesignFixedBinaryMatch` already
hard-assert `prob_T == 0.5` for their own unsupported-configuration reasons,
§1), so a `factors` list implying more than two total combinations (e.g.
the original `list(A=2, B=2)` 2×2 example) is now rejected at construction
time instead of silently producing wrong results downstream. Regression
coverage was added directly (`test-fixed-designs-greedy.R`): the two-arm
`{-1,+1}` balance property, the `>2`-combination rejection, and an
end-to-end check that `InferenceAllSimpleMeanDiff`/`InferenceContinOLS`
both recover the correct, non-`NA`, matching-sign estimate against a
`DesignFixedFactorial` design with a known injected effect.

**What this means for the rest of this document:** `DesignFixedFactorial`
is no longer live evidence of a leak — it's now just an ordinary two-arm
design that happens to be built on an `expand.grid`-over-`factors`
representation. It's still the right structural precedent to generalize
*from* in Phase 1 (§3), but it no longer demonstrates what breaks when `w`
escapes the binary convention, and it is not, itself, a source of `K`-arm
capability today — that capability is exactly what §3 still needs to add,
on top of lifting the `prod(unlist(factors)) == 2` assertion this fix just
introduced.

**Forward reference:** the fix above kept `DesignFixedFactorial` on
`Design`'s existing per-`K`-conditional `{-1,+1}` convention (it inherits
`get_w()` unmodified, which still applies `2L * private$w - 1L`). §3c below
proposes going further for the actual `K`-arm work — dropping the `{-1,+1}`
transform from `Design` entirely, for every `K` including `K == 2` — which
would mean this fix is an intermediate state, not the final one: lifting
the `prod(unlist(factors)) == 2` assertion in a future `K`-arm pass would
also mean revisiting `get_w()` on every design, not just `DesignFixedFactorial`.

## 3. Phase 1, design side: which `Design*` subclasses generalize cheaply

Excluding the KK matching family (Phase 2, see §5), EDI's ~26 concrete
design subclasses fall into three buckets by how their assignment mechanism
currently encodes "binary":

### 3a. Free via `randomizr`'s native multi-arm support

Four fixed designs already delegate their core randomization to the
`randomizr` package, whose `*_ra()` functions take `num_arms`/`conditions`/
`prob_each` arguments for multi-arm designs out of the box (verified
against the installed package: `args(randomizr::block_ra)`,
`args(randomizr::cluster_ra)`, and `args(randomizr::complete_ra)` all list
`num_arms`, `conditions`, and `prob_each` alongside the binary `prob`):

| Design | Call site | Current (binary) call |
|---|---|---|
| `DesignFixedBlocking` | `design_fixed_blocking.R:117` | `randomizr::block_ra(blocks = strata_keys, prob = private$prob_T)` |
| `DesignFixedCluster` | `design_fixed_cluster.R:77` | `randomizr::cluster_ra(clusters = cluster_ids, prob = private$prob_T)` |
| `DesignFixedOptimalBlocks` | `design_fixed_optimal_blocks.R:115` | `randomizr::block_ra(blocks = block_ids, prob = private$prob_T)` |
| `DesignFixedBlockedCluster` | `design_fixed_blocked_cluster.R:80` | `randomizr::block_and_cluster_ra(...)` |

For these four, `K`-arm support is close to a one-line change per
`draw_ws_raw()`: swap `prob = private$prob_T` for
`num_arms = private$K` (or `conditions`/`prob_each` if unequal allocation is
wanted), and change `as.numeric(as.character(...))` — which currently
coerces `randomizr`'s `"0"`/`"1"` factor output to a numeric `{0,1}` vector
— to preserve `randomizr`'s `K`-level factor labels instead. No new C++ is
needed for this bucket.

### 3b. Generalizes via a documented K-arm algorithm, but needs new C++

The remaining non-KK designs implement their own biased-coin/urn/
minimization mechanic in C++, and every one of them currently returns/
accepts a scalar `{0,1}` decision:

| Design | Mechanism | Current binary call |
|---|---|---|
| `DesignFixedBernoulli` | i.i.d. Bernoulli(`prob_T`) | `generate_permutations_bernoulli_cpp(n, r, prob_T)` (`design_fixed_bernoulli.R:45-49`) |
| `DesignSeqOneByOneBernoulli` | same, one at a time | analogous |
| `DesignFixediBCRD` / `DesignSeqOneByOneiBCRD` | permuted-block complete randomization | `draw_ws_raw` targeting exact counts (`design_fixed_ibcrd.R:45-50`, `design_seq_one_by_one_ibcrd.R:52-77`) |
| `DesignSeqOneByOneEfron` | Efron's (1971) biased coin | `assign_wt`/`draw_ws_raw` (`design_seq_one_by_one_efron.R:46,65`) |
| `DesignSeqOneByOneAtkinson` | Atkinson's (1982) `D_A`-optimal biased coin | `atkinson_assign_weight_cpp(private$w[1:(t-1)], ...)` (`design_seq_one_by_one_atkinson.R:51-57`) — takes the *history* of past `{0,1}` assignments directly |
| `DesignSeqOneByOneUrn` | Wei's (1978) urn model | `prob_T = (alpha + beta*nC) / (2*alpha + beta*(nT+nC))` then `rbinom` (`design_seq_one_by_one_urn.R:62-64`) — the formula itself is binary-count-shaped |
| `DesignSeqOneByOnePocockSimon` | minimization on stratifying covariates | `assign_wt`/`draw_ws_raw` (`design_seq_one_by_one_pocock_simon.R:68,87`) |
| `DesignSeqOneByOneSPBR`, `DesignSeqOneByOneRandomBlockSize` | (random) permuted block sizes | block-size arithmetic assumes `n_T = round(block_size * prob_T)` (`design_seq_one_by_one_spbr.R:64`) |
| `DesignFixedRerandomization` | rerandomize until covariate balance criterion met | `n_T = round(n * prob_T)` rejection sampler (`design_fixed_rerandomization.R:70,89`) |
| `DesignFixedGreedy` | greedy sequential-swap balance search | `draw_ws_raw` (`design_fixed_greedy.R:68`), hard-requires `prob_T == 0.5` |
| `DesignFixedAOptimal` / `DesignFixedDOptimal` | exchange-search optimal design | `a_optimal_search_cpp(P, H, r, round(n*prob_T))` / `d_optimal_search_cpp(P, r, round(n*prob_T))` (`design_fixed_a_optimal.R:78`, `design_fixed_d_optimal.R:67`), hard-require `prob_T == 0.5` |

Every design in this bucket has an established `K`-arm generalization in
the literature — Wei's biased-coin/urn models, Pocock-Simon minimization,
and `D_A`-optimal biased coins all have published multi-arm variants, and
permuted-block/rerandomization/exchange-search designs generalize to `K`
cells by construction (block randomization and rerandomization don't care
how many labels a block has; an exchange search just needs a `K`-way
allocation-count target instead of `n/2`). The work here is therefore real
but mechanical: each design needs (a) its scalar `prob_T`/binary count
generalized to a `K`-vector, and (b) — for the six designs that call an
`_cpp` kernel with a `{0,1}` decision or history vector baked into its
signature — a new `K`-arm sibling kernel (`generate_permutations_bernoulli_cpp`,
`atkinson_assign_weight_cpp`, `generate_permutations_atkinson_cpp`,
`a_optimal_search_cpp`, `d_optimal_search_cpp`, and the analogous kernels
for Efron/Urn/PocockSimon/SPBR/RandomBlockSize/Rerandomization/Greedy).
Budget C++ work accordingly — this is not purely an R-level change.

### 3c. Shared base-class changes needed regardless of bucket

Whichever designs are tackled first, the shared machinery in `Design`
(`design_abstract.R`), `DesignFixed` (`design_fixed_abstract.R`), and
`DesignSeqOneByOne` (`design_seq_one_by_one_abstract.R`) all need the same
handful of changes, made once and inherited everywhere (mirroring how
`sequential_inference.md` §7 found one shared-base-class change serving 13
subclasses rather than a 13-file change):

- A `K`-arm-aware replacement for `overwrite_all_subject_assignments()`'s
  `{-1,+1}`-only assert (`design_abstract.R:236-241`).
- **`get_w()`/`draw_ws_according_to_design()`/`overwrite_all_subject_assignments()`
  should drop the `{-1,+1}` transform entirely and always speak raw
  `{0,...,K-1}` labels** — for every `K`, including `K == 2`
  (`design_abstract.R:236-241,339-341,346-349`). Concretely this **reverts
  `Design`'s public API to its pre-2026-06-22 form** (see §3c-i below) and
  generalizes it to `K` arms, rather than introducing a new convention.

  This revises an earlier draft of this section twice over. The first draft
  recommended keeping `get_w()`'s `{-1,+1}` behavior for `K == 2` only,
  which keeps two conventions alive side-by-side (a `K == 2` special case
  plus a `K > 2` raw-label case) — exactly the kind of split that caused
  §2's `DesignFixedFactorial` bug. The second draft over-corrected the other
  way: it assumed *no* in-package code actually depended on `get_w()`'s
  signed output, and recommended dropping it everywhere with no
  replacement. That assumption was wrong — investigation (§3c-i) found two
  concrete classes, `InferenceIncidCMH` and `InferenceIncidExtendedRobins`,
  whose standard-error formulas require exactly-signed `{-1,+1}` input, one
  of them enforced by a hard C++-level validation, not just a convenience.
  Dropping the sign without a replacement would have silently broken both.
  The resolution kept here: `Design` speaks one convention, always
  (`{0,...,K-1}`, no exceptions, no `K`-conditional branching) — the
  *consumer*-side need for a signed contrast is met by a small `Inference`-level
  helper used only by the two classes that need it (§3c-i), not by `Design`
  carrying a second convention for their benefit.

  This is still a real, named breaking change to `get_w()`'s current
  (post-2026-06-22) public contract (`design_abstract.R:339`'s roxygen
  documents "A `{-1,+1}` vector... (+1 = treated, -1 = control)"), and
  belongs in a major-version changelog entry, not a silent point-release
  change — but every in-package consumer is either unaffected or actively
  simplified by it; see §3c-i for the complete, verified accounting.

#### 3c-i. Which code actually depends on `get_w()`'s sign, and what replaces it

A precise pass over every caller of `Design$get_w()`/`draw_ws_according_to_design()`
(not just a description of the mechanism) turns up exactly two genuine
dependents, everything else being indifferent or actively wasteful:

- **`InferenceIncidCMH`** (`inference_incidence_cmh.R:113-121`, non-blocking
  branch only — the blocking branch, `compute_cmh_block_se_cpp`, takes only
  `y` and block IDs, no `w` at all, and is unaffected by any of this).
  `get_standard_error()` draws `w_mat` via
  `private$des_obj$draw_ws_according_to_design(private$se_est_num_vectors)`,
  computes `ytw = private$y %*% w_mat`, and uses
  `SE = (2/n) * sqrt(sum(ytw^2) / K)` — valid because `E_w[y'w] = 0`
  *exactly*, for any randomization scheme, only under `{-1,+1}` coding. The
  file's own roxygen (added in the same 2026-06-22 commit that introduced
  `get_w()`'s sign) is explicit that this is "for a balanced design," citing
  Azriel et al. (2026), Eq. 3.
- **`InferenceIncidExtendedRobins`** (`inference_incidence_extended_robins.R:82-87`).
  `get_standard_error()` calls `compute_extended_robins_block_se_cpp(...,
  private$des_obj$get_w(), ...)` — and that C++ kernel *hard-validates* its
  input, `EDI/src/cmh_speedups.cpp:106`: `if (w[i] != -1.0 && w[i] != 1.0)
  return NA_REAL;`. This is a compiled-kernel contract requirement, not a
  convenience.
- **Everyone else is indifferent or actively wasteful.** `simulations_framework.R`'s
  consumers only ever check `w[i] == 1` to identify the treated group
  (`apply_treatment_and_noise_cpp`, `EDI/src/simulation_dgp.cpp:50`) — that
  check is correct regardless of what the *other* value is, so `{0,1}`
  serves it exactly as well as `{-1,+1}`. `inference_all_abstract_rand.R`
  and `inference_all_abstract_rand_bootstrap.R` call
  `draw_ws_according_to_design()` and then immediately convert the result
  back to `{0,1}` themselves (`(w_mat + 1L) / 2L`, 5+ call sites) — for
  these files the current round-trip through `{-1,+1}` is pure overhead
  that this change removes.

**Correctness check on the estimate/SE pairing (not just the encoding):**
both classes inherit `compute_estimate()` from `InferenceAllSimpleMeanDiff`,
which computes `mean(y[w==1]) - mean(y[w==0])` off `private$w` (already
`{0,1}`, untouched by anything on `Design`'s side, before or after this
change). Algebraically, `(2/n) * y'w_signed = mean(y_T) - mean(y_C)`
*exactly* when the design is balanced (`n_T = n_C = n/2`) — precisely the
precondition each class's own documentation/constructor assumes.
`InferenceIncidExtendedRobins` enforces this unconditionally at construction
(`is_blocking_design()` + `get_prob_T() == 0.5` + equal block sizes, always
required, `inference_incidence_extended_robins.R:48-58`), so its estimate
and SE are always describing the same quantity. `InferenceIncidCMH` only
enforces the equivalent checks for its *blocking* branch
(`inference_incidence_cmh.R:72-81`) — its non-blocking,
`draw_ws_according_to_design()`-based branch has no analogous guard, so a
non-blocking design without exact realized balance (e.g. plain Bernoulli
randomization, which doesn't guarantee `n_T = n_C` even at `prob_T = 0.5`)
could see its SE formula's "for a balanced design" precondition violated.
That is a **pre-existing property of today's code**, unrelated to and
unaffected by the `{0,1}`-vs-`{-1,+1}` question — flagged here as a
byproduct of this investigation, not something this document is proposing
to fix.

**Second correctness check: does the true simulated `betaT` itself depend on
the encoding?** No — `apply_treatment_and_noise_cpp` (the `SimulationFramework`
DGP, `EDI/src/simulation_dgp.cpp:46-116`) applies `betaT` via a *conditional*
in every one of its six response-type branches — `const double bt_i = (w[i]
== 1) ? betaT : 0.0;` — never a product `betaT * w[i]`. `betaT` therefore
always means exactly `E[y|treated] - E[y|control]` on the model scale,
independent of what numeric value represents "control." There is no `w`-as-
numeric-regressor use anywhere that a `{0,1}`-vs-`{-1,+1}` switch could
rescale relative to the simulated ground truth.

**Third: is `InferenceIncidExtendedRobins`'s signed-input requirement even
load-bearing mathematically, or only as an enforced contract?** Rereading
`compute_extended_robins_block_se_cpp` in full
(`EDI/src/cmh_speedups.cpp:100-162`): the kernel never uses `w[i]`'s numeric
*value* multiplicatively — `w[i] == 1.0` (line 113) is used purely to route
`y[i]` into a `sum_t`/`sum_c` accumulator, and the returned variance is built
from `p_hat_T = sum_t/n_T`, `p_hat_C = sum_c/n_C` — plain group proportions,
computed exactly the way `mean(y[w==1])`/`mean(y[w==0])` would be under
`{0,1}` coding. Unlike CMH's `y'w` formula (which mathematically *requires*
the zero-mean-under-any-randomization property that only holds for signed
coding), this kernel's `w[i] != -1.0 && w[i] != 1.0 -> NA_REAL` validation
(line 106) is a contract choice, not a mathematical necessity — it would
compute the same correct answer fed `{0,1}` directly. Not a reason to change
the C++ today (the existing contract is tested and correct as-is), but worth
noting as a real future simplification: relaxing that one validation line
would let `InferenceIncidExtendedRobins` drop `get_w_signed()` entirely and
pass `private$w` straight through, leaving `InferenceIncidCMH` as the only
class with a genuine mathematical dependency on the sign.

**The replacement mechanism** is a small private helper on the shared
`Inference` base class, inherited by all ~100 concrete classes but called
by only these two:

```r
# Inference (private):
get_w_signed = function(){
  ifelse(private$w == 1, 1, -1)   # private$w is already {0,1}; identical
                                    # to today's 2*w-1 for w in {0,1}
}
```

This is a pure identity-preserving refactor for `K == 2`: for every subject,
`ifelse(w_raw == 1, 1, -1)` computes exactly the same value today's
`2L * private$des_obj_priv_int$w - 1L` (inside `Design$get_w()`) does — the
sign computation moves, the numbers it produces don't.

One follow-on wrinkle worth naming for whenever `K > 2` designs actually
exist: both classes currently call `private$des_obj$get_w()`/
`draw_ws_according_to_design()` **directly**, bypassing `Inference`'s own
already-arm-selected `private$w` (§4) entirely — CMH needs *fresh*
reference draws from the design's randomization mechanism for its empirical
SE, not just the one observed assignment, which is why it reaches past
`Inference`'s snapshot into `des_obj` directly today. Once `des_obj` can
carry more than two arms, both call sites need to route through `Inference`'s
own (already `arm_treated`/`arm_control`-filtered) copies — or a
`K`-arm-aware sibling that filters `draw_ws_according_to_design()`'s raw
`{0,...,K-1}` draws down to the two selected arms before recoding — rather
than continuing to read the full, unfiltered `K`-arm design directly. This
is a small, mechanical follow-on to §4's arm-selection work, not a new
open question, but it's a change these two specific classes need that most
others don't, since most others only ever consume `Inference`'s own
snapshot.

- `prob_T` (scalar) needs a `K`-vector sibling, e.g. `prob_T_vec`, that
  defaults to `rep(1/K, K)` and collapses to today's `prob_T`/`1-prob_T`
  pair when `K == 2`.
- A new characterization method, e.g. `is_a_multi_arm_capable()`, mirroring
  the existing pattern of `is_a_kk_matching_capable()`
  (`design_abstract.R:43`), `is_a_cluster_capable()` (`:47`), and
  `is_a_bernoulli_capable()` (`:51`) — each of which is a `FALSE`-by-default
  method on `Design` overridden to `TRUE` on the capable subclasses. This
  gives `Inference*` classes and `SimulationFramework` a single dispatch
  point to check design capability, the same way `is_a_kk_matching_capable()`
  already gates KK-only `Inference*` classes (e.g.
  `InferenceAbstractKKWilcoxBaseIVWC$initialize`,
  `inference_all_KK_wilcox_ivwc.R:190-192`, stops unless
  `des_obj$is_a_kk_matching_capable()`).

### 3d. Internal storage representation for `w`: integer, not factor

A separate question from the encoding *values* above: what R type should
`private$w` (and the `n x r` matrices `draw_ws_raw()`/
`draw_ws_according_to_design()` produce) actually be stored as? Today it's
plain `numeric` (double) — `Design`'s private field is declared
`w = numeric()` (`design_abstract.R`), and every concrete design's
`draw_ws_raw()`/`assign_wt()` writes doubles into it. Two alternatives are
worth ruling in/out explicitly, since this storage choice is shared by
every design and every `Inference*` class and therefore multiplies across
the whole package rather than being a local decision.

**Rule out `factor`.** R factors are integer vectors under the hood plus a
`levels`/`class` attribute, which makes them a poor fit for this field for
several concrete reasons rather than a general stylistic preference: (1)
factor levels are 1-indexed by construction, which fights the `{0,...,K-1}`
convention §3c just proposed and every C++ kernel below is written against;
(2) arithmetic on a factor (`2L * w - 1L`-style transforms, `sample()`,
`tabulate()`-style balance bookkeeping used throughout `draw_ws_raw()`
implementations, e.g. `design_fixed_bernoulli.R`,
`design_seq_one_by_one_urn.R:62`) either errors or silently operates on the
underlying integer *codes* rather than the intended values — a well-known R
footgun (`as.integer(factor(2))` is `1`, the position in `levels()`, not
`2`); (3) every C++ boundary crossing (next paragraph) would need an
explicit `as.integer()`/`unclass()` strip before the call regardless,
so factor buys none of integer's coercion savings and adds its own
attribute-handling overhead on top. Factor labels remain the *right* tool
one layer up, for human-readable, presentation-only output — exactly what
`DesignFixedFactorial$get_w_factorial()` already returns today (a
1-indexed data frame of factor-level combinations for display), which is
appropriately kept separate from `private$w`'s internal numeric code.

**Prefer plain `integer` over `numeric` (double) for the internal/hot-path
storage.** This isn't just a style preference — it's grounded in the actual
C++ kernel signatures `w`/`w_mat` cross today. A scan of `EDI/src/*.cpp`
(excluding `RcppExports.cpp`) for `w`/`w_mat` parameter types finds the
large majority declared `IntegerVector`/`IntegerMatrix`, including the
kernels sitting in the hottest paths — inner loops called once per
permutation/bootstrap replicate, not once per `Inference*` construction:
`compute_simple_mean_diff_parallel_cpp`'s `w_mat`
(`simple_mean_diff_parallel.cpp:17`), every `generate_permutations_*_cpp`
kernel's internal `w_mat` (`generate_permutations.cpp`, 9 sites),
`rand_bootstrap_mean_diff_parallel.cpp:57`, and the regression/survival
kernels `fast_robust_regression.cpp:341`, `fast_coxph_regression.cpp:971,1015`,
`fast_weibull_regression.cpp:225`, `fast_survival_stats.cpp:450,526`,
`fast_stereotype_logit.cpp:1157`, `fast_adjacent_category_logit.cpp:426`,
`fast_jonckheere_terpstra.cpp:234`, `fast_gehan_wilcox.cpp:155`,
`fast_logrank.cpp:143,174`, `ridit_distr_parallel.cpp:103,167,236`. A
smaller set of lower-call-frequency kernels (typically invoked once per
design/inference build, not once per replicate) take `NumericVector`/
`NumericMatrix` instead — `cmh_speedups.cpp:69`,
`compute_extended_robins_block_se_cpp` (`bootstrap_match_indices.cpp`-adjacent
files), the `kk21_weights.cpp` response-model weight estimators (5 sites),
`match_stats_from_indices_cpp`, `matching_bootstrap_loop_cpp`,
`base_bootstrap_loop_cpp`.

Since `private$w` is currently `numeric`, every one of the `IntegerVector`/
`IntegerMatrix`-typed calls above already pays an implicit double→integer
coercion at the Rcpp boundary on *every single call* — and because the
`IntegerVector`/`IntegerMatrix` group specifically contains the
permutation/bootstrap inner-loop kernels (called potentially thousands of
times per `Inference*` object, not once), that coercion cost is realized at
exactly the call frequency this package's own performance work
(`package_metadata/audits/perf_experiments_final.md`) has spent the most
effort optimizing around. Storing `private$w` (and the matrices
`draw_ws_raw()` produces) as plain R `integer` instead:
- Eliminates that coercion for the majority (by both count and call
  frequency) of kernel call sites.
- Shifts the coercion burden to the smaller, colder-path `NumericVector`-typed
  group instead — a net win, not a wash, given the frequency skew above.
- Halves the memory footprint per element (4 bytes vs. 8), which is
  negligible for a length-`n` vector but not negligible for the `n x r`
  permutation matrices `draw_ws_raw(r)` builds, where `r` (number of
  permutations/bootstrap draws) commonly runs into the hundreds or
  thousands.
- Costs nothing on the R side: `NA_integer_` is exactly as capable an
  "unassigned" sentinel as `NA_real_` for a field that only ever holds
  small non-negative whole numbers, and `w`'s only cooperating fields that
  need to interoperate element-wise with it (`y`, `dead`) are independently
  double for reasons specific to *them* (continuous responses, etc.), not
  because they need to match `w`'s type.

**Recommendation:** store `private$w` and `draw_ws_raw()`/
`draw_ws_according_to_design()`'s output matrices as `integer`
(`{0L,...,(K-1)L}`), reserve `factor` exclusively for user-facing display
methods analogous to `get_w_factorial()`, and treat this as a change worth
making independently of — and probably before — the `K`-arm work itself,
since it's a real (if modest) performance win for every design shipped
today, not just future `K`-arm ones.

## 4. Phase 1, inference side: choosing two arms in `Inference$initialize()`

The user's proposed Phase-1 inference story — "demand the user decide which
two arms" — has a real advantage: it requires **no changes to any concrete
`Inference*` class**. §3c's revised architecture (`Design` always speaks
raw `{0,...,K-1}` labels, with no per-`K` special case) changes *where*
that selection logic should live relative to an earlier draft of this
section: not as a new `Design$restrict_to_two_arms()` method that builds a
second `Design` object, but as new optional constructor arguments on the
one shared `Inference` base class every concrete class already calls
through — `Inference$initialize()` (`inference_all_abstract.R:34-107`).

**Why the shared base class, not a new `Design` method.**
`Inference$initialize()` already builds its own working copies of
everything it needs from `des_obj` — `private$y`, `private$w`,
`private$dead` as plain value copies (`inference_all_abstract.R:52-55`),
and `private$X` via its own model-matrix construction path
(`:61-91`, either reusing the design's precomputed matrix or building a
fresh one from `model_formula`). Selecting two arms is a row-filter over
exactly those same fields, at exactly the point they're already being
copied — it doesn't need a second `Design` object, a second
imputation/model-matrix pass through
`covariate_impute_if_necessary_and_then_create_model_matrix()`
(`design_abstract.R:510-584`), or a second round-trip through
`add_all_subjects_to_experiment()`/`overwrite_all_subject_assignments()`/
`add_all_subject_responses()`. Concretely:

```r
initialize = function(des_obj, arm_treated = NULL, arm_control = NULL, ...) {
  ...
  w_raw = private$des_obj_priv_int$w   # {0,...,K-1}
  arms_present = sort(unique(w_raw))
  if (is.null(arm_treated) && is.null(arm_control)) {
    if (length(arms_present) != 2L) {
      stop("This design has ", length(arms_present), " arms; ",
           "arm_treated/arm_control must be specified explicitly.")
    }
    arm_control = arms_present[1L]; arm_treated = arms_present[2L]
  }
  keep = w_raw %in% c(arm_control, arm_treated)
  private$w = as.numeric(w_raw[keep] == arm_treated)   # recoded {0,1}
  private$y = private$des_obj_priv_int$y[keep]
  ...                                                   # dead, X rows, etc.
}
```

For every design shipped today (`K == 2`), `arms_present` always has length
2, so `arm_treated`/`arm_control` default correctly and **every existing
call site — internal or user-written — is unaffected**: this is additive,
not a signature break, unlike §3c's `get_w()` change. `K > 2` designs
simply require the caller to say which two arms they mean, which is exactly
the "force a choice" behavior the user described.

**What this buys over the original `Design`-method proposal:**
- No second `Design` build, so no second imputation pass and no
  `assertClass`/`assert_all_responses_recorded()` re-satisfaction concern —
  `Inference$initialize()` already runs those checks against the *original*
  `des_obj` (`inference_all_abstract.R:39,44`), which is never mutated or
  duplicated.
- The blocking/matching-substructure concern from an earlier draft of this
  section gets strictly easier, not harder: `private$m` (block/cluster IDs)
  can still be read off the *intact* original `des_obj` rather than a
  truncated derivative, so a stratified class (CMH etc.) can still see the
  true block structure while contrasting only two of the `K` arms within
  it. The underlying statistical caveat is unchanged, though — a block that
  only ever contained arms other than the chosen pair contributes nothing
  to that particular contrast, and stratified-inference classes need to
  handle that (drop it, or treat its within-block contribution as zero)
  regardless of which architecture implements the selection.
- One shared change (`Inference$initialize()`) rather than a new public
  `Design` method that then needs `should_run_asserts()`-gated validation
  of its own, mirroring how §3c's other base-class changes are also
  concentrated in one place per `sequential_inference.md §7`'s precedent.
- `private$w` lands here already recoded to `{0,1}` (line 544 above), so the
  shared `get_w_signed()` helper (§3c-i) — needed only by
  `InferenceIncidCMH`/`InferenceIncidExtendedRobins` — can be defined once
  on `Inference` and just reads `private$w`, with no knowledge of `K` or of
  which two arms were selected. Those two classes still need their own
  follow-on update to stop reading `des_obj$get_w()`/
  `draw_ws_according_to_design()` directly (§3c-i's closing paragraph),
  since today they bypass this arm-selected snapshot entirely.

**Statistical caveat that Phase 1 should document but not solve:** if a
user constructs multiple `Inference*` objects against the same `K`-arm
`des_obj` with different `arm_treated`/`arm_control` pairs (e.g. every
pairwise comparison among `K = 4` arms, or all `K-1` comparisons against a
reference arm), the resulting p-values/CIs share subjects across
comparisons and standard familywise error rate inflation applies. Phase 1
should not attempt an automatic correction — that's squarely Phase 3a's job
(§6a) — but `arm_treated`/`arm_control`'s documentation should say this
explicitly so it isn't discovered the hard way.

## 5. Phase 2: KK matching-on-the-fly for `K` arms

The KK family is deferred to its own phase because its core mechanism is
**pairwise by construction**, not merely binary by convention — the
distinction that matters is between "stores a `{0,1}` scalar" (mechanical,
§3) and "the matching *algorithm itself* only knows how to produce pairs"
(a genuine new-algorithm problem).

- **`DesignSeqOneByOneKK14`** (`design_seq_one_by_one_KK14.R:61-108`):
  `assign_wt()` either randomizes into a "reservoir" (`private$m[t] = 0`,
  line 64) or finds the nearest already-reservoired subject by Mahalanobis
  distance and assigns the new subject **the opposite arm**:
  ```r
  1 - private$w[reservoir_indices[min_sqd_dist_index]]
  ```
  (`design_seq_one_by_one_KK14.R:96`). `1 - x` is a binary-only operation.
  Generalizing this to `K` arms means replacing "assign the opposite of my
  match" with something like "assign whichever of the `K` arms is furthest
  behind its target allocation among the matched cohort" — which changes
  the *statistical* character of the matching, not just its arithmetic: KK14
  matches exactly one new subject to exactly one reservoir subject and
  balances covariates between exactly two arms per match. A `K`-arm version
  needs either (a) `K`-tuple matching — group `K` mutually-close reservoir
  subjects, one per arm, generalizing non-bipartite matching to a "non-`K`-
  partite" partition problem with no off-the-shelf R equivalent to the
  `nbpMatching` package (see next bullet) — or (b) recursive pairwise
  matching (match into pairs on covariates as today, then randomly split
  each pair's two members across two of the `K` arms), which is easy to
  implement but does not guarantee the "all `K` arms are covariate-balanced
  against each other" property that motivates matching in the first place —
  only guarantees balance between whichever two arms a given pair happened
  to be split across. **This is the primary open research question in this
  document** — it should be resolved (or at least prototyped both ways) via
  simulation before committing to an implementation.
- **`DesignFixedBinaryMatch`** (`design_fixed_binary_match.R:13-174`) is the
  fixed-sample sibling of the same idea: it calls
  `compute_binary_match_structure(X, mahal_match)`
  (`design_fixed_binary_match.R:142`), documented as using "non-bipartite
  matching via `nbpMatching`" (roxygen header, lines 3-6), which by
  construction produces exactly `n/2` pairs — `nbpMatching::nonbimatch()`
  solves a minimum-weight perfect-matching problem on a graph, a
  fundamentally 2-partite formulation. There is no drop-in `K`-partite
  generalization of `nbpMatching` on CRAN; a `K`-arm fixed matched design
  would need either a custom combinatorial-optimization kernel (e.g.
  iterative pairwise matching across arm-pairs, or a relaxation solved via
  integer programming similar to what `DesignFixedDOptimal`/`AOptimal`
  already do via exchange search) or accept the "match into pairs, then
  randomly assign pair-members across arms" approximation from the KK14
  discussion above.
- **`DesignFixedMatchingGreedyPairSwitching`**
  (`design_fixed_matching_greedy_pair_switching.R`) does a greedy swap
  search *within* `DesignFixedBinaryMatch`'s pairs
  (`InferenceNonParamBootstrap`'s validity table describes it as flipping
  "assignments only ever... within binary-match pairs, so every pair has
  exactly one treated subject",
  `inference_all_abstract_non_param_boot.R:85-92`) — it inherits
  `DesignFixedBinaryMatch`'s pairwise dependency directly and would only
  become `K`-arm-relevant once that design is.
- **`DesignSeqOneByOneKK21`** (`design_seq_one_by_one_KK21.R:11`) inherits
  from `DesignSeqOneByOneKK14` and adds response-adaptive matching *weights*
  estimated via `num_boot` bootstrap draws from response models fit on the
  history of `w == 1` vs. `w == 0` subjects so far (`initialize` docs at
  `design_seq_one_by_one_KK21.R:49-71` list per-response-type weight
  estimators — negative binomial for count, beta for proportion, Weibull
  AFT for survival, proportional odds for ordinal — each implicitly fit as
  a two-group contrast). Generalizing KK21 to `K` arms requires **both**
  the KK14 matching-structure generalization above **and** re-deriving each
  of these response-type-specific weight estimators as a `K`-group (ANOVA-
  style) contrast instead of a two-group one. Sequencing matters here:
  don't start KK21 `K`-arm work before KK14 `K`-arm matching is validated,
  since KK21 is strictly additive complexity on top of it.
  `DesignSeqOneByOneKK21stepwise` (`design_seq_one_by_one_KK21_stepwise.R`)
  is a stepwise-weight variant of the same KK21 machinery and inherits the
  same dependency.

Given the research risk concentrated in the matching-structure question
above, Phase 2 should begin with a simulation-only prototype comparing the
`K`-tuple-matching vs. recursive-pairwise-matching approaches on covariate
balance and effective sample size before any package code is written —
analogous to how `sequential_inference.md`'s Phase 0 (§8 of that document)
recommends a zero-package-changes proof of concept before investing in
plumbing.

## 6. Phase 3: multi-arm `Inference`

Two paths, matching the two options the user raised, with a recommended
sequencing.

### 6a. Path (a): orchestrate existing two-arm `Inference*` objects

Build a thin layer above `Design`/`Inference` — analogous in spirit to how
`SimulationFramework` (`simulations_framework.R:237`) already sits above
that pair without either class knowing about replication, and directly
mirroring the `SequentialMonitor` pattern proposed in
`sequential_inference.md §5`. A `MultiArmInference`-style class (name TBD)
would:

1. Accept a `K`-arm `des_obj`, a set of arm-pairs to compare (all
   `choose(K, 2)` pairs, or `K-1` "vs. reference arm" Dunnett-style pairs),
   and an `inference_class` to use for each comparison.
2. For each pair, construct a throwaway
   `inference_class$new(des_obj, arm_treated = a, arm_control = b, ...)`
   directly against the original `K`-arm `des_obj` using Phase 1's
   `arm_treated`/`arm_control` selection (§4) — never storing or reusing the
   `Inference*` object, following exactly the construct-use-discard
   discipline `sequential_inference.md §6` recommends for interim
   `Inference*` objects for the same reason: each one is cheap to
   reconstruct and there is no dependency-tracking mechanism to safely keep
   a stale one around.
3. Collect the resulting per-pair estimates/SEs/p-values and apply a
   multiplicity correction — Bonferroni/Holm at minimum, ideally a proper
   Dunnett-type correction for the "`K-1` vs. common reference" case, since
   that's the structure most trials with `K > 2` arms actually use. This
   correction is new statistical work, not just orchestration plumbing —
   nothing in the codebase currently implements a multi-comparison
   correction.

This path's chief virtue is that it touches **none** of the ~100 existing
`Inference*` classes or their C++ kernels. It's the safe, high-coverage
option and should be built first.

### 6b. Path (b): native `K`-arm `Inference` classes

This is a substantially larger lift than anything else in this document,
because it changes the shape of the `Inference` contract itself, not just
what feeds it. Currently every concrete class's estimand is a numeric
scalar — `extending-edi-r6.md:23-25`'s documented contract
(`estimate`/`se`/`df`, all scalars) is mirrored in every concrete
implementation, e.g. `InferenceAllSimpleMeanDiff`'s
`private$cached_values$beta_hat_T`/`s_beta_hat_T` pair
(`inference_all_mean_diff.R:86,91`). A native `K`-arm class needs either a
`K-1`-vector of contrasts (each arm vs. a reference) plus a `(K-1)x(K-1)`
covariance matrix, or an omnibus test statistic (an `F`-test / likelihood-
ratio test across all `K` arms simultaneously, with no single "the
estimate" to report). Concretely, this would require:

- A parallel base-class contract (`InferenceMultiArm` or similar) rather
  than modifying `Inference`/`InferenceAsymp` in place — retrofitting the
  scalar contract that ~100 classes, the roxygen-documented public API, and
  `extending-edi-r6.md`'s external extension contract all depend on is not
  a safe in-place change.
- `Inference$initialize()`'s `private$w`/`private$X` construction
  (`inference_all_abstract.R:52-91`) would need a `private$W` `K-1`-column
  dummy-matrix sibling.
- New C++ kernels: every kernel enumerated in §1's last bullet
  (`compute_simple_mean_diff_parallel_cpp`,
  `compute_rand_bootstrap_mean_diff_parallel_cpp`, and their equivalents
  across the likelihood-based classes) takes a `{0,1}` vector/matrix by
  signature and would need a `K`-arm-aware sibling, response-type by
  response-type.

**Recommendation:** do not attempt to retrofit the existing `Inference`
hierarchy for this. If native multi-arm inference turns out to be worth the
investment, build a small number of purpose-written classes for the
highest-value response types only (continuous mean-difference/ANOVA first,
then incidence), and leave path (a) as the general-purpose fallback for
every other response type and estimator this package currently supports.
Path (b)'s scope is comparable to writing a second, parallel `Inference`
hierarchy — it should only be started once there's demonstrated user demand
that path (a)'s pairwise-comparison output genuinely doesn't serve, not
speculatively.

## 7. A cross-cutting note on randomization-test validity

The rand-CI/randomization-test family (`inference_all_abstract_rand*.R`)
derives its validity from replaying the *actual* design mechanism via
`draw_ws_according_to_design()` (`design_abstract.R:346-349`) conditional on
the observed covariate/response history — this is the mechanism `compute_fast_randomization_distr`
implementations consume directly (e.g.
`inference_all_mean_diff.R:190-195`, which feeds `permutations$w_mat`
straight into `compute_simple_mean_diff_parallel_cpp`). Once §3's Phase 1
generalizes `draw_ws_according_to_design()` to `K` arms, this validity
argument does *not* automatically carry over for free to the two-arm
selection an `Inference*` object makes via `arm_treated`/`arm_control`
(§4) — the reference distribution for a two-arm contrast drawn from within
a `K`-arm randomization is not the same object as "replay a genuinely
two-arm design's mechanism." Getting this right (conditioning correctly on
"the observed two-arm subset came from a `K`-arm mechanism, restricted to
two arms after the fact") is new statistical work belonging to whichever of
Phase 1/3a it's scoped into — flagged here so it isn't assumed to be free
plumbing once §3/§4 land.

## 8. Phased plan

- **Phase 1a (Design, non-KK):** switch `private$w`'s storage type from
  `numeric` to `integer` (§3d) — an independent, lower-risk change worth
  landing first, since it benefits every design shipped today, not just
  future `K`-arm ones. Then drop `get_w()`/`draw_ws_according_to_design()`/
  `overwrite_all_subject_assignments()`'s `{-1,+1}` transform in favor of
  always speaking raw `{0,...,K-1}` labels (§3c) — a named, documented
  breaking change to `Design`'s public API, landed in a major version, with
  a migration note for any external caller of `get_w()` — together with the
  shared `Inference$get_w_signed()` helper (§3c-i) and the follow-on update
  to `InferenceIncidCMH`/`InferenceIncidExtendedRobins` so they stop reading
  `des_obj$get_w()`/`draw_ws_according_to_design()` directly, since those
  are the only two classes in the package with a genuine (verified, §3c-i)
  dependency on the sign. Then add `K`/`prob_T_vec` to the shared `Design`/
  `DesignFixed`/`DesignSeqOneByOne` base classes (§3c) with `K = 2` as the
  default so every design's *internal* assignment mechanics are unaffected;
  wire `K`-arm support into the four `randomizr`-backed designs first (§3a,
  near-zero new C++); then the remaining non-KK designs one at a time (§3b,
  each needs a new C++ kernel); add `is_a_multi_arm_capable()` mirroring the
  existing `is_a_*_capable()` characterization pattern. (The
  `DesignFixedFactorial` signed-encoding seam identified in §2 has already
  been fixed as an immediate hotfix, ahead of this phased plan, on the old
  per-`K`-conditional `{-1,+1}` convention — it now needs its
  `prod(unlist(factors)) == 2` stopgap assertion lifted once this bullet's
  `get_w()` change and §3's `K`-arm mechanics land.)
- **Phase 1b (Inference, arm selection):** add `arm_treated`/`arm_control`
  parameters to `Inference$initialize()` (§4), defaulting to `NULL` and
  auto-resolving when exactly two arms are present so every existing call
  site is unaffected; ship with documentation and an example showing
  "`K`-arm design → `SomeInference$new(des_obj, arm_treated=.., arm_control=..)`
  → identical behavior to a design that was always two-arm" as the
  supported Phase-1 workflow; document (don't yet solve) the multiplicity
  caveat from §4 and the randomization-test validity gap from §7.
- **Phase 2 (KK matching, K arms):** simulation-only prototype comparing
  `K`-tuple vs. recursive-pairwise matching (§5) before writing package
  code; then `DesignSeqOneByOneKK14`; then `DesignFixedBinaryMatch`/
  `DesignFixedMatchingGreedyPairSwitching`; then `KK21`/`KK21stepwise`'s
  response-adaptive weight generalization, strictly after KK14 lands.
- **Phase 3a (Inference, orchestration):** `MultiArmInference`-style
  pairwise/Dunnett orchestration layer over Phase 1b's `arm_treated`/
  `arm_control` selection (§6a), plus a real multiplicity-correction
  implementation — the first place in this plan that needs genuinely new
  statistical code rather than generalizing an existing mechanism.
- **Phase 3b (Inference, native — demand-gated):** purpose-written
  native `K`-arm `Inference` classes for continuous and incidence responses
  only (§6b), gated on demonstrated need beyond what Phase 3a's pairwise
  output provides.

## 9. Non-goals

- Not proposing to change any *internal* two-arm design's assignment
  mechanics, `Inference*` class behavior, or default `K = 2` — every
  existing design and every existing `Inference*$new(des_obj)` call
  (without `arm_treated`/`arm_control`) must return numerically identical
  estimates/SEs/p-values after this work lands, verified in §3c-i down to
  the exact estimator/SE algebra for the two classes that touch the sign.
  The one deliberate, explicitly-named exception is `Design`'s public
  `get_w()`/`draw_ws_according_to_design()`/`overwrite_all_subject_assignments()`
  convention itself (§3c), which this document *does* propose changing from
  `{-1,+1}` to raw `{0,...,K-1}` labels for every `K`, with the displaced
  `{-1,+1}` need met by a small, narrowly-used `Inference`-level helper
  (§3c-i) rather than by `Design` carrying two conventions. This is a real
  breaking change to `Design`'s public contract, scoped to a major version
  with a migration note, not something to sneak in silently.
- Not selecting between `K`-tuple matching and recursive-pairwise matching
  for Phase 2 (§5) — that is the primary open research question in this
  document and needs simulation evidence, not a document-level decision.
- Not proposing an automatic multiplicity correction as part of Phase 1 —
  `arm_treated`/`arm_control` selection is a manual, user-driven operation
  in Phase 1; automatic correction is scoped to Phase 3a only.
- Not proposing to modify the existing `Inference`/`InferenceAsymp` scalar
  contract, or `extending-edi-r6.md`'s external extension contract, under
  any phase of this plan. Phase 3b's native multi-arm classes are additive
  (a new parallel hierarchy), never a retrofit of the existing one.
- `DesignFixedFactorial`'s `SimulationFramework`-degeneracy encoding bug
  (§2) has already been fixed and is out of scope going forward — the
  remaining `DesignFixedFactorial`-related work is exactly the `K > 2`
  generalization §3/§8 already scope under Phase 1, not a bugfix.
