# Save/Load API for `Design` Objects

Generated: 2026-08-09

Related: [sequential_inference.md](sequential_inference.md) (the
`analysis_log` ledger discussed there is the concrete trigger for this
document — it's the first new private field proposed to be added to
`Design` after users may already have serialized objects in the wild),
[extending-edi-r6.md](extending-edi-r6.md).

## Assessment (why this document exists)

The question that prompted this: EDI has no dedicated save/load API. Is it
sufficient for a "professional" package to implicitly rely on the user
calling base R's `saveRDS()`/`readRDS()` directly on a `Design` object
(persisting the trial), while treating `Inference*` objects purely as
disposable, reconstructible-on-demand wrappers around it (never persisted)?

**Mostly yes — the mental model is correct and already matches how the
code behaves** (see [sequential_inference.md](sequential_inference.md) §3:
`Inference$initialize()` takes value-copy snapshots of `y`/`w`/`X` and
caches derived quantities against them; there is no notion in the codebase
today of an `Inference*` object outliving the `Design` state it was built
from). Making `Design` the single unit of persistence, with `Inference*`
always rebuilt fresh from it, is the right shape.

But "just document it" is not, by itself, sufficient to call this a
professional API. There are two structural gaps that documentation alone
cannot close, both because they only manifest *after* the package has
moved on from the version a given `Design` object was serialized under:

**1. No version/compatibility stamp — and this bites immediately.**
`saveRDS()` on an R6 object serializes the actual method closures as they
existed at save time. If a user saves a `des_obj`, upgrades EDI, and
`readRDS()`s it back, the reloaded object's methods are frozen to the old
package version — and any *new* private field added to the class
definition since the object was saved simply will not exist on it.
`Design`'s private field list
(`EDI/R/design_abstract.R:471-512`) has no version marker today. Since the
package is under active development and this very session already proposed
growing that list (`analysis_log` in
[sequential_inference.md](sequential_inference.md) §7), any method that
reads a field added after some users' objects were already saved needs to
tolerate that field being absent — and there's currently no stamped
signal on the object itself that would let such a method (or the user) know
it's looking at "old-shape" data versus discovering it by accident via a
`NULL`.

**2. Serialization safety hasn't actually been verified, only inferred.**
`private$all_subject_data_cache` and `private$m` (KK-family designs) are
populated from Rcpp calls (`compute_all_subject_data_cpp()`, matching
kernels). Tracing the R-side usage suggests everything landing in those
fields is a plain R matrix/vector, not an external pointer — but that's an
inference from call-site usage, not a verified fact about the C++ return
types. "An Rcpp function secretly returns something holding an `XPtr`" is a
classic silent-corruption failure mode for exactly this kind of save/load
claim (the object serializes without error, then misbehaves or crashes
after reload) — worth a real smoke test before documenting this as
supported, not just documenting it on the strength of code reading.

Given those two gaps, the recommended shape stays minimal and matches the
original proposal — no dedicated `save_edi_design()`/`load_edi_design()`
wrapper functions; that's more surface than this needs — but with a small,
concrete punch list completed *before* the roxygen documentation is written
as a claim of support, not after:

1. Add a version stamp to `Design$initialize()`.
2. Make every newly-added private field self-initializing (lazy-init on
   first use) rather than assuming `initialize()` set it, so objects
   serialized before a given feature shipped don't break when new methods
   are called on them post-upgrade.
3. Actually run a save/reload/continue-enrolling smoke test on at least one
   plain design and one KK-matching design before writing the roxygen
   claim.
4. Only then document the contract explicitly: persist `des_obj` via
   `saveRDS()`; `Inference*` objects are disposable and must not be
   persisted.

## TODOs

Ordered roughly by dependency (later items assume earlier ones are done),
grouped by concern.

### A. Core correctness / forward-compatibility

- [ ] **Add a version stamp to `Design`.** In `Design$initialize()`
  (`EDI/R/design_abstract.R:82`), set a new private field, e.g.
  `private$edi_version_created = as.character(utils::packageVersion("EDI"))`.
  This is the single piece of infrastructure everything else in this
  document depends on: without it, there is no way for code (or a user) to
  tell, after `readRDS()`, whether an object predates a given schema
  change. Cheap, additive, no behavior change for existing code paths.

- [ ] **Add a public accessor** `get_edi_version_created()` (mirroring the
  existing `get_t()`/`get_n()`-style accessor pattern already used
  throughout `Design`, e.g. `EDI/R/design_abstract.R:303`) so users and
  diagnostic tooling can introspect this without reaching into `private`.

- [ ] **Make every *newly added* private field self-initializing.**
  Concretely: any method that reads a private field introduced after the
  version-stamp TODO above must not assume `initialize()` set it. Pattern:
  ```r
  record_analysis_event = function(payload){
    if (is.null(private$analysis_log)) private$analysis_log = list()
    private$analysis_log[[length(private$analysis_log) + 1L]] = ...
  }
  ```
  This is the direct fix for gap #1 above and is what actually makes "old
  objects still work after a package upgrade" true, rather than merely
  documented. Apply this pattern to `analysis_log` /
  `record_analysis_event()` / `get_analysis_events()` from
  [sequential_inference.md](sequential_inference.md) §7 specifically, since
  that is the next field slated to be added.

- [ ] **Decide version-mismatch *behavior*, not just detection.** Once
  objects carry `edi_version_created`, decide what happens when
  `packageVersion("EDI")` at load time differs from the stamped value:
  silently proceed (current implicit behavior, now at least detectable),
  emit a one-time `message()`/`warning()` on first method call after
  detecting a mismatch, or (only for a *major* version bump, if the package
  ever declares a breaking schema change) refuse with a clear error naming
  the mismatch. Recommend: warn on major-version mismatch only, silent
  otherwise — most field additions will be additive/backward-compatible
  given `lock_objects = FALSE`, so erroring by default would be overly
  aggressive for routine minor/patch upgrades.

- [ ] **Audit `Design`'s own private fields for serialization safety.**
  Confirm none of `private$all_subject_data_cache`, `private$permutations_cache`,
  `private$lin_centered_covariates` (`EDI/R/design_abstract.R:471-512`) can
  ever hold a non-serializable value (external pointer, connection, socket,
  environment that isn't meant to be duplicated). These are populated from
  C++ kernels (`compute_all_subject_data_cpp()` and friends) — trace the
  actual `Rcpp::` return types, don't just infer from call-site usage.

- [ ] **Extend the audit to `DesignMatching`/`DesignSeqOneByOne` subclass
  private fields**, not just base `Design`. In particular
  `private$m` (match/reservoir index, `DesignSeqOneByOneKK14`:
  `EDI/R/design_seq_one_by_one_KK14.R`), and any per-subclass caches in
  Atkinson/SPBR/PocockSimon (e.g. `strata_cols`-keyed caches). Same failure
  mode as above, different files — this needs to be per-subclass because
  each concrete `DesignSeqOneByOne*` class adds its own private state on
  top of the shared abstract class.

- [ ] **Document (and verify) RNG/reproducibility semantics across a
  save/reload cycle.** `private$seed` is only consumed once, inside
  `maybe_set_seed()` at construction time
  (`EDI/R/design_abstract.R:638`) — it is *not* re-applied on `readRDS()`.
  This means continuing to enroll subjects after a reload draws from
  whatever the *global* `.Random.seed` happens to be in the new session,
  not a deterministic continuation of the original stream. That's almost
  certainly the right behavior for a real trial (you don't want
  determinism in production), but it should be stated explicitly in the
  roxygen so nobody expects bit-for-bit-reproducible continuation across a
  restart, and it's worth double-checking it doesn't silently break any
  *test* that relies on reload+continue determinism.

- [ ] **Check `Design$duplicate()` against the same concerns.**
  `duplicate()` (`EDI/R/design_abstract.R:441`) already has interesting,
  possibly-relevant precedent: it explicitly clears
  `d$.__enclos_env__$private$seed = NULL` on the copy. Confirm whether
  `readRDS()`-based reload should follow the same convention (leave `seed`
  as-is, since it's inert after construction — see previous TODO — or
  proactively clear it for consistency with `duplicate()`'s existing
  choice). Also confirm `clone()` (which `duplicate()` wraps) correctly
  copies newly-added fields like `analysis_log` without special-casing —
  R6's default shallow `clone()` copies all bindings present in `private`
  at clone time, so this should be automatic, but verify once the new
  field exists.

### B. Testing

- [ ] **Round-trip smoke test (plain design).** Build a
  `DesignSeqOneByOneBernoulli`, enroll several subjects with responses,
  `saveRDS()`/`readRDS()`, continue enrolling, run an `Inference*` class
  against it, and assert identical structure/estimates to an unsaved
  control run with the same seed/data. This is the direct verification for
  gap #2 in the Assessment above — do this *before* writing the roxygen
  claim, not after.

- [ ] **Round-trip smoke test (KK-matching design).** Same as above but for
  `DesignSeqOneByOneKK14`, specifically exercising `private$m` and
  `private$all_subject_data_cache` mid-trial (i.e. save/reload while the
  reservoir/matching state is non-trivial, not just at `t = 0`).

- [ ] **Version-mismatch behavior test.** Construct an object, manually
  strip/rewrite its stamped `edi_version_created` (simulating an
  object saved by an older package version, including the "field doesn't
  exist at all" case for objects saved before the version stamp itself was
  added), and assert the chosen behavior from the "decide version-mismatch
  behavior" TODO above actually fires (warns / stays silent / errors as
  specified) rather than crashing on a missing field.

- [ ] **Cross-subclass field-audit test.** A single parameterized test
  iterated over all 13 `DesignSeqOneByOne*` subclasses (matching the list
  already enumerated in `SimulationFramework..default_design_classes`,
  `EDI/R/simulations_framework.R:4105`) that does a save/reload round trip
  and asserts `identical()` (or an appropriate near-equality check for any
  floating-point cache) between pre-save and post-reload state.

### C. Documentation

- [ ] **Add a `@section Saving and loading:` block to `Design`'s roxygen**
  (and/or `DesignSeqOneByOne` specifically, since that's the primary
  multi-session use case) stating the supported contract explicitly:
  persist the `Design` object via `saveRDS()`/`readRDS()`; `Inference*`
  objects are disposable, cheaply reconstructed from the design, and must
  never be persisted directly (reconstruct them fresh each session).

- [ ] **Include a worked example** in that section: save mid-trial, reload
  in a fresh session, continue enrolling, run inference — mirroring the
  round-trip tests in section B so the documented example is literally the
  tested path, not just an untested illustration.

- [ ] **State the RNG/reproducibility caveat explicitly** (see the A
  section TODO above) in the same roxygen block, so it's discoverable at
  the point a user would look for save/load guidance, not buried in an
  internal design doc.

- [ ] **State explicitly that `Inference*` objects must not be
  `saveRDS()`'d.** Nothing currently prevents a user from trying it (no
  guard exists), and it would "work" mechanically (R6 objects serialize
  regardless of whether that's a sane thing to do) while producing a frozen
  snapshot the user might mistakenly believe stays live against the
  design. Cheapest fix here is documentation-only; a runtime guard is
  probably not worth the complexity unless this turns out to be a common
  support question in practice — noted as a fallback, not a committed TODO.

### D. Open decisions (need a call before implementation, not just research)

- [ ] **`saveRDS()` `version=` argument: recommend one, or leave to the
  user?** The only existing internal precedent in the package,
  `SimulationFramework`'s replication cache, explicitly pins
  `saveRDS(cache_record, tmp, version = 2)`
  (`EDI/R/simulations_framework.R:2446,2479`). Decide whether the new
  `Design` save/load documentation should recommend the same explicit
  `version = 2` (for cross-R-version portability) for consistency with that
  existing internal convention, or leave it to R's current default and not
  editorialize. Leaning toward recommending consistency with the existing
  internal precedent, but flagging as a decision rather than assuming it.

- [ ] **Confirm the "no wrapper functions" decision stands** after the
  audit/testing items above are done, rather than being assumed up front.
  If the serialization audit (A) turns up something that genuinely needs
  wrapping (e.g. a field that must be stripped/rebuilt on load, not just
  lazily defaulted), a thin `load_edi_design()` that does that
  transformation may become justified even though it isn't today. Revisit
  this decision after, not instead of, doing the audit — don't build the
  wrapper speculatively.

## Non-goals

- Not building a general-purpose object versioning/migration framework for
  the whole package — this is scoped to `Design`/`DesignSeqOneByOne`
  persistence specifically.
- Not proposing that `Inference*` objects gain any persistence support —
  the recommendation is the opposite: make their disposability explicit and
  documented, not add save/load capability to them.
- Not committing to the `save_edi_design()`/`load_edi_design()` wrapper
  functions in section D — that's an open decision contingent on what the
  audit in section A finds, not a planned deliverable.
