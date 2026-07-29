# SEXP Removal / RcppEigen-Only Core Conversion Spec

Generated: 2026-07-28

## Scope

This spec defines an implementation plan for removing `SEXP` (and other
R-specific `Rcpp::` types) from the *numeric core* of `EDI/src`, so that every
kernel's actual computation is expressed purely in terms of `Eigen` types and
plain C++, with R-glue confined to a thin, mechanically generated boundary
layer. This is prerequisite groundwork for exposing the same numeric cores to
Python (via pybind11) without duplicating logic, and is independently useful
because it makes the C++ cores unit-testable outside of R.

This spec does **not** implement Python bindings. It only produces
Rcpp-independent cores plus thin R adapters. See the earlier conversation
in this session for the follow-up pybind11 spec.

Related documents:

- [perf_experiments_final.md](perf_experiments_final.md) — do not duplicate
  perf work tracked there; this spec is structural, not a performance pass.

## Motivation

An audit of `EDI/src` (2026-07-28) found:

- 29 exported functions across 13 files declared with `SEXP` as the return
  type (`fast_survival_models_optim.cpp`, `fast_poisson_glmm.cpp`,
  `fast_cpoisson_combined.cpp`, `fast_clogit_plus_glmm.cpp`,
  `fast_zero_one_inflated_beta.cpp`, `fast_stereotype_logit.cpp`,
  `fast_ordinal_regression.cpp`, `fast_ordinal_probit_regression.cpp`,
  `fast_ordinal_cloglog_regression.cpp`, `fast_ordinal_cauchit_regression.cpp`,
  `fast_coxph_regression.cpp`, `sample_mode.cpp`, `fast_logrank.cpp`). In every
  case checked, the `SEXP` return type is cosmetic — the function body only
  ever returns `List::create(...)` or `wrap(...)`, both of which convert
  implicitly to `SEXP`. `sample_mode.cpp` is the one exception: it switches on
  `TYPEOF(data)` to dispatch across R vector types, which is a genuine
  R-boundary concern (see Non-Goals).
- 54 files call `List::create(...)` to build a named return list — this is
  the dominant SEXP-touching pattern, and in every sampled case it appears
  only in the last few lines of a function whose body up to that point is
  pure `Eigen::VectorXd`/`MatrixXd` arithmetic.
- 51 files accept `SEXP ..._sexp` parameters purely to build a zero-copy
  `Eigen::Map` over the R object (pattern: `SEXP X_sexp` ->
  `Rcpp::NumericMatrix X_r(X_sexp)` -> `Eigen::Map<const Eigen::MatrixXd>
  X(X_r.begin(), ...)`). RcppEigen already performs this exact conversion
  automatically when a parameter is declared as `Eigen::Map<const
  Eigen::MatrixXd>` directly — the manual `_sexp` + `NumericMatrix` + `Map`
  dance is unnecessary boilerplate, not a required pattern.
- 15 files depend on R's RNG (`unif_rand`, `GetRNGstate`/`PutRNGstate`,
  `R::runif`, etc.): `bootstrap_indices.cpp`, `bootstrap_match_indices.cpp`,
  `generate_permutations.cpp`, `pocock_simon_assign.cpp`,
  `design_fixed_greedy.cpp`, `efron_redraw.cpp`, `kk14_redraw.cpp`,
  `rerandomization_helpers.cpp`, `binary_match_search.cpp`,
  `atkinson_redraw_batch.cpp`, `atkinson_assign.cpp`, `fast_sample_int.cpp`,
  `simulation_dgp.cpp`, `sample_bootstrap_distr_weighted_distances.cpp`. This
  is the one place with a genuine design decision (see RNG Handling below),
  not a mechanical rewrite.
- Zero raw R C API calls (`Rf_allocVector`, `PROTECT`/`UNPROTECT`,
  `SET_VECTOR_ELT`) and zero callbacks into R (`Rcpp::Function`,
  `Rcpp::Environment`, `Rf_eval`) anywhere in `EDI/src`. This rules out the
  hardest class of R-coupling; everything remaining is either `Rcpp::` sugar
  types or RNG.

Net picture: the refactor is large in file count but low in algorithmic risk.
The numeric logic does not need to change; only the boundary code does.

## Non-Goals

- Do not change any numeric algorithm, tolerance, or default argument.
  Behavior must be bit-identical before/after for every migrated function
  (see Testing).
- Do not touch `perf_experiments_final.md`-tracked optimizations; if a
  migrated file happens to also need a perf fix, track that separately.
- Do not attempt to make RNG-dependent files produce cross-language-identical
  streams unless explicitly decided (see RNG Handling) — that is a distinct,
  larger effort (reimplementing R's Mersenne-Twister-with-inversion
  generator) and is out of scope here.
- Do not migrate `sample_mode.cpp`'s `TYPEOF`-dispatch logic into the
  Rcpp-free core. That function's job is literally "dispatch on R vector
  type," which is an R-boundary concern by definition — keep it entirely
  in the R-facing wrapper layer, unconverted.
- Do not build the pybind11 module in this pass. This spec stops at
  "Rcpp-free cores + thin Rcpp wrappers that produce byte-identical R
  output." Exposing the cores to Python is a follow-up spec.

## Target Architecture

### Canonical Result Container (replaces per-function `List::create`)

Rather than defining a bespoke plain-C++ struct per migrated function (which
would mean ~50 one-off types), use one generic ordered key/value container in
the Rcpp-free core, plus exactly one converter per target language. This
keeps the core's return-building code nearly line-for-line identical to the
current `List::create(Named("x") = a, Named("y") = b)` calls, so the mechanical
diff per function is small and easy to review.

Add `EDI/src/result_map.h` (no `Rcpp` or `pybind11` include — this is the
Rcpp-free core boundary):

```cpp
#ifndef EDI_RESULT_MAP_H
#define EDI_RESULT_MAP_H

#include <RcppEigen.h>  // for Eigen::VectorXd/MatrixXd only — not for SEXP use
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace edi {

// std::monostate represents R_NilValue / Python None (typed "not available").
using ResultValue = std::variant<
	std::monostate,
	bool,
	int,
	double,
	std::string,
	Eigen::VectorXd,
	Eigen::MatrixXd
>;

class ResultMap {
public:
	ResultMap& set(std::string name, ResultValue value) {
		entries_.emplace_back(std::move(name), std::move(value));
		return *this;
	}

	const std::vector<std::pair<std::string, ResultValue>>& entries() const {
		return entries_;
	}

private:
	std::vector<std::pair<std::string, ResultValue>> entries_;
};

} // namespace edi

#endif
```

Add `EDI/src/result_map_rcpp.h` (the only place that includes both
`result_map.h` and `Rcpp.h` — this header is never included by a core `.cpp`
file, only by the thin `// [[Rcpp::export]]` wrapper functions):

```cpp
#ifndef EDI_RESULT_MAP_RCPP_H
#define EDI_RESULT_MAP_RCPP_H

#include "result_map.h"
#include <Rcpp.h>

namespace edi {

inline Rcpp::List to_rcpp_list(const ResultMap& m) {
	Rcpp::List out;
	Rcpp::CharacterVector names;
	for (const auto& [name, value] : m.entries()) {
		out.push_back(std::visit([](auto&& v) -> SEXP {
			using T = std::decay_t<decltype(v)>;
			if constexpr (std::is_same_v<T, std::monostate>) {
				return R_NilValue;
			} else {
				return Rcpp::wrap(v);
			}
		}, value));
		names.push_back(name);
	}
	out.names() = names;
	return out;
}

} // namespace edi

#endif
```

A parallel `result_map_pybind.h` (pybind11 `py::dict` converter, same
`std::visit` shape) is added in the follow-up pybind11 spec, not here — but
because the core already returns `edi::ResultMap`, that follow-up file is
the *only* new code that spec needs to write per function; no core changes.

### Migration Pattern Per Function

For a function currently shaped like:

```cpp
// [[Rcpp::export]]
SEXP fast_poisson_glmm_cpp(SEXP X_sexp, SEXP y_sexp, /* ... */) {
	Rcpp::NumericMatrix X_r(X_sexp);
	Eigen::Map<const Eigen::MatrixXd> X(X_r.begin(), X_r.nrow(), X_r.ncol());
	// ... pure Eigen computation ...
	return List::create(
		Named("b") = par.head(p),
		Named("converged") = converged
	);
}
```

migrate to:

```cpp
// fast_poisson_glmm_core.h / .cpp — no Rcpp include
edi::ResultMap fast_poisson_glmm_core(
	const Eigen::Map<const Eigen::MatrixXd>& X,
	const Eigen::Map<const Eigen::VectorXd>& y,
	/* ... */
) {
	// ... identical Eigen computation, unchanged ...
	return edi::ResultMap{}
		.set("b", par.head(p))
		.set("converged", converged);
}
```

```cpp
// fast_poisson_glmm.cpp — thin Rcpp wrapper, kept in the R package build only
#include "fast_poisson_glmm_core.h"
#include "result_map_rcpp.h"

// [[Rcpp::export]]
Rcpp::List fast_poisson_glmm_cpp(
	const Eigen::Map<const Eigen::MatrixXd>& X,
	const Eigen::Map<const Eigen::VectorXd>& y,
	/* ... */
) {
	return edi::to_rcpp_list(fast_poisson_glmm_core(X, y, /* ... */));
}
```

Notes:

- Declaring the parameter directly as `Eigen::Map<const Eigen::MatrixXd>`
  lets RcppEigen's `Rcpp::as<>` perform the zero-copy conversion at the
  export boundary automatically — this deletes the manual `_sexp` +
  `NumericMatrix`/`NumericVector` + `Map` boilerplate in the 51 files that
  currently do it by hand, it does not just relocate it.
- The wrapper function is intentionally still declared with
  `// [[Rcpp::export]]` and still lives in the same `.cpp` file as before
  (or a small `_core.cpp` split — see File Structure) so `Rscript
  fast_roxygenize.R` continues to regenerate `RcppExports.cpp`/`RcppExports.R`
  unchanged in signature.
- `NA_REAL`/`NA_INTEGER` scalars stored as `double`/`int` payload values pass
  through `Rcpp::wrap` unchanged, since they're just bit patterns — no special
  case needed in `to_rcpp_list`.

### File Structure

Two layout options, in order of preference:

1. **Same file, split top/bottom** (preferred for smaller files — most of the
   50 files in scope): keep `fast_poisson_glmm.cpp` as one file, with the
   Rcpp-free core function(s) at the top (no `Rcpp.h` include needed there
   because `RcppEigen.h` already pulls in the Eigen types without requiring
   any `SEXP`/`List` usage) and the thin `// [[Rcpp::export]]` wrappers at
   the bottom. This avoids a wave of new files for simple cases.
2. **`_core.cpp` split** (use only for files that are large — e.g.
   `fast_survival_models_optim.cpp`, `fast_clogit_plus_glmm.cpp`,
   `fast_cpoisson_combined.cpp`, all with 3+ SEXP-returning entry points and
   substantial shared internal state/objective classes): extract to
   `fast_survival_models_optim_core.h`/`.cpp` + keep
   `fast_survival_models_optim.cpp` as the thin wrapper file. Do this only
   where the existing file has grown large enough that the split materially
   helps readability — don't split trivially small files just for
   consistency.

Either way, the invariant is: no `.cpp`/`.h` file that is part of a "core" is
allowed to `#include <Rcpp.h>` or reference `SEXP`, `Rcpp::List`,
`Rcpp::NumericVector`, etc. This is mechanically checkable (see Tests).

## RNG Handling

The 15 RNG-dependent files are not part of the mechanical SEXP-removal pass
in this spec. They need an explicit decision because R's RNG stream
(`GetRNGstate`/`unif_rand`/`PutRNGstate`) is what makes `set.seed()`
reproducible today, and swapping to `std::mt19937` changes that.

Decision for this spec: **do not migrate RNG-dependent files in Phase 1–4
below.** Leave them SEXP-touching and R-coupled for now. Track the RNG
question as a separate follow-up decision (own spec) with two options to
choose between at that time:

- keep R's RNG as the source of randomness (call `unif_rand()` from both the
  R build and, later, a small embedded-R-free reimplementation of R's
  Mersenne-Twister + inversion normal generator for the Python build), which
  preserves cross-language seed reproducibility at implementation cost, or
- accept that Python-side draws use a different generator (e.g.
  `std::mt19937` seeded independently), document the divergence, and give up
  cross-language bit-identical reproducibility for resampling-based
  inference.

This spec's Phase list therefore only covers the ~50 non-RNG files.

## Error Handling

Replace the 19 `Rcpp::stop(...)` call sites in migrated core functions with
a plain `std::runtime_error` (or a small `edi::CoreError` subclass if a
distinguishable exception type proves useful during migration — decide when
the first case is hit, don't pre-build it speculatively). The thin Rcpp
wrapper does not need to catch and re-wrap: Rcpp already translates
`std::exception` into an R condition at the export boundary, matching
current `Rcpp::stop` behavior.

## Testing

For every migrated function, migration must be a no-op from R's point of
view. Verify with:

1. Build the package before migration
   (`R CMD INSTALL --no-docs EDI`) and capture output of the existing test
   file that exercises the function (from `package_tests/` — identify via
   `grep -rl <function_name> package_tests/`) as a baseline.
2. Migrate the function per the pattern above.
3. Run `Rscript fast_roxygenize.R` to regenerate `RcppExports.cpp`/`.R`
   ([[feedback_roxygenize]]).
4. Rebuild and re-run the same test file; diff against the baseline. Outputs
   must be bit-identical (not just "close") since no numeric code changed.
5. Add a compile-time check that core files stay Rcpp-free: a small script
   (`scripts/check_core_no_rcpp.sh`) that greps every file matching
   `*_core.cpp`/`*_core.h` plus every file with a "core function" region for
   `#include <Rcpp` or `\bSEXP\b` and fails if found. Wire this into whatever
   CI/pre-commit check already runs `fast_roxygenize.R` if one exists;
   otherwise leave it as a standalone script invoked manually until a CI
   pass is added.

Do not skip step 4's bit-identical diff even for functions that look like
pure refactors — this is the entire safety net for a change that touches
~50 files without altering intended behavior.

## Implementation Phases

### Phase 1: Infrastructure

- Add `EDI/src/result_map.h` and `EDI/src/result_map_rcpp.h` as specified
  above.
- Migrate exactly one representative function end-to-end
  (`fast_poisson_glmm_cpp` in `fast_poisson_glmm.cpp`, chosen because it was
  already inspected in this spec's audit) to validate the pattern, including
  the bit-identical test-diff step.
- Add `scripts/check_core_no_rcpp.sh`.

### Phase 2: Zero-copy Input Cleanup (51 files)

- For every file using the manual `SEXP ..._sexp` -> `Rcpp::NumericMatrix`/
  `NumericVector` -> `Eigen::Map` pattern, change the parameter type directly
  to `const Eigen::Map<const Eigen::MatrixXd>&` / `const Eigen::Map<const
  Eigen::VectorXd>&` and delete the manual conversion lines. This phase does
  not require the `ResultMap` change yet and can proceed independently/in
  parallel with Phase 3 — do it first since it is lower-risk and mechanically
  identical across all 51 files.
- Exclude the 15 RNG files (Non-Goals) and `sample_mode.cpp` (needs its
  `SEXP`/`TYPEOF` dispatch, Non-Goals).

### Phase 3: Output Migration — `List::create` Files (54 files)

- For each file, split the tail `List::create(...)` into an `edi::ResultMap`
  build, per the Migration Pattern.
- Batch by response-type family for reviewability (matches the existing
  `Collate:` grouping in `EDI/DESCRIPTION`): continuous/OLS family first
  (smallest, already partly audited), then count/ordinal, then
  incidence/survival/proportion, then the KK-combined/GLMM files last (they
  tend to be the largest and most state-heavy, e.g.
  `fast_clogit_plus_glmm.cpp`, `fast_cpoisson_combined.cpp`).
- Apply the File Structure decision (same-file split vs `_core` split) per
  file as encountered, not all upfront.

### Phase 4: Output Migration — Raw `SEXP`-Return Files (13 files, 29 functions)

- Same pattern as Phase 3; these are functionally identical to Phase 3 cases
  (return type is cosmetic — see Motivation) except the declared return type
  changes from `SEXP` to `Rcpp::List` on the wrapper.
- `sample_mode.cpp` is excluded (Non-Goals) — leave entirely as-is.

### Phase 5: Verification Sweep

- Run `scripts/check_core_no_rcpp.sh` across all of `EDI/src` and fix any
  remaining leaks.
- Run the full `package_tests/` suite once (not per-file) to catch any
  cross-file interaction missed by the per-function bit-identical diffs.
- Update `[[feedback_roxygenize]]`-style docs if any exported function
  signature's argument order changed (it should not have — flag if it did).

## TODO Checklist

- [ ] TODO-1: Add `EDI/src/result_map.h` (`edi::ResultValue`, `edi::ResultMap`).
- [ ] TODO-2: Add `EDI/src/result_map_rcpp.h` (`edi::to_rcpp_list`).
- [ ] TODO-3: Add `scripts/check_core_no_rcpp.sh`.
- [ ] TODO-4: Migrate `fast_poisson_glmm_cpp` end-to-end as the pilot; confirm
      bit-identical test diff before proceeding further.
- [ ] TODO-5: Phase 2 — zero-copy input cleanup across the 51 `_sexp`-pattern
      files (excludes the 15 RNG files and `sample_mode.cpp`).
- [ ] TODO-6: Phase 3 — migrate continuous/OLS-family `List::create` files.
- [ ] TODO-7: Phase 3 — migrate count/ordinal-family `List::create` files.
- [ ] TODO-8: Phase 3 — migrate incidence/survival/proportion-family
      `List::create` files.
- [ ] TODO-9: Phase 3 — migrate KK-combined/GLMM `List::create` files
      (`fast_clogit_plus_glmm.cpp`, `fast_cpoisson_combined.cpp`, and
      similarly state-heavy files).
- [ ] TODO-10: Phase 4 — migrate the remaining raw-`SEXP`-return files:
      `fast_survival_models_optim.cpp`, `fast_zero_one_inflated_beta.cpp`,
      `fast_stereotype_logit.cpp`, `fast_ordinal_regression.cpp`,
      `fast_ordinal_probit_regression.cpp`,
      `fast_ordinal_cloglog_regression.cpp`,
      `fast_ordinal_cauchit_regression.cpp`, `fast_coxph_regression.cpp`,
      `fast_logrank.cpp`.
- [ ] TODO-11: Phase 5 — full-repo `check_core_no_rcpp.sh` sweep +
      full `package_tests/` run.
- [ ] TODO-12 (separate follow-up spec, not this one): decide RNG-stream
      strategy for the 15 excluded files, then migrate them.
- [ ] TODO-13 (separate follow-up spec, not this one): write the pybind11
      binding layer against the now-Rcpp-free cores, adding
      `result_map_pybind.h`.

## Acceptance Criteria

The feature is complete when:

- No file under `EDI/src` that is part of a migrated function's core
  contains `#include <Rcpp` or the token `SEXP`, except the excluded RNG
  files and `sample_mode.cpp`.
- Every migrated function's R-visible output is bit-identical to its
  pre-migration output on the existing `package_tests/` suite.
- `RcppExports.cpp`/`RcppExports.R` regenerate via `Rscript fast_roxygenize.R`
  with unchanged exported function signatures (name, argument order, argument
  types as seen from R).
- `scripts/check_core_no_rcpp.sh` passes across `EDI/src`.
- `edi::ResultMap`/`edi::ResultValue` is the only return-building mechanism
  used by migrated core functions — no new bespoke per-function result
  structs were introduced.
