# Plan: Professionally Comprehensive Public API Test Suite

## Goal

Make the package test suite substantially more comprehensive by testing more legal combinations of public API arguments, without turning the suite into an unbounded runtime sink or coupling it to private implementation details.

The testing target is the public behavior of the package:

- Exported R6 classes, especially design and inference classes.
- Public R6 constructors and public methods.
- Exported functions listed in `EDI/NAMESPACE`.
- Documented package workflows represented in `EDI/man/*.Rd`.
- Package-level test runners under `package_tests/`.

Internal functions may be exercised indirectly through public calls. Direct tests of private methods, non-exported objects, or `EDI:::` implementation kernels should not be part of this plan unless they are promoted to a supported public surface. Test helpers are acceptable when they generate legal fixtures, instantiate public objects, run public methods, or summarize coverage.

## Current State

The existing suite already has useful breadth:

- `package_tests/comprehensive_tests.R` sweeps response type, design type, dataset, treatment effect, inference class, model formula, and many inference method families.
- `run_tests_for_response()` constructs many legal design and inference combinations and records per-call outcomes.
- `run_inference_checks_impl()` already exercises many public inference methods, including asymptotic, bootstrap, Bayesian bootstrap, parametric bootstrap, Bartlett, jackknife, randomization, randomization-bootstrap, and custom-randomization paths.
- `run_exhaustive_remaining_inference_classes()` discovers exported `Inference*` R6 generators and attempts legal construction through the public API.
- `EDI/tests/testthat/test-argument-permutations.R` contains an early argument-permutation idea, but it mixes public behavior with direct implementation-level calls.

The main gap is not just more loops. The suite needs an explicit legal-argument model: for each public API method, define valid values, constraints between arguments, expected support by class/response/design, runtime tier, and required invariants.

## Principles

1. Test public contracts first.

   A failure should tell us that a user-visible behavior regressed. Avoid direct tests of private fields, private methods, unexported helpers, or C++ implementation functions unless they are exported and documented as supported user APIs.

2. Generate legal combinations from metadata.

   Do not hand-write every permutation. Maintain registries describing legal argument values and constraints, then generate a bounded set of combinations.

3. Use pairwise coverage by default, targeted higher-order coverage where risk demands it.

   Exhaustive Cartesian products will become impractical. Pairwise combinations catch most argument interaction defects at a fraction of the cost. Use 3-way or exhaustive sweeps only for small, high-risk surfaces.

4. Separate "does it run" from "is it statistically correct".

   Combination tests should verify that legal arguments are accepted, return structurally valid outputs, respect documented shapes/ranges, and preserve deterministic behavior under fixed seeds. Statistical calibration and numerical accuracy should live in focused tests with tighter fixtures.

5. Make coverage measurable.

   Every generated call should produce a coverage row recording the public class/function, method, argument set, support decision, status, runtime, output class/shape, and failure reason.

6. Keep runtime tiers explicit.

   A professional suite needs fast CI tests, broader nightly tests, and expensive release/audit tests. Legal argument coverage should be selectable by tier.

## Scope

### In Scope

- Public design constructors:
  - `DesignSeqOneByOne*`
  - `DesignFixed*`
  - documented required arguments such as `n`, `response_type`, `design_formula`, `strata_cols`, `cluster_col`, block sizes, matching options, optimization/search options, and assignment parameters.

- Public design methods:
  - adding subjects
  - assigning treatments
  - recording responses
  - retrieving assignments, covariates, responses, and design summaries where public methods exist.

- Public inference constructors:
  - all exported `Inference*` classes that are concrete user-facing classes.
  - constructor arguments such as `design`, `model_formula`, `verbose`, and family-specific documented options.

- Public inference methods:
  - estimate methods.
  - asymptotic p-values and confidence intervals.
  - bootstrap, Bayesian bootstrap, parametric bootstrap, Bartlett, jackknife, randomization, randomization-bootstrap, and custom-randomization methods where the public capability exists.

- Exported package functions:
  - high-level exported helpers, public fitting wrappers, public matrix/data builders, simulation framework APIs, and exported diagnostics/configuration helpers.

- Package test infrastructure:
  - `package_tests/comprehensive_tests.R`
  - `package_tests/analyze_comprehensive_tests.R`
  - any new helper files used only to generate legal test cases and coverage summaries.

### Out of Scope

- Direct `EDI:::` tests of unexported implementation functions.
- Direct tests of private R6 methods or private fields.
- Tests that assert implementation details such as cache internals, workspace allocation strategy, or specific optimized code paths, except through public behavior and documented outputs.
- Exhaustive testing of illegal arguments. Illegal-argument tests are still useful, but this plan is specifically about expanding legal-argument coverage.

## Proposed Architecture

### 1. Public API Inventory

Create a generated inventory of public API surfaces:

- `api_kind`: `r6_class`, `r6_public_method`, `exported_function`, `package_runner`.
- `class_name`: class for R6 methods.
- `function_or_method`: constructor/function/method name.
- `source`: `NAMESPACE`, `R6 public_methods`, or Rd usage.
- `documented`: whether the item has an Rd entry or inherited Rd documentation.
- `response_family`: continuous, incidence, proportion, count, survival, ordinal, or all.
- `capabilities`: estimate, asymp, bootstrap, bayesian_bootstrap, parametric_bootstrap, bartlett, jackknife, rand, rand_bootstrap, design_assignment, simulation, fitting.

Implementation approach:

- Add a helper script such as `package_tests/public_api_inventory.R`.
- Use `getNamespaceExports("EDI")` for exported functions/classes.
- Use each exported R6 generator's `public_methods` and inherited public methods for method inventory.
- Use `tools::Rd_db("EDI")` or parsed `EDI/man/*.Rd` only to enrich documentation coverage, not as the sole source of truth.

Deliverable:

- `package_tests/public_api_inventory.csv`.
- A testthat check that every concrete exported R6 class has inventory rows for its constructor and public methods.

### 2. Legal Argument Registry

Create a registry describing legal arguments for public APIs. Store it as R code, not CSV, because many values are functions over fixtures.

Suggested file:

- `package_tests/legal_argument_registry.R`

Each entry should describe:

- `target`: class/function/method.
- `arg`: argument name.
- `values`: named legal values or fixture builders.
- `constraints`: predicates that determine whether a value is legal in context.
- `tier`: `smoke`, `ci`, `nightly`, `release`.
- `risk`: `low`, `medium`, `high`.
- `notes`: why the value matters.

Example registry concept:

```r
list(
  compute_bootstrap_confidence_interval = list(
    alpha = c(default = 0.05, narrow = 0.01, wide = 0.10),
    B = list(smoke = 5L, ci = 25L, nightly = 101L),
    type = c("percentile", "basic", "bca", "studentized"),
    na.rm = c(TRUE, FALSE),
    show_progress = FALSE,
    min_number_usable_samples = c(3L, 5L, 10L)
  )
)
```

This registry should distinguish:

- Constructor arguments.
- Method arguments.
- Fixture arguments used to create valid designs and responses.
- Runtime controls such as `B`, `r`, `show_progress`, root iterations, and tolerances.

### 3. Fixture Registry

Create public-API fixtures that build valid data/design/inference objects.

Suggested file:

- `package_tests/public_api_fixtures.R`

Fixture dimensions:

- `response_type`: continuous, incidence, proportion, count, survival, ordinal.
- `n`: tiny, small, medium.
- `p`: intercept-only, one covariate, multiple covariates, factor covariates, mixed numeric/factor covariates.
- `formula`: `NULL`, `~ 1`, `~ .`, reduced covariate formula, formula with factor terms.
- `design_family`: sequential, fixed, stratified, clustered, blocked, matched, optimal.
- `assignment_balance`: balanced, mildly imbalanced where legal.
- `survival_censoring`: none, light, moderate.
- `edge data`: finite but numerically challenging legal cases, such as rare incidence, overdispersed counts, boundary-near proportions not exactly 0/1 unless the public method supports it.

Fixtures must build objects only through public constructors and public methods:

- instantiate exported design class with legal constructor arguments.
- add subjects through public design methods.
- assign treatment through public design methods.
- add responses through public design methods.
- instantiate exported inference class through its public constructor.

### 4. Combination Generator

Add a bounded generator that turns the registry into test cases.

Suggested file:

- `package_tests/generate_legal_argument_cases.R`

The generator should support:

- pairwise coverage for most methods.
- targeted 3-way coverage for high-risk interactions.
- exhaustive coverage for tiny finite spaces.
- seeded randomized sampling for very large spaces.
- deterministic case IDs derived from target and argument values.
- skip reasons for unsupported legal combinations.

Recommended first implementation:

- Use a simple in-repo pairwise generator rather than adding a dependency.
- For each target, enumerate candidate values, then greedily choose cases until all legal pairs are covered.
- Apply constraints before pair accounting.
- Always include the documented default case.
- Always include current comprehensive harness cases so the new coverage report can compare old versus new.

### 5. Public Method Support Detection

Avoid hard-coded assumptions where the public API can answer support questions.

Use public methods when available:

- `get_supported_testing_types()`.
- documented public capability methods.
- exported capability table helpers if they are intended as public support.

When support can only be determined by attempting construction or method call, record that as:

- `status = "unsupported"`
- `skip_reason = <public error message>`

Do not inspect private capability methods in the combination tests. If the package needs a reliable way to answer support questions, add or document a public capability API separately before depending on it in the suite.

### 6. Public Output Invariants

For legal public calls, assert generic public invariants:

- no unexpected error.
- result has documented class/type.
- numeric outputs are finite unless documented missing/non-estimable behavior applies.
- p-values are in `[0, 1]`.
- confidence intervals are length 2 or documented structured output.
- confidence interval bounds are ordered.
- estimates are scalar or documented vector/matrix shape.
- bootstrap/randomization distributions have expected length after legal `B` or `r`, allowing documented `na.rm` behavior.
- repeated calls with the same seed and same public object state are deterministic where the API promises or implies seed determinism.
- calls with `show_progress = FALSE` do not write progress output.

Class-specific invariants should be added sparingly and only when they are public contract checks, not implementation checks.

### 7. Coverage Reporting

Extend the comprehensive result schema or add a sibling CSV:

- `api_kind`
- `class_name`
- `method`
- `case_id`
- `argument_signature`
- `fixture_id`
- `response_type`
- `design_type`
- `inference_class`
- `capability`
- `tier`
- `status`: ok, error, unsupported, skipped_slow, skipped_invalid_registry
- `duration_time_sec`
- `output_type`
- `output_shape`
- `invariant_status`
- `error_message`

Add an analyzer report showing:

- public API surfaces discovered.
- public API surfaces covered.
- uncovered exported classes/functions/methods.
- argument values covered per public method.
- pairwise coverage percentage per method.
- legal combinations skipped and why.
- slowest public calls by target and argument signature.
- newly failing combinations since the last baseline.

The coverage report should make gaps visible. A comprehensive suite is not credible unless it can say what it did not cover.

## Runtime Tiers

### Smoke Tier

Purpose:

- Fast developer feedback.

Target runtime:

- minutes, not hours.

Coverage:

- documented default case for each major public constructor/method family.
- one legal non-default argument case per high-risk method.
- tiny `B`/`r` values.
- no expensive confidence interval root searches unless already cheap.

### CI Tier

Purpose:

- Pull request regression detection.

Coverage:

- pairwise legal argument coverage for high-priority public methods.
- all response families.
- representative design families.
- representative inference families.
- strict output invariants.

### Nightly Tier

Purpose:

- broad legal argument coverage.

Coverage:

- pairwise coverage for most public methods.
- targeted 3-way coverage for resampling methods.
- more datasets and formulas.
- larger but still bounded `B`/`r`.
- slow-call accounting.

### Release Tier

Purpose:

- pre-release confidence and auditability.

Coverage:

- all concrete exported R6 classes.
- all documented public methods.
- maximal legal argument registry coverage.
- selected exhaustive sweeps for small methods.
- long-running inference paths enabled with explicit time budgets.

## Priority Order

### Phase 1: Inventory and Reporting

1. Generate public API inventory from `NAMESPACE`, R6 public methods, and Rd files.
2. Add an analyzer that reports public API coverage from existing comprehensive results.
3. Identify currently uncovered exported classes, functions, and public methods.
4. Mark internal-only tests that should be demoted, removed, or rewritten through public APIs.

Outcome:

- We know the real coverage gap before adding large amounts of work.

### Phase 2: Public Fixture Layer

1. Create small deterministic fixtures for each response family.
2. Create legal design fixtures through public constructors and public design methods.
3. Create legal inference fixtures through public constructors.
4. Add fixture validation tests.

Outcome:

- Every later argument-combination test has stable, reusable public objects.

### Phase 3: Core Argument Registry

1. Add registry entries for the most important public inference methods:
   - `compute_estimate`
   - asymptotic p-values and confidence intervals
   - bootstrap p-values and confidence intervals
   - Bayesian bootstrap p-values and confidence intervals
   - randomization p-values and confidence intervals
2. Include constructor arguments for the main design and inference classes.
3. Implement pairwise case generation and case IDs.
4. Record coverage rows without failing the suite on expected unsupported combinations.

Outcome:

- The suite begins testing legal argument interactions systematically.

### Phase 4: Full Public Inference Coverage

1. Add registry entries for parametric bootstrap, Bartlett, jackknife, m-out-of-n bootstrap, subsampling, and randomization-bootstrap methods.
2. Add targeted 3-way coverage for high-risk combinations:
   - interval `type` x `alpha` x resampling size.
   - randomization `delta` x `transform_responses` x response family.
   - bootstrap `type` x `na.rm` x minimum usable samples.
   - formula x design type x inference class.
3. Add method-specific invariants only where public documentation supports them.

Outcome:

- The expensive and interaction-heavy parts of the public API are covered without exhaustive explosion.

### Phase 5: Exported Function Coverage

1. Inventory exported non-R6 functions.
2. Classify each exported function:
   - high-level user function.
   - public fitting wrapper.
   - public data/model helper.
   - exported low-level computational routine.
3. Add legal argument registry entries for documented high-level and fitting functions first.
4. For exported low-level routines, test them only as public exported routines, with documented legal inputs and public output invariants.

Outcome:

- Exported functions receive the same argument-combination discipline as R6 methods.

### Phase 6: Quality Gates

Add gates that can be run in CI or nightly:

- no new exported public API without inventory entry.
- no public method with zero legal test cases unless explicitly exempted.
- no registry argument value that is never selected by the generator.
- no test case with unknown support status.
- no persistent new `error` rows in legal combinations.
- no undocumented direct private/internal tests in the comprehensive suite.

## Initial High-Value Argument Spaces

Start with these because they are broad, user-facing, and likely to expose interactions:

- Confidence interval methods:
  - `alpha`: `0.01`, `0.05`, `0.10`.
  - `type`: documented interval types per method.
  - resampling size: tiny CI value, normal CI value, larger nightly value.
  - `na.rm`: `TRUE`, `FALSE` where public method supports it.
  - `show_progress`: `FALSE` in automated tests, one smoke check for `TRUE` only if output is controlled.

- P-value methods:
  - `delta`: null value, small nonzero legal effect.
  - `type`: documented p-value types.
  - transform arguments: legal response transformations by response family.
  - resampling size: tier-specific.

- Constructors:
  - `response_type`: all legal response families.
  - `design_formula` / `model_formula`: `NULL`, `~ 1`, `~ .`, reduced formula, factor-containing formula.
  - stratification arguments: one stratum, two strata, factor strata.
  - cluster arguments: balanced cluster IDs and blocked cluster IDs.

- Simulation framework:
  - response family.
  - design generator.
  - inference generator.
  - treatment effect.
  - number of replications.
  - seed.
  - parallel/sequential execution if the public API exposes it.

## Handling Legal But Unsupported Combinations

Some combinations can be legal at the argument level but unsupported by a particular class/design/response context. Treat these as first-class outcomes:

- If support can be determined through public support APIs, skip before calling and record `unsupported`.
- If support is discovered through a documented public error, record `unsupported` with the message.
- If the call errors for a reason that is not a documented unsupported case, record `error` and fail the appropriate tier.

This distinction prevents the suite from pretending unsupported combinations are bugs while still preserving an audit trail.

## Helper Functions Allowed by This Plan

Helpers are allowed when they support public API testing. Examples:

- public API inventory builders.
- legal argument registries.
- pairwise case generators.
- deterministic fixture builders.
- public output invariant checkers.
- result recorders and coverage analyzers.
- skip classifiers based on public support methods or documented public error messages.

Helpers should not:

- inspect private R6 fields to decide what to test.
- call private R6 methods.
- call unexported package functions directly.
- assert on internal caches, workspaces, or implementation path choices.

## Concrete Deliverables

1. `package_tests/public_api_inventory.R`
2. `package_tests/legal_argument_registry.R`
3. `package_tests/public_api_fixtures.R`
4. `package_tests/generate_legal_argument_cases.R`
5. `package_tests/run_public_api_argument_matrix.R`
6. `package_tests/analyze_public_api_argument_coverage.R`
7. `package_tests/public_api_argument_coverage.csv`
8. `package_tests/public_api_argument_failures.csv`
9. CI/nightly scripts that select smoke, CI, nightly, or release tiers.

## Definition of Done

The test suite can be considered professionally comprehensive when:

- Every exported class/function is inventoried.
- Every concrete exported R6 class has at least one public constructor test.
- Every documented public method has either legal coverage or an explicit exemption.
- High-priority public inference methods have pairwise legal-argument coverage.
- High-risk methods have targeted higher-order interaction coverage.
- Coverage reports show argument values and argument pairs covered per public method.
- Unsupported combinations are recorded separately from failures.
- Runtime tiers are stable and documented.
- Direct internal implementation tests are not required to claim public API coverage.

