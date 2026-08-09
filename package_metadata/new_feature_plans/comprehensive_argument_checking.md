# Spec: Comprehensive Legal Argument Combination Checking

## Objective

Implement systematic tests for legal combinations of arguments passed into exported public API functions and public R6 constructors/methods. The package already has many strong one-argument checks through `checkmate` and custom assertions, but those mostly prove that individual argument values are legal in isolation. This spec adds multidimensional checking: individually legal values must also work together in real public API calls.

Example target behavior:

- `alpha = 0.01` is legal.
- `type = "studentized"` is legal.
- `B = 25L` is legal.
- `na.rm = FALSE` is legal.
- The joint call `compute_bootstrap_confidence_interval(alpha = 0.01, type = "studentized", B = 25L, na.rm = FALSE, show_progress = FALSE)` is legal for a supported fixture and returns a valid public result or a documented unsupported/non-estimable status.

## Non-Goals

- Do not exhaustively test illegal arguments.
- Do not blindly run Cartesian products over all legal values.
- Do not treat direct private/internal tests as public API coverage.
- Do not assert incidental private representation, cache shape, or branch choices unless the implementation surface has an explicit internal contract.

## Terminology

- `local contract`: A one-argument legality rule, usually extracted from `checkmate`, a custom assertion, a formal default, or Rd documentation.
- `combination constraint`: A predicate over two or more arguments and fixture context.
- `fixture`: A deterministic object/data setup that makes a public API call meaningful.
- `case`: One generated public API call with a concrete fixture and argument set.
- `coverage scope`: `public_contract` or `internal_safety_net`.
- `tier`: `smoke`, `ci`, `nightly`, or `release`.

## Architecture

The implementation has five layers:

1. Public API inventory.
2. Local argument contract extraction.
3. Registry enrichment with values, constraints, tiers, and invariants.
4. Combination generation and execution.
5. Coverage/failure analysis.

The new runner should complement `package_tests/comprehensive_tests.R`. The comprehensive harness remains the broad end-to-end statistical workflow runner. The new runner focuses on public API argument-combination coverage and produces joinable result files.

## Artifacts

Implement these files in `package_tests/`:

- `public_api_inventory.R`
- `extract_checkmate_argument_contracts.R`
- `checkmate_argument_contracts.csv`
- `public_argument_contract_registry.R`
- `public_argument_combination_constraints.R`
- `public_argument_combination_fixtures.R`
- `generate_public_argument_combinations.R`
- `run_public_argument_combinations.R`
- `analyze_public_argument_combinations.R`
- `public_argument_combination_cases.csv`
- `public_argument_combination_results.csv`
- `public_argument_combination_coverage.csv`
- `public_argument_combination_failures.csv`

## Phase 0: Baseline Audit

Goal: establish the current gap before implementing a new runner.

TODO:

- [x] Generate a list of exported symbols with `getNamespaceExports("EDI")`.
- [x] Classify each exported symbol as `r6_class`, `function`, `data`, or `other`.
- [x] For exported R6 classes, list inherited and class-local public methods.
- [x] Count configurable arguments per exported function, constructor, and public method.
- [x] Identify APIs with more than one configurable argument.
- [x] Identify APIs with `checkmate` or custom assertion calls.
- [x] Identify APIs currently reached by `package_tests/comprehensive_tests.R`.
- [x] Identify APIs currently only covered by one-off `testthat` tests.
- [x] Produce a baseline table of APIs that have only unidimensional checks.

Outputs:

- `package_tests/public_api_inventory.csv`
- `package_tests/public_argument_baseline_gap_report.csv`

Acceptance criteria:

- [x] Every exported symbol appears in the inventory.
- [x] Every concrete exported R6 class has constructor and public-method rows.
- [x] The gap report can answer: "which public APIs have multiple configurable arguments and no legal-combination coverage?"

Completed outputs:

- `package_tests/public_api_inventory.R`
- `package_tests/public_api_inventory.csv`
- `package_tests/public_argument_baseline_gap_report.csv`

Completion notes:

- Inventory rows: 9216.
- Exported symbol rows: 249.
- Exported R6 classes: 128.
- R6 public method rows: 8967.
- Exported R6 classes missing `initialize` rows: 0.
- Baseline gap rows: 3514.
- Baseline gap rows with detected assertion/check calls: 1891.

## Phase 1: Checkmate Contract Extraction

Goal: extract local argument contracts from source assertions.

Implementation notes:

- Use `utils::getParseData(parse(file, keep.source = TRUE))`.
- Use regex only as a fallback for constructs the parser pass cannot recover.
- Capture source file and line for every extracted contract.
- Detect calls inside and outside `should_run_asserts()` blocks.
- Support checkmate assertions and package-specific wrappers.

Initial assertions to support:

- `assertChoice`
- `assertSubset`
- `assertNumeric`
- `assertIntegerish`
- `assertCount`
- `assertFlag`
- `assertLogical`
- `assertString`
- `assertFormula`
- `assertList`
- `assertClass`
- `assertDataFrame`
- `assertMatrix`
- `assertResponseType`
- `assertNoCensoring`

TODO:

- [x] Implement parser that maps R source files to enclosing exported function/R6 method.
- [x] Extract formal arguments and default expressions.
- [x] Extract assertion call, asserted argument, and named assertion parameters.
- [x] Extract literal `choices` vectors from `assertChoice()` and `assertSubset()`.
- [x] Extract `lower`, `upper`, `len`, `min.len`, `max.len`, `positive`, `null.ok`, and `any.missing`.
- [x] Mark contracts as `public_contract` or `internal_safety_net`.
- [x] Preserve file and line for each contract.
- [x] Add a smoke test for the extractor using representative files with `assertChoice`, `assertNumeric`, `assertCount`, and `assertFormula`.

Contract schema:

- `api_name`
- `api_kind`
- `class_name`
- `method_name`
- `arg`
- `default_expr`
- `assertion`
- `choices`
- `lower`
- `upper`
- `len`
- `min_len`
- `max_len`
- `positive`
- `null_ok`
- `any_missing`
- `source_file`
- `source_line`
- `coverage_scope`

Outputs:

- `package_tests/checkmate_argument_contracts.csv`

Acceptance criteria:

- [x] Contracts are extracted for high-volume assertion families used by inference/design APIs.
- [x] The extractor records enough source location to manually inspect every contract.
- [x] Re-running extraction is deterministic.

Completed outputs:

- `package_tests/extract_checkmate_argument_contracts.R`
- `package_tests/checkmate_argument_contracts.csv`
- `package_tests/smoke_test_extract_checkmate_argument_contracts.R`

Completion notes:

- Parsed R source files with `keep.source = TRUE`: 179 of 179.
- Public function/method inventory rows scanned: 9088.
- Extracted contract rows: 6517.
- Extracted assertion families: `assertChoice`, `assertClass`, `assertCount`, `assertFlag`, `assertFormula`, `assertIntegerish`, `assertList`, `assertLogical`, `assertNoCensoring`, `assertNumeric`, `assertResponseType`, `assertSubset`.
- `assertChoice` rows: 628.
- `assertNumeric` rows: 2278.
- `assertCount` rows: 2257.
- `assertFormula` rows: 30.
- Rows inside `should_run_asserts()` blocks: 6506.
- Rows with source lines: 6517.
- Coverage scope rows marked `public_contract`: 6517.
- Deterministic CSV hash after two runs: `e0df4ef5885d7772a21e9958ddafeabdee79af9b85a911f34ad40e5e32fd0e60`.

## Phase 2: Registry Enrichment

Goal: turn local contracts into runnable domains and cross-argument constraints.

Checkmate-derived contracts seed the registry, but the registry is the reviewed source for test generation.

TODO:

- [x] Create `public_argument_contract_registry.R`.
- [x] Add domain constructors for common assertion types.
- [x] Add representative values for numeric ranges:
  - default.
  - normal interior value.
  - legal value near lower bound.
  - legal value near upper bound where useful.
- [x] Add tier-specific values for runtime arguments such as `B`, `r`, root iterations, and progress flags.
- [x] Add full choice sets for small `assertChoice` domains.
- [x] Add CI/nightly sampling rules for large choice domains.
- [x] Add `NULL`, `~ 1`, `~ .`, and reduced fixture-specific formulas for legal formula arguments.
- [x] Add registry metadata for optional dependency requirements.
- [x] Add registry metadata for known slow paths.
- [x] Add registry metadata for documented unsupported contexts.

Registry entry shape:

```r
list(
  target = "compute_bootstrap_confidence_interval",
  api_kind = "r6_public_method",
  domains = list(
    alpha = c(default = 0.05, narrow = 0.01, wide = 0.10),
    B = list(smoke = 5L, ci = 25L, nightly = 101L),
    type = c("percentile", "basic", "bca", "studentized"),
    na.rm = c(TRUE, FALSE),
    show_progress = FALSE,
    min_number_usable_samples = c(3L, 5L, 10L)
  ),
  constraints = list(
    "min_number_usable_samples <= B",
    "type != 'bca' || fixture$supports_jackknife"
  ),
  invariants = c("ci_length_2", "ci_ordered", "ci_numeric_or_documented_missing")
)
```

Acceptance criteria:

- [x] Every high-priority public inference method with extracted argument contracts has a registry entry.
- [x] Every registry entry includes domains, constraints, invariants, and tier metadata.
- [x] Registry values can be traced back to checkmate, defaults/Rd, or explicit hand-written rationale.

Completed outputs:

- `package_tests/public_argument_contract_registry.R`
- `package_tests/public_argument_contract_registry.rds`
- `package_tests/public_argument_contract_registry.csv`
- `package_tests/smoke_test_public_argument_contract_registry.R`

Completion notes:

- Registry entries: 2443.
- Flattened registry value rows: 23394.
- Distinct target/argument domains: 6468.
- Tiers represented: `smoke`, `ci`, `nightly`, `release`.
- Value sources represented: `checkmate_bound`, `checkmate_choices`, `checkmate_count`, `checkmate_flag`, `fixture`, `formal_default`, `formula_rule`, `hand_rule`, `runtime_tier`, `unsupported_context`.
- High-priority methods represented: `compute_asymp_confidence_interval`, `compute_asymp_two_sided_pval`, `compute_bayesian_bootstrap_confidence_interval`, `compute_bayesian_bootstrap_two_sided_pval`, `compute_bootstrap_confidence_interval`, `compute_bootstrap_two_sided_pval`, `compute_rand_confidence_interval`, `compute_rand_two_sided_pval`.
- `compute_estimate` has no Phase 1 argument-contract rows because it has no extracted configurable checked arguments.
- Registry CSV hash: `896bdad2f185a312c1766831d169725721a7019890d173e2c6e82256e036b968`.
- Registry RDS hash: `72ee18ca87b9b217cc4b44a64c3a1981e48f77940b947668b919c1d88115ab1f`.

## Phase 3: Cross-Argument Constraint Library

Goal: encode legal interactions explicitly instead of silently skipping them.

Implement `public_argument_combination_constraints.R`.

TODO:

- [x] `min_number_usable_samples <= B`.
- [x] `type = "bca"` requires jackknife support and enough exchangeable units.
- [x] `type = "studentized"` requires finite standard errors or documented non-estimable behavior.
- [x] `use_rcpp = FALSE` requires the documented optional backend package.
- [x] `show_progress = TRUE` is not used in automated CI unless output is intentionally captured.
- [x] `transform_responses` values are compatible with response family and data support.
- [x] formulas reference only columns present in the fixture.
- [x] `strata_cols` reference existing covariates and are factor/discretized when the public design requires that.
- [x] `cluster_col` references valid cluster IDs.
- [x] `strata_cols` and `cluster_col` do not conflict unless the design class supports the combination.
- [x] fixed/blocking designs receive compatible `n`, block sizes, and strata levels.
- [x] survival methods with censoring use methods that publicly support censoring.
- [x] optional dependency paths are marked `skipped_dependency` when unavailable.
- [x] known slow legal combinations are marked `skipped_slow` below their allowed tier.

Constraint result schema:

- `is_valid`: TRUE/FALSE.
- `status`: `valid`, `unsupported`, `skipped_slow`, `skipped_dependency`, `invalid_registry`.
- `reason`: human-readable reason.
- `source`: registry key or constraint function name.

Acceptance criteria:

- [x] No generated case is discarded without a recorded constraint result.
- [x] Impossible combinations do not count against pairwise coverage.
- [x] Unsupported legal contexts are distinguished from bugs.

Completed outputs:

- `package_tests/public_argument_combination_constraints.R`
- `package_tests/smoke_test_public_argument_combination_constraints.R`

Completion notes:

- The constraint library evaluates every candidate against every Phase 3 constraint and returns auditable rows with `constraint_name`, `is_valid`, `status`, `reason`, and `source`.
- `constraint_context_result()` summarizes detailed rows into the case-level status that later phases can write into generated case and result tables.
- The smoke test covers `valid`, `unsupported`, `skipped_slow`, `skipped_dependency`, and `invalid_registry`.

## Phase 4: Public Fixture Layer

Goal: create deterministic public API fixtures that make legal combinations executable.

Implement `public_argument_combination_fixtures.R`.

Fixture families:

- continuous.
- incidence.
- proportion.
- count.
- survival without censoring.
- survival with light censoring.
- ordinal.

Design fixture dimensions:

- sequential design.
- fixed design.
- stratified design.
- clustered design.
- blocked-cluster design.
- matched or pair-based design where legal.
- optimal/search design where runtime permits.

Data fixture dimensions:

- `n`: tiny and small.
- `p`: intercept-only, one numeric covariate, mixed numeric/factor covariates.
- formula support: `~ 1`, `~ .`, reduced formula, factor-containing formula.
- legal edge cases: rare incidence, overdispersed counts, boundary-near proportions, and censored survival where supported.

TODO:

- [x] Build fixtures only through exported constructors and public methods.
- [x] Assign treatments through public design methods.
- [x] Add responses through public design methods.
- [x] Expose fixture metadata needed by constraints.
- [x] Add validation functions for fixture consistency.
- [x] Keep fixture sizes small enough for CI.

Fixture metadata:

- `fixture_id`
- `response_type`
- `design_type`
- `n`
- `p`
- `has_strata`
- `has_cluster`
- `has_matching`
- `has_censoring`
- `supports_jackknife`
- `supports_randomization`
- `available_columns`
- `tier`

Acceptance criteria:

- [x] Every response family has at least one smoke fixture.
- [x] High-priority inference methods have at least one compatible fixture.
- [x] Fixture creation itself is tested and deterministic.

Completed outputs:

- `package_tests/public_argument_combination_fixtures.R`
- `package_tests/smoke_test_public_argument_combination_fixtures.R`

Completion notes:

- Smoke fixtures cover continuous, incidence, proportion, count, survival without censoring, survival with light censoring, and ordinal responses.
- Smoke design fixtures cover fixed, sequential, stratified, clustered, blocked-cluster, matched, and low-iteration greedy/search designs.
- Fixed fixtures use exported R6 constructors plus public `add_all_subjects_to_experiment()`, `assign_w_to_all_subjects()`, and `add_all_subject_responses()` methods.
- Sequential fixtures use exported R6 constructors plus public `add_one_subject_to_experiment_and_assign()` and `add_one_subject_response()` methods.
- All fixture metadata is validated against public getters such as `get_response_type()`, `get_n()`, `get_t()`, `get_w()`, and `any_censoring()`.

## Phase 5: Combination Generator

Goal: generate bounded legal combinations with measurable coverage.

Implement `generate_public_argument_combinations.R`.

Generation algorithm:

1. Build candidate domains from registry.
2. Expand fixture contexts for each target.
3. Always include default call.
4. Include one non-default legal value per configurable argument.
5. Generate pairwise cases over legal values.
6. Generate targeted 3-way cases for high-risk APIs.
7. Use exhaustive generation only for tiny finite domains.
8. Apply constraints before counting coverage.
9. Assign deterministic `case_id`.

TODO:

- [x] Implement default-case generation.
- [x] Implement one-non-default-per-argument generation.
- [x] Implement greedy pairwise generation.
- [x] Implement targeted 3-way generation hooks.
- [x] Implement deterministic case IDs from target, fixture, and argument signature.
- [x] Emit all cases and rejected candidate summaries.
- [x] Add unit tests for pairwise coverage accounting.

Case schema:

- `case_id`
- `api_kind`
- `class_name`
- `function_name`
- `method_name`
- `fixture_id`
- `tier`
- `argument_signature`
- `argument_coverage_kind`
- `source_contracts`
- `constraint_status`
- `constraint_reason`

Outputs:

- `package_tests/public_argument_combination_cases.csv`

Acceptance criteria:

- [x] Default cases are always present.
- [x] Every generated case has a deterministic ID.
- [x] Pairwise coverage can be recomputed from the generated case table.
- [x] Constraint-rejected candidates are auditable.

Completed outputs:

- `package_tests/generate_public_argument_combinations.R`
- `package_tests/public_argument_combination_cases.csv`
- `package_tests/public_argument_combination_rejected_candidates.csv`
- `package_tests/public_argument_combination_coverage.csv`
- `package_tests/smoke_test_generate_public_argument_combinations.R`

Completion notes:

- The generator reads the Phase 2 registry RDS, expands compatible Phase 4 fixtures, and applies Phase 3 constraints before coverage accounting.
- It emits default, one-non-default, pairwise, targeted 3-way, and tiny-exhaustive candidate rows.
- Constrained-out candidates remain auditable through `constraint_status`, `constraint_reason`, `constraint_sources`, and the rejected-candidate CSV.

## Phase 6: Runner

Goal: execute generated cases against exported public APIs and record structured results.

Implement `run_public_argument_combinations.R`.

Runner behavior:

- Load installed/package-under-test `EDI`.
- Load fixture builders, registry, constraints, and generated cases.
- Select cases by tier and optional filters.
- Execute public constructors/functions/methods only.
- Use timeouts for each case.
- Capture warnings and messages.
- Record output type/shape and invariant results.
- Distinguish `ok`, `error`, `unsupported`, `skipped_slow`, `skipped_dependency`, and `invalid_registry`.

Suggested CLI:

```sh
Rscript package_tests/run_public_argument_combinations.R smoke
Rscript package_tests/run_public_argument_combinations.R ci
Rscript package_tests/run_public_argument_combinations.R nightly
Rscript package_tests/run_public_argument_combinations.R release
Rscript package_tests/run_public_argument_combinations.R ci InferencePropBetaRegr compute_bootstrap_confidence_interval
```

TODO:

- [ ] Implement tier argument parsing.
- [ ] Implement optional target filters.
- [ ] Implement per-case timeout.
- [ ] Implement warning/message capture.
- [ ] Implement resumable output keyed by `case_id`.
- [ ] Implement invariant checks.
- [ ] Implement final nonzero exit when CI-tier unexpected errors exist.

Result schema:

- `case_id`
- `api_kind`
- `class_name`
- `function_name`
- `method_name`
- `fixture_id`
- `argument_signature`
- `tier`
- `status`
- `duration_time_sec`
- `warning_message`
- `message_text`
- `output_type`
- `output_shape`
- `invariant_status`
- `error_message`

Outputs:

- `package_tests/public_argument_combination_results.csv`
- `package_tests/public_argument_combination_failures.csv`

Acceptance criteria:

- Smoke tier runs end-to-end.
- Results are resumable by `case_id`.
- Unexpected errors are separated from unsupported or dependency-skipped cases.

## Phase 7: Invariant Library

Goal: validate public outputs consistently.

Implement reusable invariant functions.

TODO:

- [ ] `no_unexpected_error`
- [ ] `numeric_finite_or_documented_missing`
- [ ] `pval_in_0_1`
- [ ] `ci_length_2`
- [ ] `ci_ordered`
- [ ] `estimate_shape_valid`
- [ ] `distribution_length_valid`
- [ ] `documented_nonestimable_ok`
- [ ] `seeded_call_deterministic`
- [ ] `output_names_stable_when_documented`

Acceptance criteria:

- Invariants return structured status and reason.
- Method-specific invariants are opt-in through registry.
- Missing/non-estimable outputs are accepted only when documented or explicitly classified.

## Phase 8: Analyzer and Coverage Report

Goal: make coverage gaps and failures visible.

Implement `analyze_public_argument_combinations.R`.

Coverage metrics:

- exported public APIs with at least one legal-combination case.
- exported public APIs with multiple configurable arguments but no combination cases.
- public APIs with only unidimensional checks.
- public APIs with checkmate contracts but no multidimensional cases.
- argument values covered by API.
- legal pairs covered by API.
- targeted 3-way interactions covered by API.
- skipped legal combinations and reasons.
- failing legal combinations.
- slowest legal combinations.
- drift between checkmate contracts and registry.

TODO:

- [ ] Join inventory, contracts, registry, generated cases, and results.
- [ ] Compute argument-value coverage.
- [ ] Compute pairwise coverage from valid candidate space.
- [ ] Compute high-risk 3-way coverage.
- [ ] Emit uncovered API report.
- [ ] Emit drift report.
- [ ] Emit slow-case report.
- [ ] Emit CI failure summary.

Outputs:

- `package_tests/public_argument_combination_coverage.csv`
- `package_tests/public_argument_combination_failures.csv`
- optional `package_tests/public_argument_combination_report.html`

Acceptance criteria:

- Analyzer identifies APIs that still only have unidimensional checks.
- Analyzer identifies registry drift from checkmate contracts.
- CI can fail on new unexpected errors without failing on documented unsupported contexts.

## Phase 9: Integration With Existing Harness

Goal: make the new runner complement, not duplicate, `comprehensive_tests.R`.

TODO:

- [ ] Reuse response/design/inference labels compatible with comprehensive results.
- [ ] Reuse or wrap fixture-building logic where practical.
- [ ] Keep new argument-combination result files separate from statistical comprehensive results.
- [ ] Add join keys for `response_type`, `design_type`, `inference_class`, `function_run`, and `dataset_name` where available.
- [ ] Add documentation explaining when to use each runner.
- [ ] Add CI/nightly entry points without changing existing long-running comprehensive commands.

Acceptance criteria:

- Existing comprehensive tests continue to run unchanged.
- The new runner can be run independently.
- Failures from the new runner can be promoted into focused `testthat` regression tests.

## Phase 10: CI and Quality Gates

Goal: enforce the new coverage discipline incrementally.

TODO:

- [ ] No new exported public API without inventory entry.
- [ ] No public API with multiple configurable arguments and zero combination coverage unless explicitly exempted.
- [ ] No checkmate-derived public argument contract without registry coverage or exemption.
- [ ] No registry value that is never selected by any tier.
- [ ] No generated case with unknown constraint status.
- [ ] No CI-tier unexpected `error` rows.

Rollout:

1. Start gates as reporting-only.
2. Enable hard failure for extractor determinism and schema validity.
3. Enable hard failure for smoke-tier unexpected errors.
4. Enable hard failure for new uncovered high-priority APIs.
5. Enable hard failure for CI-tier unexpected errors.

Acceptance criteria:

- CI has a small, stable smoke/CI tier.
- Nightly/release tiers can be run manually or scheduled.
- Coverage gates improve over time without blocking unrelated work prematurely.

## High-Priority Targets

Start with APIs that have many arguments or high interaction risk:

- Public inference constructors with `model_formula`, `verbose`, `use_rcpp`, `optimization_alg`, `smart_cold_start_default`, and family-specific options.
- Public design constructors with `response_type`, `n`, `design_formula`, `strata_cols`, `cluster_col`, block parameters, matching parameters, optimization/search controls, and seeds.
- Public inference methods:
  - `compute_asymp_two_sided_pval`
  - `compute_asymp_confidence_interval`
  - `compute_wald_*`
  - `compute_score_*`
  - `compute_lik_ratio_*`
  - `compute_gradient_*`
  - `compute_bootstrap_*`
  - `compute_bayesian_bootstrap_*`
  - `compute_lik_ratio_bootstrap_*`
  - `compute_param_bootstrap_*`
  - `compute_lik_ratio_bartlett_*`
  - `compute_jackknife_*`
  - `compute_rand_*`
  - `compute_rand_bootstrap_*`
- Exported public fitting and helper functions with multiple optional arguments.
- `SimulationFramework` and related report APIs.

## Definition of Done

- [ ] Checkmate-derived argument contracts are extracted for exported public APIs.
- [ ] Every exported public API with more than one configurable argument has legal-combination coverage or an explicit exemption.
- [ ] High-priority APIs have pairwise legal-argument coverage.
- [ ] High-risk APIs have targeted 3-way legal-argument coverage.
- [ ] Cross-argument constraints are explicit, versioned, and reported.
- [ ] Failures distinguish unsupported legal contexts from bugs.
- [ ] Coverage reports identify APIs that still only have unidimensional argument checks.
- [ ] Smoke/CI tiers run with stable runtime and deterministic output.
- [ ] Nightly/release tiers provide broader coverage without changing the CI contract.
