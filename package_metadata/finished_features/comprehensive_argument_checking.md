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
- `public_argument_combination_integration.R`
- `check_public_argument_combination_quality_gates.R`
- `run_public_argument_combination_ci.R`
- `run_public_argument_combination_nightly.R`
- `public_argument_combination_runner.md`
- `public_argument_combination_quality_gates.md`
- `public_argument_combination_cases.csv`
- `public_argument_combination_results.csv`
- `public_argument_combination_coverage.csv`
- `public_argument_combination_failures.csv`
- `public_argument_combination_cases_integrated.csv`
- `public_argument_combination_results_integrated.csv`
- `public_argument_combination_coverage_integrated.csv`

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

- [x] Implement tier argument parsing.
- [x] Implement optional target filters.
- [x] Implement per-case timeout.
- [x] Implement warning/message capture.
- [x] Implement resumable output keyed by `case_id`.
- [x] Implement invariant checks.
- [x] Implement final nonzero exit when CI-tier unexpected errors exist.

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

- [x] Smoke tier runs end-to-end.
- [x] Results are resumable by `case_id`.
- [x] Unexpected errors are separated from unsupported or dependency-skipped cases.

Completed outputs:

- `package_tests/run_public_argument_combinations.R`
- `package_tests/public_argument_combination_results.csv`
- `package_tests/public_argument_combination_failures.csv`
- `package_tests/smoke_test_run_public_argument_combinations.R`

Completion notes:

- The runner records constrained-out cases without executing them and preserves their `unsupported`, `skipped_slow`, `skipped_dependency`, or `invalid_registry` status.
- Valid cases are executed through exported functions or exported R6 class generators and public methods only.
- Warning and message text, output type/shape, invariant status, duration, and errors are recorded per `case_id`.
- Existing result files are resumed by skipping already-recorded `case_id`s.

## Phase 7: Invariant Library

Goal: validate public outputs consistently.

Implement reusable invariant functions.

TODO:

- [x] `no_unexpected_error`
- [x] `numeric_finite_or_documented_missing`
- [x] `pval_in_0_1`
- [x] `ci_length_2`
- [x] `ci_ordered`
- [x] `estimate_shape_valid`
- [x] `distribution_length_valid`
- [x] `documented_nonestimable_ok`
- [x] `seeded_call_deterministic`
- [x] `output_names_stable_when_documented`

Acceptance criteria:

- [x] Invariants return structured status and reason.
- [x] Method-specific invariants are opt-in through registry.
- [x] Missing/non-estimable outputs are accepted only when documented or explicitly classified.

Completed outputs:

- `package_tests/public_argument_combination_invariants.R`
- `package_tests/smoke_test_public_argument_combination_invariants.R`

Completion notes:

- Invariants return one row per requested check with `invariant`, `status`, and `reason`.
- `evaluate_public_argument_invariants()` runs only the semicolon-delimited invariant names registered on each generated case.
- `summarize_public_argument_invariants()` preserves the Phase 6 runner's compact `invariant_status` output while retaining reusable structured checks.
- Documented missing/non-estimable numeric outputs are accepted only when context explicitly marks them as documented.

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

- [x] Join inventory, contracts, registry, generated cases, and results.
- [x] Compute argument-value coverage.
- [x] Compute pairwise coverage from valid candidate space.
- [x] Compute high-risk 3-way coverage.
- [x] Emit uncovered API report.
- [x] Emit drift report.
- [x] Emit slow-case report.
- [x] Emit CI failure summary.

Outputs:

- `package_tests/public_argument_combination_coverage.csv`
- `package_tests/public_argument_combination_failures.csv`
- optional `package_tests/public_argument_combination_report.html`

Acceptance criteria:

- [x] Analyzer identifies APIs that still only have unidimensional checks.
- [x] Analyzer identifies registry drift from checkmate contracts.
- [x] CI can fail on new unexpected errors without failing on documented unsupported contexts.

Completed outputs:

- `package_tests/analyze_public_argument_combinations.R`
- `package_tests/public_argument_combination_coverage.csv`
- `package_tests/public_argument_combination_failures.csv`
- `package_tests/public_argument_combination_uncovered_apis.csv`
- `package_tests/public_argument_combination_registry_drift.csv`
- `package_tests/public_argument_combination_slowest_cases.csv`
- `package_tests/public_argument_combination_ci_failure_summary.csv`
- `package_tests/public_argument_combination_report.html`
- `package_tests/smoke_test_analyze_public_argument_combinations.R`

Completion notes:

- The analyzer joins the Phase 0 inventory, Phase 1 contracts, Phase 2 registry, Phase 5 cases, and Phase 6 results.
- Coverage rows include target-level case counts, valid candidate counts, pairwise-or-higher counts, targeted 3-way counts, executed-ok counts, skip counts, and covered argument-value/pair counts.
- Drift rows distinguish missing registry arguments, registry arguments without extracted contracts, registry targets absent from inventory, and registry targets not represented in the current generated case table.
- CI failure summary treats only unexpected execution errors and invariant failures as blocking; unsupported, dependency-skipped, and slow-skipped cases remain non-blocking.

## Phase 9: Integration With Existing Harness

Goal: make the new runner complement, not duplicate, `comprehensive_tests.R`.

TODO:

- [x] Reuse response/design/inference labels compatible with comprehensive results.
- [x] Reuse or wrap fixture-building logic where practical.
- [x] Keep new argument-combination result files separate from statistical comprehensive results.
- [x] Add join keys for `response_type`, `design_type`, `inference_class`, `function_run`, and `dataset_name` where available.
- [x] Add documentation explaining when to use each runner.
- [x] Add CI/nightly entry points without changing existing long-running comprehensive commands.

Acceptance criteria:

- [x] Existing comprehensive tests continue to run unchanged.
- [x] The new runner can be run independently.
- [x] Failures from the new runner can be promoted into focused `testthat` regression tests.

Completed outputs:

- `package_tests/public_argument_combination_integration.R`
- `package_tests/public_argument_combination_cases_integrated.csv`
- `package_tests/public_argument_combination_results_integrated.csv`
- `package_tests/public_argument_combination_coverage_integrated.csv`
- `package_tests/run_public_argument_combination_ci.R`
- `package_tests/run_public_argument_combination_nightly.R`
- `package_tests/public_argument_combination_runner.md`
- `package_tests/smoke_test_public_argument_combination_integration.R`

Completion notes:

- Integrated case rows: 10.
- Integrated result rows: 10.
- Integrated coverage rows: 9216.
- Added fixture-derived join keys: `response_type`, `design_type`, `design`, `dataset_name`, `dataset`, and `dataset_source`.
- Added public API target-derived join keys: `inference_class` and `function_run`.
- Comprehensive-compatible design labels use the existing `comprehensive_tests.R` naming style, e.g. `DesignFixedBernoulli` becomes `FixedBernoulli`.
- The new CI and nightly entry points write only `public_argument_combination_*` artifacts and do not modify `package_tests/comprehensive_tests.R` or its long-running commands.
- Documented when to use `comprehensive_tests.R` versus the argument-combination runner and how to promote argument-combination failures into focused `testthat` regression tests.

## Phase 10: CI and Quality Gates

Goal: enforce the new coverage discipline incrementally.

TODO:

- [x] No new exported public API without inventory entry.
- [x] No public API with multiple configurable arguments and zero combination coverage unless explicitly exempted.
- [x] No checkmate-derived public argument contract without registry coverage or exemption.
- [x] No registry value that is never selected by any tier.
- [x] No generated case with unknown constraint status.
- [x] No CI-tier unexpected `error` rows.

Rollout:

1. Start gates as reporting-only.
2. Enable hard failure for extractor determinism and schema validity.
3. Enable hard failure for smoke-tier unexpected errors.
4. Enable hard failure for new uncovered high-priority APIs.
5. Enable hard failure for CI-tier unexpected errors.

Acceptance criteria:

- [x] CI has a small, stable smoke/CI tier.
- [x] Nightly/release tiers can be run manually or scheduled.
- [x] Coverage gates improve over time without blocking unrelated work prematurely.

Completed outputs:

- `package_tests/check_public_argument_combination_quality_gates.R`
- `package_tests/public_argument_combination_quality_gates.csv`
- `package_tests/public_argument_combination_quality_gate_summary.csv`
- `package_tests/public_argument_combination_gate_exemptions.csv`
- `package_tests/public_argument_combination_quality_gates.md`
- `package_tests/smoke_test_public_argument_combination_quality_gates.R`
- Updated `package_tests/run_public_argument_combination_ci.R`
- Updated `package_tests/run_public_argument_combination_nightly.R`

Completion notes:

- Quality gate rows in the current smoke/CI artifact set: 23391.
- Active hard gate rows in `ci` mode: 0.
- Implemented reporting gates for uncovered multi-argument APIs, checkmate contracts without registry coverage, and registry values not selected by the current artifact set.
- Implemented hard integrity gates for missing exported inventory rows, schema validity, unknown generated constraint statuses, and smoke/CI unexpected errors or invariant failures.
- Added `strict` mode to promote future hard gates such as uncovered high-priority APIs after coverage has matured.
- Added an empty exemption CSV with exact-match keys for temporary documented exceptions.
- The new CI and nightly entry points now run the quality gates after generation, execution, analysis, and integration.
- Roxygen assessment: no roxygen documentation is required for this phase because the new argument-checking paradigm added package test scripts, generated artifacts, and package metadata only; it did not add or modify exported `EDI` package functions/classes or user-facing package APIs. User-facing runner documentation is in `package_tests/public_argument_combination_runner.md` and gate documentation is in `package_tests/public_argument_combination_quality_gates.md`.

## Phase 11: Promote Stable Combination Rules Into Runtime Assertions

Goal: rotate durable cross-argument rules discovered by the argument-combination runner into actual public runtime validation, while keeping broad coverage generation in `package_tests/`.

Rationale:

- The current runner is a mining, coverage, and regression-generation layer.
- Public users benefit when stable, general cross-argument failures are caught immediately in constructors, public methods, and exported functions.
- Assertion runtime is not a primary blocker for promotion because users can disable assertion checks through the package's existing assertion controls.
- Not every generated constraint belongs in runtime assertions; fixture-only, slow-path, optional-backend, and testing-tier constraints should remain in the test/registry layer.

TODO:

- [x] Keep the `package_tests/` runner as the broad pairwise, targeted 3-way, and tiny-exhaustive legal-combination coverage system.
- [x] Audit generated constraints and quality-gate reports for stable, general rules that should become runtime API checks.
- [x] Do not promote fixture-only rules into runtime checks.
- [x] Do not promote slow-path labels into runtime checks.
- [x] Do not promote tier-specific test-generation rules into runtime checks.
- [x] Do not promote optional-dependency availability rules unless the public API already commits to checking that dependency at call time.
- [x] Add helper assertions for bootstrap and resampling arguments, such as `assertBootstrapArgs(type, B, min_number_usable_samples, ...)`.
- [x] Add helper assertions for model formula/data context, such as `assertFormulaContext(model_formula, data, ...)`.
- [x] Add helper assertions for strata and cluster relationships, such as `assertStrataClusterArgs(strata_cols, cluster_col, data, ...)`.
- [x] Add helper assertions for durable response/design/inference compatibility rules where they are not already enforced.
- [x] Call the helper assertions from public constructors, exported functions, and public R6 methods after existing one-dimensional `checkmate` checks.
- [x] Ensure every new promoted runtime check is inside the existing `should_run_asserts()` wrapper.
- [x] Preserve existing `should_run_asserts()` behavior so promoted checks run only when package assertions are enabled.
- [x] Ensure promoted assertions produce direct, user-actionable error messages that name the conflicting arguments.
- [x] Add focused `testthat` tests for each promoted runtime rule.
- [x] Update roxygen/Rd documentation only if a promoted rule changes documented user-facing behavior or a new exported helper is introduced.

Initial candidate runtime rules:

- `min_number_usable_samples <= B` for bootstrap-like APIs.
- Incompatible `type`/method combinations where the incompatibility is general and independent of a fixture.
- Invalid `strata_cols`/`cluster_col` relationships, including missing columns, overlap rules, and factor requirements where the design requires factor strata.
- Formula/context mismatches where `model_formula` references columns unavailable in the public call context.
- Response/design/inference compatibility rules that are currently rediscovered only by fixture execution but are true for all users.

Acceptance criteria:

- [x] Stable cross-argument failures are caught during normal public API use, not only by offline test runners.
- [x] Runtime assertions are deterministic and every new promoted check is guarded by `should_run_asserts()`; they may be thorough when assertions are enabled.
- [x] Test-generation constraints remain the source for broad coverage and fixture-dependent discovery.
- [x] Promoted runtime helpers are covered by focused `testthat` tests and by at least one generated legal-combination case.
- [x] No private representation details or fixture-only assumptions leak into user-facing errors.

Completed outputs:

- `EDI/R/additional_asserts.R`
- Updated `EDI/DESCRIPTION`
- Updated `EDI/R/inference_all_abstract_non_param_boot.R`
- Updated `EDI/R/inference_all_abstract.R`
- Updated `EDI/R/design_fixed_abstract.R`
- Updated `EDI/R/design_seq_one_by_one_abstract.R`
- Updated `EDI/R/design_fixed_blocked_cluster.R`
- Updated `EDI/man/InferenceNonParamBootstrap.Rd`
- `EDI/tests/testthat/test-runtime-argument-combination-assertions.R`

Completion notes:

- Added unexported runtime helper assertions in `EDI/R/additional_asserts.R`.
- Promoted `min_number_usable_samples <= B` into `compute_bootstrap_two_sided_pval()` and `compute_bootstrap_confidence_interval()`.
- Consolidated model-formula variable availability checks into `assertFormulaContext()` and wired it into `Inference$initialize()`.
- Added `assertStrataClusterArgs()` for durable strata/cluster relationships, including constructor-time `strata_cols`/`cluster_col` overlap checks and data-time missing-column checks.
- Kept fixture-only, slow-path, optional-backend, and tier-specific registry constraints in the `package_tests/` runner rather than promoting them to runtime behavior.
- Every promoted runtime helper call is inside an existing `should_run_asserts()` block, and each helper also returns immediately when assertions are disabled.
- Existing response/design/inference compatibility checks were already enforced by `assertResponseType()`, `assertNoCensoring()`, and design-specific assertions; no broader compatibility rule was promoted beyond the stable cases listed above.
- Updated roxygen/Rd documentation for nonparametric bootstrap `min_number_usable_samples`: default is 5 and it must be less than or equal to `B`.

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
