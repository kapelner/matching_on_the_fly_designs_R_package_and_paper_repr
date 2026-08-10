# Spec: Professionally Comprehensive Test Suite

## Objective

Implement a professionally comprehensive test suite for the EDI package. The suite must cover public API behavior, legal argument combinations, statistical workflow paths, high-fan-in internal safety nets, runtime tiers, coverage reporting, and quality gates.

This spec is intentionally broader than `package_metadata/new_feature_plans/comprehensive_argument_checking.md`, but it depends on that spec.

## Strict Dependency

`package_metadata/new_feature_plans/comprehensive_argument_checking.md` is a hard prerequisite.

Do not begin implementation of this comprehensive suite spec until the argument-checking spec is fully implemented and accepted. In particular, the following must already exist and be working:

- checkmate-derived public argument contract extraction.
- public API inventory for exported functions/classes/public R6 methods.
- legal argument registry seeded from checkmate contracts.
- cross-argument constraint library.
- deterministic public fixtures.
- public argument-combination case generator.
- public argument-combination runner.
- public argument-combination analyzer.
- smoke/CI tiers for legal argument combinations.
- coverage reports that distinguish unidimensional argument checks from pairwise/higher-order legal-combination checks.

Reason: the comprehensive suite should build on a reliable legal-argument matrix rather than inventing a second, incompatible coverage model. Argument-combination coverage is the foundation for all later public API, fixture, workflow, and internal safety-net coverage.

## Non-Goals

- Do not duplicate the argument-combination implementation from `comprehensive_argument_checking.md`.
- Do not replace `package_tests/comprehensive_tests.R` in the first rollout.
- Do not make private/internal tests count as public contract coverage.
- Do not run unbounded Cartesian products.
- Do not add broad CI gates before the smoke and CI tiers are stable.

## Inputs From The Dependency

The completed argument-checking implementation must provide these inputs:

- `package_tests/public_api_inventory.csv`
- `package_tests/checkmate_argument_contracts.csv`
- `package_tests/public_argument_contract_registry.R`
- `package_tests/public_argument_combination_constraints.R`
- `package_tests/public_argument_combination_fixtures.R`
- `package_tests/public_argument_combination_cases.csv`
- `package_tests/public_argument_combination_results.csv`
- `package_tests/public_argument_combination_coverage.csv`
- `package_tests/public_argument_combination_failures.csv`

This comprehensive suite spec will consume those files rather than redefine them.

## Architecture

The comprehensive suite has seven layers:

1. Dependency intake and validation.
2. Existing comprehensive harness hardening.
3. Public API workflow coverage.
4. Statistical method-family coverage.
5. Internal high-fan-in safety-net coverage.
6. Unified coverage analysis.
7. Tiered CI/nightly/release gates.

## Artifacts

New or extended files under `package_tests/`:

- `validate_argument_checking_dependency.R`
- `comprehensive_suite_registry.R`
- `comprehensive_suite_fixtures.R`
- `comprehensive_suite_internal_surfaces.R`
- `run_comprehensive_suite.R`
- `analyze_comprehensive_suite.R`
- `comprehensive_suite_coverage.csv`
- `comprehensive_suite_failures.csv`
- `comprehensive_suite_exemptions.csv`
- `comprehensive_suite_report.html`

Existing files to integrate with, not replace initially:

- `package_tests/comprehensive_tests.R`
- `package_tests/analyze_comprehensive_tests.R`
- `package_tests/dedupe_comprehensive_tests_results.R`
- `package_tests/run_comprehensive_tests_sunk.R`

## Phase -1: Dependency Completion Gate

Goal: block comprehensive-suite implementation until `comprehensive_argument_checking.md` is complete.

TODO:

- [x] Run the smoke tier from `run_public_argument_combinations.R`.
- [x] Run the CI tier from `run_public_argument_combinations.R`.
- [x] Run `analyze_public_argument_combinations.R`.
- [x] Verify `public_api_inventory.csv` exists and includes exported R6 classes/functions.
- [x] Verify `checkmate_argument_contracts.csv` exists and is deterministic across two runs.
- [x] Verify `public_argument_combination_cases.csv` includes default, one-non-default, and pairwise cases.
- [x] Verify `public_argument_combination_coverage.csv` reports APIs that still only have unidimensional checks.
- [x] Verify CI-tier unexpected errors are either fixed or explicitly exempted.
- [x] Record accepted dependency commit/hash/date in this spec's implementation notes when implementation starts.

Acceptance criteria:

- Argument-combination smoke and CI tiers pass.
- Coverage and failure files are generated.
- The analyzer can identify unidimensional-only public APIs.
- No comprehensive-suite phase below starts before this phase is complete.

Implementation notes:

- Accepted dependency verification date: 2026-08-10.
- Accepted dependency commit: `c2fbd8158f4c12f1ef5a8b1073edc6bffb20cb82` (commit date `2026-08-10T11:26:28+03:00`).
- `run_public_argument_combinations.R smoke` passed with 10 aggregate result rows and 0 unexpected errors.
- `run_public_argument_combinations.R ci` passed with 29 aggregate result rows and 0 unexpected errors.
- `analyze_public_argument_combinations.R` generated 9,216 coverage rows, 0 failure rows, 2,809 drift rows, and 1,875 uncovered API rows.
- `public_api_inventory.csv` contains 9,216 rows: 121 exported functions, 128 exported R6 classes, and 8,967 public R6 method rows.
- `checkmate_argument_contracts.csv` was deterministic across two consecutive extractor runs; SHA-256 `df52b166be68a090c4799c9346cb051ecd4ef3622c6bf66fbc07bba0b47145b0`.
- `public_argument_combination_cases.csv` contains smoke and CI rows with default, one-non-default, pairwise, and targeted 3-way case kinds; SHA-256 `ea13e0c79f90c584fa8c474592d27e82c6e1aefd898007096454a409923b7016`.
- `public_argument_combination_ci_failure_summary.csv` reports 0 unexpected errors, 0 invariant failures, and `ci_should_fail = FALSE`.
- `check_public_argument_combination_quality_gates.R ci` passed with 23,391 quality-gate rows and 0 active hard gates.

## Phase 0: Baseline Integration Audit

Goal: map existing comprehensive tests and result files onto the completed public argument-combination inventory.

TODO:

- [x] Load `public_api_inventory.csv`.
- [x] Load `public_argument_combination_coverage.csv`.
- [x] Load current comprehensive result CSVs.
- [x] Identify public APIs covered by `package_tests/comprehensive_tests.R`.
- [x] Identify public APIs covered only by argument-combination tests.
- [x] Identify public APIs not covered by either.
- [x] Identify existing testthat files that directly test internals.
- [x] Classify internal tests as `keep`, `rewrite_public`, `internal_safety_net`, or `remove`.
- [x] Create initial `comprehensive_suite_exemptions.csv`.

Outputs:

- `package_tests/comprehensive_suite_baseline_audit.csv`
- `package_tests/comprehensive_suite_exemptions.csv`

Acceptance criteria:

- Every exported public API has one of: argument-combination coverage, comprehensive workflow coverage, focused testthat coverage, or exemption.
- Internal tests have an explicit classification.

Implementation notes:

- Added `package_tests/audit_comprehensive_suite_baseline.R` to reproduce the Phase 0 audit artifacts from the public API inventory, argument-combination coverage, current comprehensive result CSVs, and testthat sources.
- Generated `package_tests/comprehensive_suite_baseline_audit.csv` with 9,338 rows: 9,216 public API rows and 122 direct-internal testthat file rows.
- Public API coverage status counts: 7,274 `comprehensive_workflow`, 643 `focused_testthat`, 1,189 `multiple_sources`, and 110 `uncovered`.
- Generated `package_tests/comprehensive_suite_exemptions.csv` with 110 initial `phase0_uncovered_public_api` rows, so every currently uncovered public API has an explicit baseline exemption.
- Classified all 122 direct-internal testthat file rows as `internal_safety_net`; no internal testthat rows are unclassified.

## Phase 1: Suite Registry

Goal: define the unified registry for comprehensive suite coverage.

Implement `package_tests/comprehensive_suite_registry.R`.

Registry dimensions:

- public API surface.
- response type.
- design family.
- inference family.
- dataset/fixture family.
- method family.
- runtime tier.
- expected invariants.
- source runner: `argument_combinations`, `comprehensive_tests`, `testthat`, `internal_safety_net`.

TODO:

- [x] Define registry schema.
- [x] Import public API surfaces from the argument-checking inventory.
- [x] Import argument-combination coverage status.
- [x] Add method-family taxonomy:
  - estimate.
  - asymptotic.
  - exact.
  - bootstrap.
  - Bayesian bootstrap.
  - parametric bootstrap.
  - Bartlett.
  - jackknife.
  - randomization.
  - randomization-bootstrap.
  - simulation framework.
- [x] Add design-family taxonomy:
  - sequential.
  - fixed.
  - stratified.
  - clustered.
  - blocked.
  - matched.
  - optimal/search.
- [x] Add response-family taxonomy.
- [x] Add tier metadata.
- [x] Add exemption hooks with reasons and expiry dates.

Outputs:

- `package_tests/comprehensive_suite_registry.csv`

Acceptance criteria:

- Registry can answer "what must be tested for this exported public API?"
- Registry can answer "which runner is responsible for this coverage?"
- Registry distinguishes public contract coverage from internal safety-net coverage.

Implementation notes:

- Added `package_tests/comprehensive_suite_registry.R` to build the unified suite registry from `public_api_inventory.csv`, `public_argument_combination_coverage.csv`, `comprehensive_suite_baseline_audit.csv`, and `comprehensive_suite_exemptions.csv`.
- Generated `package_tests/comprehensive_suite_registry.csv` with 9,338 rows: 9,216 `public_contract` rows and 122 `internal_safety_net` rows.
- Registry schema includes public API surface, coverage scope, response/design/inference/dataset fixture/method families, runtime tier, required coverage, expected invariants, source runner, argument-combination coverage metrics, and exemption metadata.
- Source-runner counts: 7,274 `comprehensive_tests`, 643 `testthat`, 1,188 `comprehensive_tests;testthat`, 1 `argument_combinations;comprehensive_tests;testthat`, 110 `exemption`, and 122 `internal_safety_net`.
- Required-coverage counts: 3,332 `argument_combinations_or_workflow`, 255 `constructor_workflow`, 5,516 `public_workflow`, 3 `focused_public_test_or_workflow`, 110 `exemption`, and 122 `internal_safety_net`.
- Validation checks passed: 0 public rows with blank `required_coverage`, 0 public rows with blank `source_runner`, 0 uncovered public rows without exemption hooks, and 0 internal rows with a non-`internal_safety_net` runner.
- `comprehensive_suite_registry.csv` SHA-256: `ae2b6af181cd7762cba8977b6ffbf9c87607520b1c0eed475f5ee3d9814cd7e4`.

## Phase 2: Fixture Consolidation

Goal: reuse the argument-checking fixtures and extend them only where comprehensive workflow testing needs broader data.

Implement `package_tests/comprehensive_suite_fixtures.R`.

TODO:

- [ ] Import or wrap `public_argument_combination_fixtures.R`.
- [ ] Add fixture aliases compatible with existing `comprehensive_tests.R` labels.
- [ ] Add larger nightly/release fixtures where needed.
- [ ] Add survival fixtures with no/light/moderate censoring.
- [ ] Add count/proportion/incidence edge fixtures.
- [ ] Add ordinal fixtures with multiple level counts.
- [ ] Add clustered/blocked/matched design fixtures.
- [ ] Add fixture validation tests.
- [ ] Ensure every fixture is deterministic under fixed seed.

Acceptance criteria:

- Smoke fixtures are fast enough for CI.
- Nightly/release fixtures are marked with runtime tier.
- Existing comprehensive harness can either reuse these fixtures or map its datasets to equivalent fixture metadata.

## Phase 3: Existing Comprehensive Harness Hardening

Goal: make `package_tests/comprehensive_tests.R` easier to analyze and integrate without rewriting it wholesale.

TODO:

- [ ] Add stable case IDs to comprehensive result rows if missing.
- [ ] Add explicit `coverage_scope`.
- [ ] Add explicit `runner = "comprehensive_tests"`.
- [ ] Add explicit `github_commit_id` to comprehensive result rows for unambiguous test-run versioning.
- [ ] Add explicit `method_family`.
- [ ] Add explicit `argument_coverage_kind` where a call came from argument-combination inputs.
- [ ] Ensure `record_result()` emits enough fields for unified analysis.
- [ ] Preserve resumability of existing result files.
- [ ] Preserve existing filters for response/design/dataset/inference/function family.
- [ ] Keep existing slow-skip behavior, but make skip reasons analyzable.

Outputs:

- updated comprehensive result schema.
- migration notes for old result CSVs.

Acceptance criteria:

- Existing comprehensive commands still run.
- Existing result files remain readable.
- New rows can be joined with suite registry rows.

## Phase 4: Public API Workflow Coverage

Goal: ensure exported public APIs are exercised in realistic workflows, not only as isolated argument-combination calls.

TODO:

- [ ] For each exported concrete design class, create at least one complete workflow:
  - construct design.
  - add subjects.
  - assign treatment.
  - add responses.
  - retrieve public summaries/assignments/responses.
- [ ] For each exported concrete inference class, create at least one compatible completed design fixture.
- [ ] For `InferenceSuite`, test discovery and named `inference_params` forwarding through public constructors.
- [ ] For `SimulationFramework`, test a small complete simulation with at least one design and inference generator.
- [ ] For exported non-R6 public functions, classify as high-level, fitting wrapper, helper, or low-level exported routine and add workflow or direct public tests.
- [ ] Record workflow coverage separately from argument-combination coverage.

Outputs:

- workflow result rows in `comprehensive_suite_coverage.csv`.

Acceptance criteria:

- Every exported concrete R6 class has constructor coverage and at least one workflow or explicit exemption.
- High-priority exported functions have public workflow coverage.
- Workflow tests use public APIs unless explicitly marked internal safety-net.

## Phase 5: Statistical Method-Family Coverage

Goal: broaden method-family coverage while respecting runtime tiers.

TODO:

- [ ] Map all public inference methods to method families.
- [ ] Confirm argument-combination coverage exists for high-priority method arguments from the dependency.
- [ ] Add workflow calls for each method family by response/design/inference class where supported.
- [ ] Include representative p-value methods.
- [ ] Include representative confidence interval methods.
- [ ] Include representative estimate methods.
- [ ] Include debug distribution checks in nightly/release where useful.
- [ ] Add explicit non-estimable classification rules.
- [ ] Add support/unsupported classification based on public capability where available.

Method families:

- [ ] estimate.
- [ ] asymptotic Wald/score/likelihood-ratio/gradient.
- [ ] exact incidence/ordinal methods.
- [ ] nonparametric bootstrap.
- [ ] Bayesian bootstrap.
- [ ] parametric bootstrap.
- [ ] Bartlett correction.
- [ ] jackknife.
- [ ] randomization.
- [ ] randomization-bootstrap.
- [ ] m-out-of-n bootstrap/subsampling.

Acceptance criteria:

- Every high-priority inference family has smoke coverage.
- Every supported method family has CI or nightly coverage.
- Slow method families are not silently skipped; they are tiered or exempted.

## Phase 6: Internal High-Fan-In Safety Nets

Goal: directly test important shared internals without confusing that coverage with public API guarantees.

TODO:

- [ ] Build `comprehensive_suite_internal_surfaces.R`.
- [ ] Rank internals by fan-in, numerical risk, argument complexity, and historical failures.
- [ ] Identify validators/canonicalizers used by many public APIs.
- [ ] Identify shared resampling helpers.
- [ ] Identify model-matrix/data-shaping helpers.
- [ ] Identify numerical kernels with exported/public wrappers.
- [ ] Add direct tests only where public-only failures would be hard to diagnose.
- [ ] Link every internal surface to public entrypoints that use it.
- [ ] Mark coverage as `internal_safety_net`.

Acceptance criteria:

- No internal test exists without registry row and rationale.
- Internal coverage is reported separately.
- High-fan-in internal failures help localize public workflow failures.

## Phase 7: Unified Runner

Goal: provide one entry point that orchestrates the already-implemented argument-combination runner, the existing comprehensive harness, workflow checks, and internal safety nets.

Implement `package_tests/run_comprehensive_suite.R`.

Suggested CLI:

```sh
Rscript package_tests/run_comprehensive_suite.R smoke
Rscript package_tests/run_comprehensive_suite.R ci
Rscript package_tests/run_comprehensive_suite.R nightly
Rscript package_tests/run_comprehensive_suite.R release
Rscript package_tests/run_comprehensive_suite.R ci InferencePropBetaRegr
```

TODO:

- [ ] Validate Phase -1 dependency before running.
- [ ] Parse tier and optional filters.
- [ ] Run or import argument-combination results.
- [ ] Run selected comprehensive harness paths.
- [ ] Run workflow checks.
- [ ] Run internal safety-net checks for the selected tier.
- [ ] Enforce per-case timeouts.
- [ ] Preserve resumability.
- [ ] Emit unified result files.

Outputs:

- `package_tests/comprehensive_suite_results.csv`
- `package_tests/comprehensive_suite_failures.csv`

Acceptance criteria:

- Smoke tier runs end-to-end.
- CI tier can run without invoking release-scale comprehensive sweeps.
- Runner refuses to run if required argument-checking artifacts are missing.

## Phase 8: Unified Analyzer

Goal: report the full coverage picture across argument combinations, workflows, comprehensive harness paths, and internal safety nets.

Implement `package_tests/analyze_comprehensive_suite.R`.

TODO:

- [ ] Load argument-combination coverage.
- [ ] Load comprehensive suite registry.
- [ ] Load comprehensive harness results.
- [ ] Load workflow/internal result rows.
- [ ] Join all rows by public API, response type, design family, inference family, method family, tier, and coverage scope.
- [ ] Report uncovered public APIs.
- [ ] Report public APIs that only have argument-combination coverage but no workflow coverage.
- [ ] Report public APIs that only have workflow coverage but no multidimensional argument-combination coverage.
- [ ] Report high-fan-in internals lacking safety-net coverage.
- [ ] Report unsupported/skipped legal contexts by reason.
- [ ] Report slowest cases by tier.
- [ ] Emit HTML report if practical.

Outputs:

- `package_tests/comprehensive_suite_coverage.csv`
- `package_tests/comprehensive_suite_report.html`

Acceptance criteria:

- Analyzer can explain every uncovered API or point to an exemption.
- Analyzer distinguishes public contract, workflow, statistical method, and internal safety-net coverage.
- Analyzer imports argument-combination coverage rather than recomputing it independently.

## Phase 9: Runtime Tiers

Goal: make comprehensive coverage practical.

Tier definitions:

Smoke:

- [ ] dependency validation.
- [ ] argument-combination smoke results imported.
- [ ] one default workflow per major design/inference family.
- [ ] highest-fan-in internal validator smoke checks.
- [ ] tiny `B`/`r` values only.

CI:

- [ ] argument-combination CI results imported.
- [ ] pairwise legal argument coverage for high-priority public APIs is already satisfied by dependency.
- [ ] representative workflow coverage across response/design/inference families.
- [ ] selected method-family checks.
- [ ] strict public output invariants.

Nightly:

- [ ] broader comprehensive harness paths.
- [ ] more datasets/fixtures/formulas.
- [ ] targeted high-risk 3-way method-family workflows.
- [ ] broader internal safety-net checks.
- [ ] slow-case reporting.

Release:

- [ ] maximal registry-driven suite coverage.
- [ ] selected exhaustive small spaces.
- [ ] larger resampling counts.
- [ ] all required high-fan-in internal safety nets.
- [ ] final coverage report with exemptions.

Acceptance criteria:

- Smoke and CI tiers are stable.
- Nightly/release tiers are intentionally broader and may be scheduled or manual.
- No tier reimplements argument-combination generation.

## Phase 10: Quality Gates

Goal: enforce coverage without making the suite brittle.

TODO:

- [ ] Gate: argument-checking dependency must pass before this suite runs.
- [ ] Gate: no exported public API missing from registry unless exempted.
- [ ] Gate: no concrete exported R6 class without constructor coverage unless exempted.
- [ ] Gate: no high-priority public API without workflow coverage unless exempted.
- [ ] Gate: no public API with multiple configurable arguments unless already covered by argument-checking dependency or exempted.
- [ ] Gate: no internal safety-net test without rationale.
- [ ] Gate: no new CI-tier unexpected errors.
- [ ] Gate: no unknown skip/support status.
- [ ] Gate: no expired exemption.

Rollout:

1. Reporting-only gates.
2. Hard-fail dependency validation and schema validity.
3. Hard-fail smoke-tier unexpected errors.
4. Hard-fail missing high-priority coverage.
5. Hard-fail CI-tier unexpected errors.

Acceptance criteria:

- Gates are tier-aware.
- Exemptions are explicit and reviewed.
- Coverage ratchets upward without blocking unrelated work prematurely.

## Definition Of Done

- [ ] `comprehensive_argument_checking.md` is fully implemented and accepted before this spec begins.
- [ ] Argument-combination artifacts are validated as an input dependency.
- [ ] Every exported public API has inventory status, coverage status, or exemption.
- [ ] Every concrete exported R6 class has constructor coverage or exemption.
- [ ] Every high-priority public API has workflow coverage or exemption.
- [ ] Public APIs with multiple configurable arguments are covered by the argument-checking dependency or explicitly exempted.
- [ ] Statistical method-family coverage is tiered and reported.
- [ ] High-fan-in internals have safety-net coverage or explicit rationale.
- [ ] Public contract coverage and internal safety-net coverage are reported separately.
- [ ] Smoke and CI tiers run deterministically.
- [ ] Nightly and release tiers produce comprehensive reports.
- [ ] Existing comprehensive harness behavior remains available during rollout.
