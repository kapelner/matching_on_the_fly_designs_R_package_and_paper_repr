# Future Inference Hierarchy

Date: 2026-08-03

## Purpose

This document specifies the target `Inference*` architecture. It is not a
patch plan for the current inheritance ladder. The goal is a shallow hierarchy
where estimator identity is represented by inheritance and optional inference
algorithms are represented by checked components.

The architecture must make these properties true by construction:

1. A class exposes an operation only when it can execute that operation.
2. Each behavior has one component owner.
3. Each mutable field has one state owner.
4. Components are inherited through metadata and never composed twice.
5. Component requirements and collisions are validated when the class generator
   is defined.
6. Discovery reads metadata and pure compatibility predicates, never
   constructors or private method names.
7. Likelihood semantics are metadata, not superclass identity.
8. There are no legacy aliases. The package is unreleased, so old names should
   be deleted instead of carried forward.

## Architectural Rule

Inheritance answers one question:

```text
Can every child be substituted for this parent as the same kind of estimator?
```

Components answer a different question:

```text
Which optional algorithms and helper behaviors does this estimator expose?
```

Do not use inheritance to distribute randomization, bootstrap, jackknife, Wald,
likelihood-test, parametric-bootstrap, KK pass-through, GEE, GLMM, cache, or
count-likelihood behavior.

The target hierarchy is shallow:

```text
Inference
  -> optional estimator-family base with genuinely shared estimator math
      -> concrete estimator
```

Allowed parents:

- `Inference`: owns root data normalization, lifecycle, root state, and the
  minimal estimate contract.
- Family bases: allowed only when every descendant is substitutable for the
  base and the base does not expose optional inference algorithms.
- Concrete estimators: provide the estimator-specific point estimate,
  likelihood/objective hooks, simulation hooks, and compatibility metadata.

Disallowed parents:

- Any base whose main purpose is to add an optional algorithm.
- Any base whose descendants must opt out of public APIs with `supports_*`
  flags or throwing methods.
- Any likelihood-tier parent such as `NoLik`, `QuasiLik`, `PartialLik`, or
  `FullLik` if that parent carries public inference APIs.

## Class Metadata

Every canonical generator must have one immutable metadata record:

```r
list(
  name = "InferenceContinOLS",
  parent = "Inference",
  abstract = FALSE,
  exported = TRUE,
  response_types = "continuous",
  design_families = c("bernoulli", "complete_randomization"),
  compatibility = is_compatible_continuous_ols,
  likelihood_tier = "full",
  direct_components = c(
    "Wald",
    "LikelihoodTests",
    "NonparametricBootstrap",
    "ParametricLikelihoodBootstrap"
  ),
  required_packages = character()
)
```

Mandatory fields:

- `name`: canonical generator name.
- `parent`: canonical parent generator name or `NULL` for `Inference`.
- `abstract`: `TRUE` only for non-instantiable infrastructure.
- `exported`: `TRUE` only for user-facing estimator classes.
- `response_types`: response domains accepted by the estimator.
- `design_families`: randomization/design domains accepted by the estimator.
- `compatibility`: pure predicate over normalized design metadata.
- `likelihood_tier`: one of `none`, `quasi`, `partial`, or `full`.
- `direct_components`: components added by this class.
- `required_packages`: optional packages needed to instantiate or execute.

There is no alias metadata. If two names currently refer to the same generator,
choose the better canonical name and delete the other assignment. Because the
package is unreleased, compatibility aliases and deprecation schedules add
complexity without protecting users.

## Likelihood Tier

`likelihood_tier` is semantic metadata:

- `none`: no likelihood is defined for inference.
- `quasi`: estimating equation, quasi-likelihood, GEE, robust, or sandwich
  objective that is not a normalized model likelihood.
- `partial`: conditional or partial likelihood.
- `full`: full model likelihood.

Classification follows the implemented objective, not the class name. A class
named with `Lik` is not automatically likelihood-bearing. A likelihood tier
does not imply any public method. It only constrains which capabilities are
eligible.

Tier rules:

| Tier | Permitted likelihood capabilities |
|---|---|
| `none` | no likelihood ratio, score, gradient, Bartlett correction, or likelihood parametric bootstrap |
| `quasi` | Wald/sandwich inference; estimating-equation tests only under explicitly named capabilities |
| `partial` | score, gradient, likelihood ratio, and parametric calibration only when hooks exist |
| `full` | score, gradient, likelihood ratio, and parametric calibration only when hooks exist |

## Components

A component is a named, registered unit of optional behavior. It may provide
public methods, private helper methods, state ownership, host requirements,
capabilities, and dependencies on other components.

Expected component families:

- `RandomizationTest`
- `RandomizationCI`
- `NonparametricBootstrap`
- `RandomizationBootstrap`
- `BayesianBootstrap`
- `Jackknife`
- `Wald`
- `LikelihoodTests`
- `ParametricLikelihoodBootstrap`
- `StandardModelCache`
- `CountLikelihoodPlumbing`
- `KKPassThrough`
- `KKCompound`
- `KKGEE`
- `KKGLMM`
- `OffOptimumLikelihoodEval`
- `BartlettApproximation`

Empty scaffolds are not components. A registered component must either provide
behavior or be marked `status = "scaffold"` and forbidden from every effective
component set.

## Component Contract

Each component record must be executable metadata, not documentation:

```r
InferenceComponent(
  name = "ParametricLikelihoodBootstrap",
  status = "active",
  dependencies = "LikelihoodTests",
  public = list(...),
  private = list(...),
  provides_public_methods = c(...),
  provides_private_methods = c(...),
  owns_state = character(),
  requires_state = c(...),
  requires_public_methods = character(),
  requires_private_methods = c(
    "get_likelihood_test_spec",
    "simulate_under_lik_null"
  ),
  optional_public_methods = character(),
  optional_private_methods = character(),
  provides_capabilities = "parametric_likelihood_bootstrap",
  allowed_likelihood_tiers = c("partial", "full"),
  conflicts = character(),
  allowed_host_overrides = list(public = character(), private = character())
)
```

Contract rules:

1. `provides_public_methods` must exactly match names in the component public
   list.
2. `provides_private_methods` must exactly match names in the component
   private list.
3. `owns_state` contains only fields the component initializes and controls.
4. `requires_state` contains host/root fields the component reads or writes
   but does not own.
5. Every `private$foo`, `self$foo`, and `super$foo` reference in component
   bodies must be classified as provided, required, optional, or forbidden.
6. `super$foo` is allowed only when the contract names the parent service that
   supplies `foo`.
7. Components cannot redeclare root-owned state.
8. Components cannot collide with the parent, host, or another component unless
   the host declares an explicit validated override.
9. A component cannot provide capabilities that are invalid for the class
   likelihood tier.

## Component Inheritance

A class lists only the components it directly adds:

```text
effective_components(class)
  = effective_components(parent)
  + resolved_dependencies(direct_components(class))
  + direct_components(class)
```

Rules:

1. Parent components are inherited automatically.
2. A child must not re-list an inherited component.
3. Dependencies are resolved once in deterministic topological order.
4. Cycles are errors.
5. Duplicate direct or transitive components are errors.
6. Capability inheritance follows effective components.

This avoids duplicated declarations such as listing `RandomizationTest` or
`KKPassThrough` in every descendant once a parent already provides it.

## Capability Model

Capabilities are derived, not manually restated across descendants:

```text
effective_capabilities(class)
  = union(component_capabilities(effective_components(class)))
  + class_owned_capabilities(class)
```

Public optional APIs follow one rule:

```text
public method exists <=> capability exists <=> component or class hook exists
```

There are no public throwing stubs. There are no duplicate public/private
`supports_*()` methods. Callers use:

```r
obj$supports("likelihood_ratio")
obj$capabilities()
```

Those calls read immutable metadata.

Capability requirements live in one table:

```r
capability_requires = list(
  likelihood_ratio = list(
    likelihood_tier = c("partial", "full"),
    private_methods = "get_likelihood_test_spec"
  ),
  parametric_likelihood_bootstrap = list(
    likelihood_tier = c("partial", "full"),
    capabilities = "likelihood_ratio",
    private_methods = c(
      "get_likelihood_test_spec",
      "simulate_under_lik_null"
    )
  ),
  bayesian_bootstrap = list(
    private_methods = "compute_estimate_with_bootstrap_weights"
  ),
  nonparametric_bootstrap = list(
    private_methods = "compute_estimate",
    metadata = "resampling_policy"
  )
)
```

Public methods also live in one table:

```r
public_methods_for_capability = list(
  likelihood_ratio = c(...),
  parametric_likelihood_bootstrap = c(...),
  bayesian_bootstrap = c(...),
  nonparametric_bootstrap = c(...)
)
```

## Class Factory

All new or migrated generators should be created through one checked factory:

```r
InferenceContinOLS <- define_inference_class(
  name = "InferenceContinOLS",
  parent = Inference,
  metadata = list(
    abstract = FALSE,
    exported = TRUE,
    response_types = "continuous",
    design_families = c("bernoulli", "complete_randomization"),
    compatibility = is_compatible_continuous_ols,
    likelihood_tier = "full",
    required_packages = character()
  ),
  components = c(
    "Wald",
    "LikelihoodTests",
    "NonparametricBootstrap",
    "ParametricLikelihoodBootstrap"
  ),
  public = list(
    initialize = initialize_continuous_ols,
    compute_estimate = compute_estimate_continuous_ols
  ),
  private = list(
    get_likelihood_test_spec = ols_likelihood_test_spec,
    simulate_under_lik_null = simulate_ols_under_null
  ),
  overrides = list(public = character(), private = character()),
  lock_objects = FALSE
)
```

The factory must:

1. Register class metadata.
2. Resolve inherited and direct components.
3. Resolve component dependencies.
4. Validate component contracts against the parent, host, and likelihood tier.
5. Validate capability requirements.
6. Validate public API presence.
7. Reject duplicate components.
8. Reject method and state collisions without explicit overrides.
9. Assemble public and private lists in one place.
10. Return an ordinary R6 generator.

`utils::modifyList()` and list concatenation may be used inside the factory,
but raw component splicing outside the factory is forbidden.

## Discovery

`InferenceSuite` and similar discovery tools must select classes from
metadata:

```text
candidate if:
  abstract == FALSE
  exported == TRUE
  required_packages are available
  compatibility(normalized_design_metadata) == TRUE
```

Discovery must not:

- instantiate a class to see if construction fails,
- parse constructor error text,
- use a manual denylist of infrastructure classes,
- infer model family from class name,
- infer capability from private method names.

Unavailable packages should be reported separately from design
incompatibility.

## Source Invariants

Add structural tests that enumerate every `Inference*` generator and enforce:

1. Every canonical generator has exactly one metadata record.
2. Every metadata record points to a real generator.
3. Every generator declares `abstract`, `exported`, `likelihood_tier`,
   compatibility, and direct components.
4. Every `likelihood_tier` is one of `none`, `quasi`, `partial`, or `full`.
5. No source creates a top-level `Inference*` alias assignment to another
   `Inference*` generator.
6. Public optional methods equal
   `public_methods_for_capability(effective_capabilities(class))`.
7. No public optional method is an unconditional error.
8. No public/private method-name duplicate exists without an explicit override.
9. No class defines duplicate public/private `supports_*` hooks.
10. Every component is registered.
11. Every component is composed at most once effectively.
12. Every component dependency is declared and ordered once.
13. Every component body reference is classified by its contract.
14. Every required component hook and state field exists on the final host.
15. No component redeclares root-owned state.
16. No method or state collision is implicit.
17. No source contains `eval(body(Inference...))`.
18. No source reads `$private` from an R6 generator.
19. No semantic classification depends on `has_private_method()` or
    `object_has_private_method()`.
20. Discovery returns exactly compatible exported concrete canonical classes.

Behavioral tests must remain in place for estimates, standard errors,
confidence intervals, p-values, likelihood tests, bootstrap distributions, and
randomization distributions.

## Implementation TODOs

### Metadata Registry

- [x] Create `EDI_INFERENCE_CLASS_REGISTRY`.
- [x] Add `register_inference_class()`.
- [x] Add `get_inference_class_metadata()`.
- [x] Add `get_direct_components()`.
- [x] Add `get_effective_components()`.
- [x] Add `get_effective_capabilities()`.
- [x] Add metadata validators for mandatory fields and allowed `likelihood_tier` values.
- [x] Add tests that every generator has one canonical metadata entry.
- [x] Add tests that top-level `Inference*` alias assignments do not exist.

### Component Registry

- [x] Create `EDI_INFERENCE_COMPONENTS`.
- [x] Implement `InferenceComponent()`.
- [x] Register each active component with public/private lists, provided
  methods, owned state, required state, required hooks, optional hooks,
  dependencies, capabilities, conflicts, and allowed likelihood tiers.
- [x] Mark empty future work as `status = "scaffold"` and forbid scaffold
  components from effective component sets.
- [x] Add helpers for `component_public_names()` and
  `component_private_names()`.
- [x] Add tests that declared provided names match actual component list names.
- [x] Add parser-backed tests that every component body reference is declared
  as provided, required, optional, owned, or forbidden.

### Component Resolution

- [x] Implement `resolve_inference_components()`.
- [x] Make parent components inherited automatically.
- [x] Resolve dependencies in deterministic topological order.
- [x] Reject dependency cycles.
- [x] Reject direct re-listing of an inherited component.
- [x] Reject duplicate transitive components.
- [x] Add tests proving a daughter class inherits parent components without
  listing them again.
- [x] Add tests proving capabilities are inherited from effective components.

### Lazy Component Loading

- [ ] Split component metadata from component implementation so class
  discovery, capability queries, compatibility checks, and method availability
  are resolved without parsing heavy component bodies.
- [ ] Add `component_loader` metadata to `InferenceComponent()` for components
  whose public/private lists should be loaded on demand.
- [ ] Keep lightweight method stubs for capability-backed public APIs on the
  assembled class; each stub should load the component implementation and then
  dispatch to the real method on first use.
- [ ] Add a per-class component implementation cache so a component is loaded
  at most once per R session and repeated bootstrap calls do not pay parse or
  assembly cost repeatedly.
- [ ] Add `load_inference_component(component_name)` with deterministic errors
  for missing source files, missing optional packages, invalid component
  objects, and contract mismatches after load.
- [ ] Make lazy loading invisible to `supports()`, `capabilities()`,
  `InferenceSuite` discovery, and migration validation; these must continue to
  use eager metadata only.
- [ ] Classify components by load policy: `eager` for root contracts and cheap
  shared methods, `lazy` for expensive optional methods such as parametric
  likelihood bootstrap, Bartlett correction, simulation-heavy randomization
  bootstrap, Bayesian bootstrap, and GLMM/ordinal/count likelihood plumbing.
- [ ] Add tests that a class advertising `parametric_likelihood_bootstrap` does
  not load the parametric-bootstrap implementation until
  `compute_lik_ratio_bootstrap_*()` or diagnostics are called.
- [ ] Add tests that lazy loading preserves public method presence before load,
  produces identical results after load, and leaves unsupported methods absent.
- [ ] Add tests that lazy component dependencies are loaded in resolved
  topological order and that dependency cycles still fail before any
  implementation body is loaded.
- [ ] Add tests that optional-package failures occur only when invoking the lazy
  capability, not during package load, registry population, or class discovery.
- [ ] Benchmark package load and `InferenceSuite` discovery before and after
  lazy loading, with explicit targets for parse time and loaded object count.

#### Performance Gates

- [ ] Keep only a small static metadata table eager at package load.
- [ ] Do not scan, instantiate, or force every R6 generator during package
  load.
- [ ] Cache `get_effective_components()` and `get_effective_capabilities()` so
  repeated object construction and discovery do not recompute component
  closure.
- [ ] Run expensive contract validation only in tests, CI, or explicit
  development mode; normal package use should validate only cheap metadata
  invariants.
- [ ] Lazy-load heavy component implementations only on first use of the
  corresponding capability.
- [ ] Keep `InferenceSuite` discovery metadata-only, with no constructors and
  no component implementation loading.
- [ ] Add package-load and `InferenceSuite` discovery benchmarks as performance
  gates before enabling the shallow hierarchy by default.

#### Runtime Performance Traps

- [ ] Cache resolved dispatch targets after first use so bootstrap,
  likelihood, and randomization hot paths do not repeat capability lookup,
  component resolution, loader checks, and wrapper dispatch on every call.
- [ ] Cache effective components and effective capabilities at class scope, not
  object scope, so registry cost is not shifted from package load into every
  object construction or method call.
- [ ] Lazy-load by component or coherent component bundle, not by individual
  method, to avoid many tiny parse/load operations.
- [ ] Keep lazy stubs small: they should reference component names and method
  names, not capture component objects, full registries, parsed bodies, or large
  environments.
- [ ] Avoid dynamic R6 method rebinding on every object or call; use stable
  stubs that resolve once and store a function pointer or class-level dispatch
  cache.
- [ ] Keep contract completeness checks, parser-backed body-reference checks,
  collision audits, and full component validation out of production paths; run
  them only in tests, CI, or explicit development mode.
- [ ] Keep metadata tables lightweight: plain scalars, vectors, small lists,
  and cheap predicates only. Do not store method bodies, parsed expressions, or
  R6 generators in eager metadata.
- [ ] Assemble component public/private lists once per class, not during every
  object construction.
- [ ] Keep `InferenceSuite` compatibility predicates pure and cheap; do not
  normalize the same design repeatedly, probe optional packages repeatedly, or
  fit/touch models during discovery.
- [ ] Cache optional-package availability so repeated `requireNamespace()`
  checks across many classes do not dominate discovery or first-use cost.

### Class Definition

- [x] Implement `define_inference_class()`.
- [x] Route all component list assembly through the factory.
- [x] Implement `assemble_public()`.
- [x] Implement `assemble_private()`.
- [x] Validate required public methods, private methods, and state.
- [x] Validate component likelihood-tier eligibility.
- [x] Validate component capability requirements.
- [x] Validate public optional method presence from capabilities.
- [x] Reject method/state collisions unless listed in `overrides`.
- [x] Reject public/private name duplication unless listed in `overrides`.
- [x] Keep `lock_objects = FALSE` until the R6 tree is stable.

### Capability Tables

- [x] Create `capability_requires`.
- [x] Create `public_methods_for_capability`.
- [x] Add `supports(capability)` to the root class as a metadata query.
- [x] Add `capabilities()` to the root class as a metadata query.
- [x] Remove duplicated public/private `supports_*` hooks from the
  factory-defined architecture once metadata is authoritative.
- [x] Remove public optional throwing methods from factory-defined classes once
  method presence is capability-driven.
- [x] Add tests that `likelihood_tier = "none"` exposes no likelihood APIs.
- [x] Add tests that `likelihood_tier = "quasi"` exposes no likelihood-ratio
  API unless represented as a separate estimating-equation capability.
- [x] Add tests that `partial` and `full` tiers do not imply parametric
  bootstrap without a null simulator.

### Component Extraction

- [x] Extract `RandomizationTest`.
- [x] Extract `RandomizationCI`.
- [x] Extract `NonparametricBootstrap`.
- [x] Extract `RandomizationBootstrap`.
- [x] Extract `BayesianBootstrap`.
- [x] Extract `Jackknife`.
- [x] Extract `Wald`.
- [x] Extract `LikelihoodTests`.
- [x] Extract `ParametricLikelihoodBootstrap`.
- [x] Extract `StandardModelCache`.
- [x] Extract `CountLikelihoodPlumbing`.
- [x] Extract `KKPassThrough`.
- [x] Extract `KKCompound`.
- [x] Extract `KKGEE`.
- [x] Extract `KKGLMM`.
- [x] Extract `OffOptimumLikelihoodEval`.
- [x] Extract `BartlettApproximation` because
  `InferenceExtBartlettApprox` provides real private behavior.

### Shallow Hierarchy Migration

- [x] Add a manifest of current generators and reviewed target metadata.
- [x] Add a migration gate that rejects any class marked migrated while it still
  descends from an algorithmic compatibility base.
- [x] Add `mark_inference_class_migrated()` that only updates manifest status
  after validating current parent, effective components, effective
  capabilities, and public optional method presence.
- [x] Add a manifest summary test that reports pending counts by current parent,
  likelihood tier, response type, and KK/non-KK status.
- [x] Add a migration-order helper that lists leaf concrete classes before
  concrete parents so family migrations do not strand subclasses.
- [x] Add a strict opt-in test flag,
  `EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY`, that fails if any concrete class is
  still pending.

#### Migration Harness

- [x] Add reusable golden-test fixtures for continuous, incidence, count,
  proportion, ordinal, and survival completed designs.
- [x] Add a helper that compares old and new class outputs for estimate, SE,
  Wald/asymptotic CI, Wald/asymptotic p-value, likelihood tests, bootstrap,
  randomization, and jackknife methods when supported.
- [x] Add a method-availability snapshot helper that records public methods
  exposed before migration and compares them to capability-driven public
  methods after migration.
- [x] Add a private-state owner snapshot helper that detects new duplicated
  state owners after a class is moved to components.
- [x] Add a temporary dual-definition convention for one migrated class at a
  time, e.g. `InferenceFooLegacy` and `InferenceFoo`, only inside tests or
  migration fixtures; do not add package-level legacy aliases.
- [x] Require every family migration PR to include before/after manifest counts
  and the list of classes newly marked `migrated`.

#### Simple No-Likelihood Estimators

- [x] Identify all `likelihood_tier = "none"` concrete classes and split them
  into exact/randomization, jackknife/asymptotic, and pure estimator groups.
- [x] Extract `ExactTest` from `InferenceExact`.
- [x] Replace `is(inf_obj, "InferenceExact")` dispatch in simulation and
  package-test paths with `supports("exact_test")` or equivalent metadata-backed
  capability checks.

Current blocker: concrete exact classes still expose bootstrap, randomization,
Bayesian-bootstrap, and jackknife methods through the old `InferenceExact`
ancestry. Before moving them to `Inference + ExactTest`, each class needs an
explicit capability decision and golden comparison for the optional methods it
keeps or drops.

##### Exact Behavior Extraction

- [x] Add a manifest table for exact incidence classes with one row per class:
  `InferenceIncidExactBinomial`, `InferenceIncidExactFisher`, and
  `InferenceIncidenceExactZhang`.
- [x] For each exact incidence class, record current public optional methods
  inherited from `InferenceExact`, including exact, bootstrap, randomization,
  Bayesian-bootstrap, and jackknife methods.
- [x] For each exact incidence class, decide which inherited optional methods
  are intentional capabilities and which are legacy accidental surface area.
- [x] Add `ExactBinomialIncidence` component for matched-pair binomial behavior
  if `InferenceIncidExactBinomial` keeps behavior not shared by all exact tests.
- [x] Add `ExactFisherIncidence` component for Fisher/mantelhaen table-building
  behavior if `InferenceIncidExactFisher` keeps behavior not shared by all exact
  tests.
- [x] Add `ExactZhangIncidence` component for Zhang helper integration if
  `InferenceIncidenceExactZhang` keeps behavior not shared by all exact tests.
- [x] Decide whether `InferenceIncidCMH` is an exact-test component, a
  Wald/asymptotic blocked-incidence component, or both; do not group it with
  exact incidence classes until its public capabilities match that decision.
- [x] Add tests that exact-specific component public/private provided names
  match their actual list names.
- [x] Add parser-backed tests that exact-specific component references are
  declared as provided, required, optional, owned, or forbidden.

##### Custom Randomization Migration

- [x] Identify every concrete or extension class that currently inherits
  `InferenceRand` directly.
- [x] Split custom randomization hosts into public extension bases and concrete
  package estimators.
- [x] Add target metadata for each custom randomization host:
  `parent = "Inference"`, `components = "RandomizationTest"`, and explicit
  class-owned capabilities if any.
- [x] Add golden tests comparing current and migrated randomization p-values and
  randomization distributions on the continuous and incidence fixtures.
- [x] Convert randomization simulation and comprehensive-test dispatch from
  `is(obj, "InferenceRand")` / `is(obj, "InferenceRandCI")` to
  `supports("randomization_test")` / `supports("randomization_ci")` before
  changing production custom-randomization inheritance.
- [x] Replace direct `InferenceRand` inheritance in one custom randomization host
  at a time.
- [x] Add tests proving a migrated custom-randomization host does not expose
  randomization confidence-interval, randomization-bootstrap, nonparametric-
  bootstrap, Bayesian-bootstrap, or jackknife APIs unless the matching component
  is listed.
- [x] Remove any inherited randomization confidence-interval or bootstrap APIs
  from migrated custom-randomization hosts unless the matching component is
  listed explicitly.
- [x] Mark each custom-randomization class migrated only after method snapshots
  and golden randomization tests pass.

##### Exact Class Migration

- [ ] Add a temporary test-only legacy definition for
  `InferenceIncidExactBinomialLegacy` and compare it with the migrated
  `InferenceIncidExactBinomial`.
- [ ] Migrate `InferenceIncidExactBinomial` to `Inference` plus `ExactTest` and
  `ExactBinomialIncidence`.
- [ ] Preserve or intentionally remove inherited bootstrap, randomization,
  Bayesian-bootstrap, and jackknife APIs on `InferenceIncidExactBinomial`
  according to the exact capability manifest.
- [ ] Add golden tests for binomial exact estimate, exact p-value, exact CI, and
  any retained optional methods.
- [ ] Mark `InferenceIncidExactBinomial` migrated after golden tests and method
  snapshots pass.
- [ ] Add a temporary test-only legacy definition for
  `InferenceIncidExactFisherLegacy` and compare it with the migrated
  `InferenceIncidExactFisher`.
- [ ] Migrate `InferenceIncidExactFisher` to `Inference` plus `ExactTest` and
  `ExactFisherIncidence`.
- [ ] Preserve or intentionally remove inherited bootstrap, randomization,
  Bayesian-bootstrap, and jackknife APIs on `InferenceIncidExactFisher`
  according to the exact capability manifest.
- [ ] Add golden tests for Fisher exact estimate, exact p-value, exact CI, and
  any retained optional methods across iBCRD, blocking, and matching fixtures.
- [ ] Mark `InferenceIncidExactFisher` migrated after golden tests and method
  snapshots pass.
- [ ] Add a temporary test-only legacy definition for
  `InferenceIncidenceExactZhangLegacy` and compare it with the migrated
  `InferenceIncidenceExactZhang`.
- [ ] Migrate `InferenceIncidenceExactZhang` to `Inference` plus `ExactTest` and
  `ExactZhangIncidence`.
- [ ] Preserve or intentionally remove inherited bootstrap, randomization,
  Bayesian-bootstrap, and jackknife APIs on `InferenceIncidenceExactZhang`
  according to the exact capability manifest.
- [ ] Add golden tests for Zhang estimate, exact p-value, exact CI, and any
  retained optional methods on Bernoulli and matching-capable incidence
  fixtures.
- [ ] Mark `InferenceIncidenceExactZhang` migrated after golden tests and method
  snapshots pass.

##### Simple Estimator Migration

- [ ] Identify simple mean-difference no-likelihood classes and record their
  current direct parent, effective components, public methods, and private-state
  owners.
- [ ] Identify Wilcoxon/rank no-likelihood classes and record their current
  direct parent, effective components, public methods, and private-state owners.
- [ ] Decide whether each simple estimator keeps randomization, randomization
  CI, nonparametric bootstrap, randomization bootstrap, Bayesian bootstrap, and
  jackknife APIs.
- [ ] Migrate one simple mean-difference class to `Inference` plus only the
  components that match its retained capabilities.
- [ ] Add golden tests for estimate, randomization, bootstrap, Bayesian
  bootstrap, and jackknife outputs for that first simple mean-difference class.
- [ ] Repeat the simple mean-difference migration class by class, using the
  migration-order helper so leaf classes move before concrete parents.
- [ ] Migrate one Wilcoxon/rank class to `Inference` plus only the components
  that match its retained capabilities.
- [ ] Add golden tests for estimate, randomization, bootstrap, Bayesian
  bootstrap, and jackknife outputs for that first Wilcoxon/rank class.
- [ ] Repeat the Wilcoxon/rank migration class by class, using the
  migration-order helper so leaf classes move before concrete parents.

##### No-Likelihood Migration Marking

- [ ] Require every no-likelihood migration PR to include before/after manifest
  counts by no-likelihood group.
- [ ] Require every no-likelihood migration PR to list newly migrated classes and
  the optional methods intentionally kept or dropped for each class.
- [ ] Require `mark_inference_class_migrated()` to pass for each newly migrated
  no-likelihood class before checking off its class-specific migration TODO.
- [ ] Require golden output comparison to pass before marking any no-likelihood
  class `migrated`.
- [ ] Require method-availability snapshots to pass before marking any
  no-likelihood class `migrated`.
- [ ] Require private-state owner snapshots to pass before marking any
  no-likelihood class `migrated`.
- [ ] After all no-likelihood classes are migrated, delete no-longer-used
  algorithmic bases in this family and remove them from
  `EDI_INFERENCE_ALGORITHM_COMPATIBILITY_BASES`.

#### Quasi And Robust Estimators

- [ ] Identify all `likelihood_tier = "quasi"` concrete classes and verify
  whether each uses GEE, robust/sandwich, quasi-likelihood, or composite
  likelihood behavior.
- [ ] Extract `RobustSandwich` from robust regression and modified Poisson paths.
- [ ] Extract `CompositeLikelihoodTests` if composite likelihood needs public
  APIs distinct from normalized likelihood tests.
- [ ] Migrate `KKGEE` users to `Inference` plus `KKGEE` and required estimator
  components after host requirements are fully declared.
- [ ] Migrate count quasi-Poisson and robust Poisson classes to `Inference` plus
  `CountLikelihoodPlumbing` and quasi-specific components.
- [ ] Ensure quasi classes do not expose normalized likelihood-ratio capability
  unless represented by `estimating_equation_likelihood_ratio`.
- [ ] Mark migrated quasi/robust classes only after estimate, SE, CI, p-value,
  and method-availability snapshots match.

#### Partial-Likelihood Estimators

- [ ] Identify all `likelihood_tier = "partial"` concrete classes.
- [ ] Extract Cox/stratified-Cox shared behavior from
  `InferenceAsympLikStdModCache`.
- [ ] Extract conditional-logit shared behavior from current conditional
  incidence and ordinal classes.
- [ ] Extract LWA Cox and survival rank-regression shared behavior before moving
  their concrete classes.
- [ ] Migrate non-KK partial-likelihood classes to `Inference` plus
  `LikelihoodTests`, `StandardModelCache`, and family-specific components.
- [ ] Migrate KK partial-likelihood classes only after `KKPassThrough` and
  `KKCompound` host contracts pass collision and dependency validation.
- [ ] Verify partial-likelihood classes do not gain
  `parametric_likelihood_bootstrap` unless they provide a null simulator.

#### Full-Likelihood Estimators

- [ ] Identify all `likelihood_tier = "full"` concrete classes and split them by
  GLM, count, ordinal, incidence, proportion, survival, and KK/IVWC families.
- [ ] Extract GLM-family standard-model-cache behavior that is currently shared
  through `InferenceAsympLikStdModCache`.
- [ ] Extract count likelihood family behavior that is currently shared through
  `InferenceCountLikelihood`, `InferenceCountLikelihoodNoParamBootstrap`, and
  `InferenceCountCompositeLikelihood`.
- [ ] Extract zero-augmented count behavior from
  `InferenceCountZeroAugmentedPoissonAbstract`.
- [ ] Extract ordinal likelihood behavior for proportional odds, adjacent
  category, cloglog, cauchit, stereotype, continuation-ratio, and ordered probit
  paths.
- [ ] Extract incidence likelihood behavior for logit, probit, log-binomial,
  modified Poisson, binomial identity, and g-computation paths.
- [ ] Extract survival likelihood behavior for Weibull, dependent-censoring
  transform, Clayton copula, and frailty paths.
- [ ] Migrate full-likelihood classes to `Inference` plus
  `LikelihoodTests`, `StandardModelCache`, `ParametricLikelihoodBootstrap` when
  warranted, and family-specific components.
- [ ] Verify every migrated full-likelihood class has finite smoke tests for
  supported likelihood, bootstrap, and Bartlett paths.

#### KK And IVWC Estimators

- [ ] Finish declaring every `KKPassThrough`, `KKCompound`, `KKGEE`, and
  `KKGLMM` host requirement as required, optional, owned, or forbidden.
- [ ] Remove all direct `InferenceMixinKKPassThrough$public` and
  `InferenceMixinKKPassThrough$private` splices from concrete classes.
- [ ] Replace every `eval(body(InferenceMixinKKPassThrough$...))` usage with a
  named component override or helper.
- [ ] Migrate KK IVWC classes to `Inference` plus `KKPassThrough`,
  `KKCompound`, and estimator-specific components.
- [ ] Migrate KK one-likelihood classes to `Inference` plus `KKPassThrough`,
  `LikelihoodTests`, `ParametricLikelihoodBootstrap` when warranted, and
  estimator-specific likelihood components.
- [ ] Migrate KK GEE and GLMM classes after GEE/GLMM component contracts reject
  missing host hooks at class definition time.
- [ ] Add focused KK regression tests for matched-set weights, IVWC weighting,
  rank reduction, nonestimable fits, and block/cluster edge cases.

#### Base Deletion

- [ ] After a current algorithmic base has no concrete descendants, convert it
  into an internal component source or delete it.
- [ ] Delete `InferenceRand`, `InferenceRandCI`, `InferenceNonParamBootstrap`,
  `InferenceRandBootstrap`, `InferenceRandBootstrapCI`,
  `InferenceBayesianBootstrap`, `InferenceJackknife`, `InferenceExact`,
  `InferenceAsymp`, `InferenceAsympLik`, `InferenceParamBootstrap`,
  `InferenceAsympLikStdModCache`, count likelihood bases, and KK compound bases
  only after no concrete class inherits from them.
- [ ] Remove the base from `EDI_INFERENCE_ALGORITHM_COMPATIBILITY_BASES` in the
  same change that deletes or converts it.
- [ ] Enable `EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY` in normal tests once all
  concrete classes are migrated.
- [ ] Add the final strict test that no concrete class descends from an
  algorithmic compatibility base.

The manifest records 106 concrete generators as `pending` because they still
inherit through algorithmic compatibility bases. The final strict gate, `no
concrete class descends from an algorithmic compatibility base`, becomes
actionable after those pending records are drained family by family.

### Discovery

- [ ] Normalize design metadata in one helper.
- [ ] Replace suite constructor probing with metadata filtering.
- [ ] Replace manual infrastructure denylists with `abstract` and `exported`
  metadata.
- [ ] Add compatibility predicates for continuous, incidence, count,
  proportion, survival, ordinal, KK, non-KK, blocked, and uncensored-only
  designs.
- [ ] Add tests that constructor failures cannot influence discovery.
- [ ] Add tests that missing optional packages are reported separately from
  incompatibility.

### Static Cleanup

- [ ] Ban raw component splicing outside `define_inference_class()`.
- [ ] Ban `eval(body(Inference...))`.
- [ ] Ban `$private` reads from R6 generator symbols.
- [ ] Ban semantic classification through private method-name sniffing.
- [ ] Ban component redeclaration of root-owned state.
- [ ] Ban scaffold components in effective component sets.
- [ ] Ban implicit method and state collisions.

### Regression Gates

- [ ] Before migrating a family, add focused golden tests for estimates,
  standard errors, confidence intervals, p-values, and applicable bootstrap or
  randomization distributions.
- [ ] Add finite smoke tests for every class with
  `parametric_likelihood_bootstrap`.
- [ ] Add focused tests for count likelihood families, standard-model-cache
  families, KK pass-through families, KK compound families, and likelihood-test
  families.
- [ ] Run `Rscript fast_roxygenize.R` after exported API, class name,
  inheritance, or roxygen changes.
- [ ] Keep package load and targeted tests passing after each migrated family.

## Definition of Done

The hierarchy is complete when:

- Inheritance is shallow and represents substitutable estimator types only.
- Every class has valid metadata.
- Every component has an enforced contract.
- Effective components are inherited and composed once.
- Effective capabilities are derived from components and class-owned hooks.
- Public optional method presence exactly equals capabilities.
- Likelihood tier constrains capabilities but does not create APIs.
- Every component dependency, requirement, and collision is validated at class
  definition time.
- Every mutable field has one owner.
- Discovery is metadata-based and side-effect free.
- No legacy aliases exist.
- No optional API is disabled through private support flags or throwing stubs.
- No component depends on accidental inheritance details.
- No source uses implementation-body copying or generator-private access.
- Numerical regression tests pass for every migrated family.
