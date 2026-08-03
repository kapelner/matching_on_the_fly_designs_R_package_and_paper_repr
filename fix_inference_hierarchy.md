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

- [ ] Extract `RandomizationTest`.
- [ ] Extract `RandomizationCI`.
- [ ] Extract `NonparametricBootstrap`.
- [ ] Extract `RandomizationBootstrap`.
- [ ] Extract `BayesianBootstrap`.
- [ ] Extract `Jackknife`.
- [ ] Extract `Wald`.
- [ ] Extract `LikelihoodTests`.
- [ ] Extract `ParametricLikelihoodBootstrap`.
- [ ] Extract `StandardModelCache`.
- [ ] Extract `CountLikelihoodPlumbing`.
- [ ] Extract `KKPassThrough`.
- [ ] Extract `KKCompound`.
- [ ] Extract `KKGEE`.
- [ ] Extract `KKGLMM`.
- [ ] Extract `OffOptimumLikelihoodEval`.
- [ ] Extract `BartlettApproximation` only if it provides real behavior.

### Shallow Hierarchy Migration

- [ ] Add a manifest of current generators and reviewed target metadata.
- [ ] Move simple no-likelihood estimators to `Inference` plus components.
- [ ] Move GEE, robust, sandwich, and quasi-likelihood estimators to shallow
  family bases plus components.
- [ ] Move Cox, conditional, and other partial-likelihood estimators to shallow
  family bases plus components.
- [ ] Move full-likelihood GLM, count, ordinal, incidence, and survival
  estimators to shallow family bases plus components.
- [ ] Move KK and IVWC estimators after KK component contracts are complete.
- [ ] Delete any algorithmic base after it has no descendants.
- [ ] Add a test that no concrete class descends from an algorithmic
  compatibility base.

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
