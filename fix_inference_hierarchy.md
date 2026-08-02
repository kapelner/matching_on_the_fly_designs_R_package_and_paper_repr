# Fix the Inference Hierarchy

Date: 2026-08-02

## Purpose

This report audits the `Inference*` class hierarchy and the
`Inference*MixIn*` / `InferenceMixin*` components. It identifies current
correctness defects, structural duplication, capability mismatches, and a
target architecture that prevents those problems by construction.

The required end state is:

1. A concrete class exposes an inference operation if and only if it can
   execute that operation correctly.
2. Each behavior and each piece of mutable state has one owner.
3. A component is composed at most once into any class.
4. Component requirements, conflicts, and semantic eligibility are checked
   when the class generator is built, not discovered during a user call.
5. Class discovery uses explicit metadata and never probes constructors.
6. Refactoring preserves numerical behavior except for the live defects
   identified below.


## Executive Finding

The present hierarchy is not an is-a hierarchy. It is a feature-accumulation
ladder. Randomization, bootstrap, Bayesian bootstrap, jackknife, asymptotic,
likelihood, and parametric-bootstrap APIs are inherited in sequence. Concrete
classes near the bottom therefore acquire methods that are unrelated to their
statistical model and then try to disable them with private flags or throwing
overrides.

That design cannot reliably express capability. It has already produced
contradictory public and private support declarations, unreachable working
code, dead overrides, duplicated mixin composition, and discovery of internal
base classes as usable estimators.

The final design should use a shallow structural class hierarchy and explicit,
validated composition for orthogonal inference algorithms. Likelihood level is
a semantic classification, not a superclass. The superseded 2026-07-28 plan's
proposal to make `InferenceNoLik`, `InferenceQuasiLik`,
`InferencePartialLik`, and `InferenceFullLik` all inherit an API-bearing
`InferenceAsympLik` would retain the central defect: no-likelihood and
quasi-likelihood classes would still expose likelihood-ratio methods.

## Current Structural Problem

The main chain is approximately:

```text
Inference
  -> InferenceRand
  -> InferenceRandCI
  -> InferenceNonParamBootstrap
  -> InferenceRandBootstrap
  -> InferenceRandBootstrapCI
  -> InferenceBayesianBootstrap
  -> InferenceJackknife
  -> InferenceAsymp
  -> InferenceMLEorKMSummaryTable
  -> InferenceAsympLik
  -> InferenceParamBootstrap
```

Representative definitions are in:

- `EDI/R/inference_all_abstract_rand.R:6`
- `EDI/R/inference_all_abstract_rand_ci.R:6`
- `EDI/R/inference_all_abstract_non_param_boot.R:124`
- `EDI/R/inference_all_abstract_rand_bootstrap.R:72`
- `EDI/R/inference_all_abstract_rand_bootstrap_ci.R:27`
- `EDI/R/inference_all_abstract_bayesian_bootstrap.R:7`
- `EDI/R/inference_all_abstract_jackknife.R:7`
- `EDI/R/inference_all_abstract_asymp.R:6`
- `EDI/R/inference_all_abstract_mle_or_km_summary_table.R:8`
- `EDI/R/inference_all_abstract_asymp_lik.R:14`
- `EDI/R/inference_all_abstract_param_boot.R:26`

A static scan found 145 `R6Class` definitions. Of those, 143 descend from
`InferenceRand`, 140 from `InferenceNonParamBootstrap`, 137 from
`InferenceBayesianBootstrap`, 136 from `InferenceJackknife`, 103 from
`InferenceAsympLik`, and 64 from `InferenceParamBootstrap`. These counts show
that inheritance is being used to distribute optional algorithms rather than
to model substitutable types.

Examples of the resulting opt-out pattern include:

- `InferenceAllSimpleMeanDiff` inherits the parametric-bootstrap branch and
  then disables likelihood support in
  `EDI/R/inference_all_mean_diff.R:24` and `:260`.
- `InferenceAllSimpleWilcox` does the same in
  `EDI/R/inference_all_simple_wilcox.R:21` and `:288`.
- `InferenceCountCompositeLikelihood` inherits the branch, sets several
  support flags to false, and adds throwing overrides in
  `EDI/R/inference_count_composite_likelihood.R:7` and `:112-145`.

This also contradicts the stated intention in
`EDI/R/inference_all_abstract_param_boot.R:20-23` that unsupported model
families remain off the parametric-bootstrap branch.


## Contract and Ownership Defects

### Current registered mixin inventory

`EDI_MIXIN_CONTRACTS` currently registers seven components. The audit result
for each is:

| Component | Current structural status |
|---|---|
| `InferenceMixinCordeiroFerrariApprox` | Unused future scaffold with only a private marker; it provides no inference behavior |
| `InferenceMixinKKGEEShared` | Reusable behavior, but enforced host validation is missing |
| `InferenceMixinKKGLMMShared` | Reusable behavior with at least one undeclared private-method dependency |
| `InferenceMixinKKPassThrough` | Reusable behavior that depends on `super`, redeclares host state, and has undeclared host dependencies |
| `InferenceMixinKKPassThroughCompound` | Reusable specialization whose composition and collision rules are incomplete |
| `InferenceMixinLemonteGradientApprox` | Unused future scaffold with only a private marker; it provides no inference behavior |
| `InferenceMixinOffOptimumLikelihoodEval` | Genuine optional behavior, but support is detected from spec shape at runtime instead of guaranteed by composition |

At audit time, only `InferenceKKPassThroughCompound` appeared in
`EDI_MIXIN_COMPOSITIONS`. The immediate duplicate-composition fix also
registered `InferenceKKPassThroughCompoundNoParamBootstrap`, but the registry
still does not describe the full set of direct component-like list splices used
by the hierarchy.

Near-term TODOs for this inventory:

- Expand `EDI_MIXIN_CONTRACTS` from documentation-style metadata into an
  enforced component registry. Each entry should classify `status` (`active`
  or `scaffold`), provided public/private methods, required public/private
  methods, required state, dependencies, conflicts, provided capabilities, and
  allowed host overrides.
- Delete `InferenceMixinCordeiroFerrariApprox` and
  `InferenceMixinLemonteGradientApprox` until they provide behavior, or mark
  them explicitly as `status = "scaffold"` and enforce that scaffold mixins are
  not used in production class composition.
- Route every direct mixin splice through `compose_inference_mixins()` or its
  successor. Raw `utils::modifyList(... InferenceMixin*$public/private ...)`
  host composition should be treated as a structural violation.
- Complete the known contract gaps: `InferenceMixinKKGLMMShared` must declare
  `clear_nonestimable_state`; `InferenceMixinKKPassThrough` must declare its
  host state requirements such as `y`, `w`, `X`, and `dead`, its
  `supports_likelihood_tests` dependency, and its dependence on a superclass
  bootstrap implementation.
- Keep `InferenceMixinKKPassThroughCompound` dependent on
  `InferenceMixinKKPassThrough` through explicit component metadata, not by
  re-splicing pass-through in descendants.
- Give `InferenceMixinOffOptimumLikelihoodEval` an explicit marker capability
  such as `off_optimum_likelihood_eval`; support must be determined by
  composition metadata, not inferred from likelihood-spec shape.
- Add host-validation tests for every mixin-using class: required host methods
  and state exist, optional references are classified, no component appears
  twice transitively, and method/state collisions occur only through declared
  overrides.
- Add static tests that parse mixin source bodies for `private$`, `self$`, and
  `super$` references and require each reference to be provided by the mixin,
  inherited from a declared base, or listed as a host requirement.

### Incomplete mixin contracts

The registry in `EDI/R/mixin_contracts.R:24-78` is documentation-like metadata,
not an enforced complete contract.

Confirmed omissions include:

- `InferenceMixinKKGLMMShared` calls
  `private$clear_nonestimable_state()` at
  `EDI/R/inference_mixin_kk_glmm_shared.R:152`, but the method is absent from
  its registered requirements.
- `InferenceMixinKKPassThrough` reads host-owned `private$y` at
  `EDI/R/inference_mixin_kk_passthrough.R:45` and calls
  `private$supports_likelihood_tests()` around lines 268-282. Neither dependency
  is declared in its contract.
- Optional references such as `compute_fast_bootstrap_distr` and
  `compute_weighted_estimate_ivwc` are not classified as optional.

The tests at `EDI/tests/testthat/test-mixin-contracts.R:11-31` validate registry
shape and collation but not contract completeness or whether a host satisfies
the contract. The collision test at lines 33-43 covers only registered mixin
compositions.

Implementation TODOs:

- Extend each `EDI_MIXIN_CONTRACTS` entry with separate required/provided and
  optional surfaces: `requires_private_methods`, `requires_private_state`,
  `requires_public_methods`, `provides_private_methods`,
  `provides_private_state`, `provides_public_methods`,
  `optional_private_methods`, and `optional_public_methods`.
- Populate those fields from the actual mixin bodies before changing tests.
  The immediate known corrections are:
  `InferenceMixinKKGLMMShared` requires `clear_nonestimable_state`;
  `InferenceMixinKKPassThrough` requires host state `y`, `w`, `X`, and `dead`;
  `InferenceMixinKKPassThrough` requires `supports_likelihood_tests`; and
  `compute_fast_bootstrap_distr` plus `compute_weighted_estimate_ivwc` are
  optional private hooks.
- Add a parser-based contract audit that scans each `InferenceMixin*` source
  body for `private$foo`, `self$foo`, and `super$foo` references. Each
  reference must be classified as provided by the mixin, required from the
  host/parent, optional, or explicitly ignored with a reason.
- Add host-satisfaction tests for every class that composes a mixin. For each
  host, combine its parent chain, directly supplied public/private lists, and
  composed mixins, then assert every required method and state field is
  available before package load can silently rely on it.
- Teach `compose_inference_mixins()` or its successor to validate required
  methods/state at class-definition time, not only in tests. A class with an
  unsatisfied component contract should fail while the generator is built.
- Replace the existing broad shape-only tests in
  `test-mixin-contracts.R:11-31` with invariant tests that compare declared
  `provides_*` names to the actual mixin list names and declared
  `requires_*`/`optional_*` names to the parsed references.
- Expand the collision test beyond `EDI_MIXIN_COMPOSITIONS` so direct
  one-mixin hosts are also checked for host/mixin method collisions and state
  redeclarations.
- Preserve a narrow allowlist for intentional references to `super`, but make
  every allowed `super$foo` call name the parent contract it depends on.
- Fail CI when a mixin gains a new `private$`, `self$`, or `super$` reference
  without updating the contract metadata in the same change.


### Duplicated state ownership

`InferenceMixinKKPassThrough$private` declares fields such as `m`, `y_temp`,
`dead`, `w`, `X`, `any_censoring`, `optimization_alg`, and `cached_mod`
around `EDI/R/inference_mixin_kk_passthrough.R:253-266`. These are already
root-owned fields in `EDI/R/inference_all_abstract.R`.

A mixin must not redeclare host state. Duplicate fields make initialization,
cache invalidation, and ownership ambiguous. Components should declare
`requires_state`; only the root or one designated state component may declare
`owns_state`.

### Behavior copied through function bodies

Code such as `eval(body(InferenceMixin...))` appears in the KK hierarchy,
including `EDI/R/inference_continuous_KK_ols_ivwc.R:35` and
`EDI/R/inference_all_abstract_KK_passthrough_compound.R:27`. This copies an
implementation without a callable contract, depends on lexical context, and
drifts when the source function changes.

Replace these copies with named shared functions or a component-provided
method. The final tree must contain no `eval(body(Inference...))` calls.

### Mixins depend on inheritance details

`InferenceMixinKKPassThrough` calls `super` in its methods near
`EDI/R/inference_mixin_kk_passthrough.R:35`. A component should not depend on
the accidental parent selected by its host. Its dependencies must be named
services or hooks supplied by the host/component contract.

## Capability Detection Defect

`InferenceRandCI` infers model semantics by testing whether private methods
named `is_a_*` exist in `EDI/R/inference_all_abstract_rand_ci.R:127-136`. It
even treats the presence of any KK pass-through marker as GLM-like.

`has_private_method()` in `EDI/R/inference_all_abstract.R:752-756` checks only
whether the name exists, not whether the method returns true or represents the
required statistical property. Private implementation details are therefore
being used as runtime type tags.

Replace all method-name sniffing with immutable class metadata. No algorithm
should infer a statistical family from the existence of a private helper.

## Aliases and Discovery

The audit found two compatibility aliases:

- `InferenceAllKKCompoundMeanDiff <- InferenceAllKKMeanDiffIVWC` in
  `EDI/R/inference_all_KK_mean_diff_IVWC.R:298`.
- `InferenceIncidExactZhang = InferenceIncidenceExactZhang` in
  `EDI/R/inference_incidence_exact_zhang.R:104`.

Aliases must be represented as alias metadata, not scanned as independent
semantic classes. They should be deprecated on a published schedule and
excluded from `InferenceSuite` results.

## Target Architecture

### Structural inheritance

Keep inheritance shallow:

```text
Inference
  -> optional response/family base containing genuinely shared estimator math
      -> concrete estimator
```

Allowed inheritance rules:

1. `Inference` owns universal data normalization, root state, lifecycle, and
   the core estimate contract.
2. A family base is allowed only when every descendant is substitutable for it
   and the base does not expose optional inference algorithms.
3. Randomization, intervals, bootstrap variants, jackknife, Wald inference,
   likelihood tests, and parametric bootstrap are components, not ancestors.
4. Likelihood tier is required class metadata, not inheritance.
5. A deprecated compatibility base may forward temporarily during migration,
   but no new concrete class may inherit it.

This avoids both the current feature ladder and a combinatorial collection of
base classes for every possible capability set.

### Required class metadata

Every generator must declare immutable metadata with at least:

```r
list(
  name = "InferenceContinOLS",
  abstract = FALSE,
  exported = TRUE,
  alias_of = NULL,
  response_types = "continuous",
  design_families = c("bernoulli", "complete_randomization"),
  compatibility = is_compatible_continuous_ols,
  likelihood_tier = "full",
  capabilities = c("wald", "likelihood_ratio", "score",
                   "nonparametric_bootstrap"),
  components = c("wald", "likelihood_tests", "nonparametric_bootstrap"),
  required_packages = character()
)
```

`likelihood_tier` must be exactly one of:

- `none`: the implementation does not define a likelihood used for inference.
- `quasi`: an estimating equation, quasi-likelihood, GEE, sandwich, or robust
  objective that is not a normalized model likelihood.
- `partial`: a conditional or partial likelihood.
- `full`: a full model likelihood.

Names such as `OneLik` are not evidence of a likelihood tier. Classification
must follow the implemented estimand and objective. A robust, quantile, IVWC,
or matching estimator remains `none` or `quasi` even if a wrapper name contains
`Lik`. Conditional logit and Cox estimators are `partial`; conventional GLM and
parametric survival/count models are `full` when their implemented objective is
the actual model likelihood.

### Capability rules

The following rules are mandatory:

| Tier | Permitted likelihood capabilities |
|---|---|
| `none` | No likelihood ratio, score, likelihood gradient, Bartlett correction, or likelihood parametric bootstrap |
| `quasi` | Wald/sandwich inference; a separately named estimating-equation score only when explicitly implemented |
| `partial` | Score, gradient, and likelihood ratio when implemented; parametric calibration only with an explicit valid null simulator |
| `full` | Score, gradient, and likelihood ratio when implemented; parametric calibration only with an explicit valid null simulator |

Observed information, Fisher information, Bartlett correction, and
off-optimum likelihood evaluation are independent capabilities. A likelihood
tier does not imply them automatically.

Additional implications:

- `parametric_likelihood_bootstrap` requires tier `partial` or `full`, a
  likelihood-test specification, and a null simulator.
- `bayesian_bootstrap` requires a weighted-estimate hook.
- `nonparametric_bootstrap` requires an estimator hook and a declared design
  resampling policy.
- `randomization` requires a statistic and a design capable of drawing valid
  assignments.
- `jackknife` requires an explicit deletion unit and duplicate/cluster policy.
- `likelihood_ratio` requires defined null and alternative objectives on the
  same likelihood scale.

### Public API rule

For every optional operation:

```text
public method exists <=> capability is declared <=> component is installed
```

There must be no public throwing stubs and no public optional method whose
private support flag is false. Generic callers can use one stable query:

```r
obj$supports("likelihood_ratio")
obj$capabilities()
```

The query reads immutable generator metadata. There are no duplicate public
and private `supports_*()` methods.

### Components

The expected orthogonal components are:

- `InferenceMixinRandomization`
- `InferenceMixinRandomizationCI`
- `InferenceMixinNonparametricBootstrap`
- `InferenceMixinRandomizationBootstrap`
- `InferenceMixinBayesianBootstrap`
- `InferenceMixinJackknife`
- `InferenceMixinWald`
- `InferenceMixinLikelihoodTests`
- `InferenceMixinParametricLikelihoodBootstrap`
- `InferenceMixinStandardModelCache`
- `InferenceMixinCountLikelihoodPlumbing`
- `InferenceMixinKKPassThrough`
- `InferenceMixinKKCompound`
- `InferenceMixinGEE`
- `InferenceMixinGLMM`

The exact names can follow repository convention, but each concept must have a
single implementation owner. Empty approximation scaffolds such as the current
Cordeiro-Ferrari and Lemonte gradient mixins should either acquire a real
contract and implementation or be removed.

### Component contract

Replace the current registry entries with executable specifications:

```r
InferenceComponent(
  name,
  dependencies = character(),
  public = list(),
  private = list(),
  owns_state = character(),
  requires_state = character(),
  requires_methods = character(),
  optional_methods = character(),
  provides_capabilities = character(),
  allowed_likelihood_tiers = character(),
  conflicts = character()
)
```

A single `define_inference_class()` factory should:

1. Resolve component dependencies topologically.
2. Reject dependency cycles.
3. Reject duplicate direct or transitive components.
4. Reject public/public, private/private, public/private, and state collisions.
5. Allow a collision only through a named, validated `overrides` declaration.
6. Validate required methods and state against the parent, host, and components.
7. Validate capabilities against likelihood tier and supplied hooks.
8. Validate that declared capabilities equal component-provided capabilities.
9. Attach immutable metadata to the R6 generator.
10. Refuse construction of an abstract generator.

Target usage, shown as architecture pseudocode:

```r
InferenceContinOLS <- define_inference_class(
  name = "InferenceContinOLS",
  parent = Inference,
  metadata = list(
    abstract = FALSE,
    exported = TRUE,
    response_types = "continuous",
    likelihood_tier = "full",
    capabilities = c(
      "wald", "likelihood_ratio", "score",
      "nonparametric_bootstrap", "parametric_likelihood_bootstrap"
    )
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
  )
)
```

The factory is preferable to repeated `modifyList()` or `c()` splicing because
R6 has single inheritance and no native checked trait system. It produces an
ordinary R6 generator after validating the composition once.

## Implementable Component Migration Spec

The migration must be executable without forcing every concrete class to list
every inherited capability. Components are inherited through class metadata in
the same direction as R6 inheritance:

```text
effective_components(class)
  = effective_components(parent)
  + direct_components(class)
  + resolved_dependencies(direct_components(class))
```

A class lists only the components it adds or intentionally overrides. For
example, if `InferenceRand` temporarily provides `RandomizationTest`, then a
daughter class does not list `RandomizationTest` again. Duplicate direct or
transitive components are errors unless a migration bridge explicitly declares
the duplicate as a temporary no-op alias.

### Component Metadata Store

Do not rely on mutating R6 generator environments as the first implementation
step. Add an ordinary package-level registry keyed by canonical class name:

```r
EDI_INFERENCE_CLASS_REGISTRY = list(
  InferenceRand = list(
    abstract = TRUE,
    exported = FALSE,
    parent = "Inference",
    likelihood_tier = "none",
    direct_components = "RandomizationTest",
    aliases = character(),
    response_types = "all",
    compatibility = "all_designs"
  )
)
```

Once this registry is stable, attaching the resolved metadata to each R6
generator is optional implementation detail. The source of truth should remain
one queryable metadata object.

TODOs:

- Create `EDI_INFERENCE_CLASS_REGISTRY` in a new collated file such as
  `inference_class_registry.R`, after component definitions and before
  concrete inference classes.
- Add `register_inference_class(name, parent, metadata, direct_components,
  aliases = character())` and use it to populate the registry.
- Add `get_inference_class_metadata(name)`,
  `get_direct_components(name)`, `get_effective_components(name)`,
  `get_effective_capabilities(name)`, and
  `get_inference_alias_target(name)`.
- Make `get_effective_components()` walk the registered parent chain, resolve
  component dependencies, and fail on duplicate component names.
- Add tests that prove `InferenceRand` can provide a randomization component
  once and that a daughter class inherits it without re-listing it.
- Add tests that a daughter re-listing an inherited component fails unless it
  declares a temporary migration override.

### Component Registry Shape

Replace the current loose `EDI_MIXIN_CONTRACTS` entries with component records
that can be validated mechanically:

```r
EDI_INFERENCE_COMPONENTS = list(
  RandomizationTest = list(
    status = "active",
    public = InferenceMixinRandomization$public,
    private = InferenceMixinRandomization$private,
    provides_capabilities = "randomization",
    dependencies = character(),
    allowed_likelihood_tiers = c("none", "quasi", "partial", "full"),
    owns_state = character(),
    requires_state = c("y", "w", "n"),
    requires_private_methods = c("run_randomization_iteration"),
    optional_private_methods = c("compute_fast_randomization_distr"),
    requires_public_methods = character(),
    optional_public_methods = character(),
    allowed_host_overrides = list(public = character(), private = character()),
    allowed_collisions = list(public = character(), private = character())
  )
)
```

TODOs:

- Rename `EDI_MIXIN_CONTRACTS` to `EDI_INFERENCE_COMPONENTS` only after all
  current tests use the new schema. Until then, add the new fields alongside
  the old fields.
- Add `component_public_names(component)` and
  `component_private_names(component)` helpers that read the actual list names
  from each component object.
- Add invariant tests that `provides_public_methods` and
  `provides_private_methods` exactly match the actual component list names.
- Add invariant tests that every `requires_*` and `optional_*` reference is a
  real symbol referenced by that component or explicitly marked as a semantic
  host requirement.
- Add a `status = "scaffold"` rule: scaffold components may be registered but
  cannot appear in any class's effective component set.

### Factory Behavior

`define_inference_class()` is the endpoint. During migration, it can be
introduced as a thin wrapper over `R6::R6Class()` while still emitting ordinary
R6 generators:

```r
define_inference_class = function(
  name,
  parent,
  metadata,
  components = character(),
  public = list(),
  private = list(),
  overrides = list(public = character(), private = character()),
  aliases = character(),
  lock_objects = FALSE
) {
  register_inference_class(name, parent$classname, metadata, components, aliases)
  resolved = resolve_inference_components(name, components, parent$classname)
  validate_component_contract(name, parent, resolved, public, private, overrides)
  R6::R6Class(
    classname = name,
    inherit = parent,
    public = assemble_public(resolved, public, overrides),
    private = assemble_private(resolved, private, overrides),
    lock_objects = lock_objects
  )
}
```

The first production version may keep `register_inference_class()` separate
from generator creation if that reduces risk. The required behavior is the
same: no class may silently acquire a component or public optional method.

TODOs:

- Implement `resolve_inference_components()` with deterministic ordering:
  inherited effective components first, then dependencies, then direct
  components.
- Implement `validate_component_contract()` with separate checks for required
  state, required private methods, required public methods, optional hooks,
  likelihood-tier compatibility, required packages, and declared overrides.
- Implement `assemble_public()` and `assemble_private()` using
  `utils::modifyList()` only inside the factory. Ban raw component
  `modifyList()` calls outside the factory after migration starts.
- Add an escape hatch `migration_allow_legacy = TRUE` only for temporary
  compatibility generators. It must require a deletion target phase in a
  comment and a test that counts remaining uses.
- Keep `lock_objects = FALSE` as the default until the whole R6 tree is stable.

### Compatibility-Bridge Recipe

Algorithmic bases should not disappear in one commit. Convert each old base
into a bridge in this order:

1. Register the old base as `abstract = TRUE`.
2. Register the component that represents the behavior currently owned by the
   base.
3. Make the old base's metadata include that component as a direct component.
4. Leave the old public/private methods byte-equivalent at first.
5. Add tests that a child inherits the effective component from the old base
   without re-listing it.
6. Convert concrete children in batches to `define_inference_class()`.
7. When no concrete class inherits the old base, delete the bridge.

Example:

```r
InferenceRand = define_inference_class(
  name = "InferenceRand",
  parent = Inference,
  metadata = list(abstract = TRUE, likelihood_tier = "none"),
  components = "RandomizationTest",
  public = legacy_rand_public,
  private = legacy_rand_private,
  migration_allow_legacy = TRUE
)

InferenceContinOLS = define_inference_class(
  name = "InferenceContinOLS",
  parent = InferenceRand,
  metadata = list(
    abstract = FALSE,
    exported = TRUE,
    response_types = "continuous",
    likelihood_tier = "full"
  ),
  components = c("Wald", "LikelihoodTests", "ParametricLikelihoodBootstrap"),
  public = contin_ols_public,
  private = contin_ols_private
)
```

`InferenceContinOLS` inherits `RandomizationTest` through `InferenceRand`; it
does not list it again.

TODOs:

- Bridge `InferenceRand` first and prove inherited randomization metadata works.
- Bridge `InferenceRandCI` only for randomization-CI behavior, not for
  bootstrap or Wald behavior.
- Bridge `InferenceNonParamBootstrap`, `InferenceRandBootstrap`, and
  `InferenceRandBootstrapCI` as separate bootstrap components.
- Bridge `InferenceBayesianBootstrap`, `InferenceJackknife`, and
  `InferenceAsymp` one at a time.
- Bridge likelihood and parametric-bootstrap bases only after
  `likelihood_tier` metadata is enforced.

### Concrete-Class Migration Recipe

Each concrete class should follow the same mechanical sequence:

1. Add a registry entry with current parent, abstract/exported status,
   response types, design compatibility predicate, required packages,
   likelihood tier, aliases, and current effective capabilities.
2. Add a characterization test that compares current public optional method
   names and current behavior for a small representative design.
3. Convert only direct list splices in that class to component references.
4. Remove duplicated components that are inherited from the parent.
5. Move throwing optional methods out of the class unless a component declares
   the capability.
6. Run the focused class/family tests.
7. Compare the manifest before and after. Differences require a defect number
   and regression test.

TODOs:

- Start with `InferenceAllSimpleMeanDiff` and `InferenceAllSimpleWilcox`
  because they expose the current opt-out pattern clearly.
- Migrate `InferenceContinOLS` next because it exercises full-likelihood,
  Wald, likelihood-test, bootstrap, and parametric-bootstrap interactions.
- Migrate one KK IVWC class only after `KKPassThrough` and `KKCompound`
  components have complete contracts and inherited component resolution works.

### Capability Metadata Rules

Capabilities are derived from effective components plus class-owned hooks. They
are not manually listed in every daughter class unless the class adds a new
component or intentionally disables a temporary bridge capability.

```text
effective_capabilities(class)
  = union(component_capabilities(effective_components(class)))
  + class_declared_capabilities
  - deprecated_bridge_suppressions
```

TODOs:

- Add `capability_requires` rules in one table, for example
  `parametric_likelihood_bootstrap` requires likelihood tier `partial` or
  `full`, `get_likelihood_test_spec`, and `simulate_under_lik_null`.
- Add `public_methods_for_capability` rules in one table, for example
  `likelihood_ratio` maps to the public LR p-value/CI methods.
- Add tests that public optional method presence equals
  `public_methods_for_capability(effective_capabilities(class))`.
- Delete public `supports_*` and private `supports_*` hooks after callers use
  `obj$supports()` backed by metadata.
- Keep temporary suppressions only in the registry, with a removal phase and a
  test counting remaining suppressions.

## Family Classification Guidance

The migration inventory must classify actual implementations, not infer type
from names. The following groups are starting points and require confirmation
against each estimator's objective:

| Tier | Typical repository families |
|---|---|
| `none` | Exact/randomization procedures, rank and Wilcoxon methods, quantile estimators, matching/IVWC summaries, g-computation summaries, KM/log-rank/RMST, Bai-style estimators |
| `quasi` | GEE, quasi-Poisson, modified Poisson with robust variance, sandwich marginal models, robust estimating equations |
| `partial` | Cox, stratified Cox, conditional logistic, conditional adjacent-category, LWA Cox |
| `full` | Gaussian OLS likelihood, logistic/probit/binomial GLMs, Poisson/negative-binomial/hurdle/zero-inflated models, beta regression, parametric ordinal links, Weibull/frailty/copula models when using the full integrated likelihood |

This table is not permission to bulk-classify by filename. Each concrete class
must be checked for the objective it actually computes and the hooks it can
satisfy.

## Automated Invariants

Add source-level tests that enumerate every `Inference*` generator and enforce:

1. Every generator has metadata and a unique canonical name.
2. Every generator declares `abstract`, `exported`, and `likelihood_tier`.
3. Abstract generators cannot be instantiated and are absent from discovery.
4. Aliases identify one canonical generator and are absent from discovery.
5. Public optional methods exactly match declared capabilities.
6. No optional public method is implemented as an unconditional error.
7. No public/private method-name duplicate exists.
8. No class defines `supports_*`, `simulate_under_*`, or `get_*_spec` publicly.
9. Every component is registered and composed exactly once.
10. Every `private$...` reference in a component is classified as provided,
    required, optional, or owned.
11. Every required component hook exists on its final host.
12. No component redeclares root-owned state.
13. Every method collision has an explicit validated override declaration.
14. No source contains `eval(body(Inference...))`.
15. No source reads `$private` from an R6 generator.
16. No semantic classification depends on private method-name presence.
17. `InferenceSuite` returns exactly the compatible exported concrete classes.
18. Parametric-bootstrap classes pass a finite null-simulation smoke test.

Also retain behavioral regression tests for estimates, standard errors,
confidence intervals, p-values, randomization distributions, and bootstrap
distributions. Structural cleanup is not allowed to change numerical results
unless a test documents one of the confirmed bug fixes.

## Concrete TODO Backlog

This backlog is the implementable form of every structural point in this
report. Keep the TODOs in this section until each item is implemented,
tested, and either removed or marked complete in the commit that performs the
work.

### Structural Hierarchy TODOs

- [ ] Add a source-audit helper that enumerates every `Inference*` R6 generator,
  its parent, descendants, public methods, private methods, and aliases.
  Suggested target: `package_tests/audit_inference_hierarchy.R` or a helper
  under `EDI/tests/testthat/helper-inference-structure.R`.
- [ ] Add a generated manifest file containing canonical class, parent,
  aliases, abstract/exported status, response types, design compatibility,
  likelihood tier, direct components, effective components, and effective
  capabilities for every generator.
- [ ] Add a test that fails when a concrete class descends from algorithmic
  compatibility bases after that base's migration phase is complete.
- [ ] Add a test that no new concrete class inherits directly from
  `InferenceRand`, `InferenceRandCI`, `InferenceNonParamBootstrap`,
  `InferenceRandBootstrap`, `InferenceRandBootstrapCI`,
  `InferenceBayesianBootstrap`, `InferenceJackknife`, `InferenceAsymp`,
  `InferenceAsympLik`, or `InferenceParamBootstrap` unless the class is listed
  as an explicit migration bridge.
- [ ] For every current compatibility base, add a registry entry with
  `abstract = TRUE`, `exported = FALSE`, and a deletion phase.
- [ ] Add a phase counter test that reports how many concrete descendants
  remain under each compatibility base.

### Mixin Inventory TODOs

- [ ] Replace the seven-entry `EDI_MIXIN_CONTRACTS` shape with an enforced
  component schema while preserving backward-compatible field names during the
  transition.
- [ ] Mark `InferenceMixinCordeiroFerrariApprox` and
  `InferenceMixinLemonteGradientApprox` as `status = "scaffold"` or delete
  them. Add a test that scaffold components cannot appear in any effective
  component set.
- [ ] Register every real direct mixin host in a composition map, not only
  hosts that combine two or more mixins.
- [ ] Register every single-host `InferenceExt*` file-split bundle as either
  a temporary extension, an inline-owner candidate, or a real component.
- [ ] Add a source-level test banning raw
  `utils::modifyList(... InferenceMixin*$public/private ...)` composition
  outside `compose_inference_mixins()` or `define_inference_class()`.
- [ ] Keep `InferenceKKPassThroughCompound` and
  `InferenceKKPassThroughCompoundNoParamBootstrap` registered with the same
  explicit dependency on `InferenceMixinKKPassThrough`.

### Mixin Contract TODOs

- [ ] Extend each component entry with `provides_public_methods`,
  `provides_private_methods`, `provides_private_state`,
  `requires_public_methods`, `requires_private_methods`,
  `requires_private_state`, `optional_public_methods`, and
  `optional_private_methods`.
- [ ] Populate provided method/state fields from actual component list names
  and fail tests when declared `provides_*` differs from source.
- [ ] Populate required and optional references by parsing component bodies for
  `private$foo`, `self$foo`, and `super$foo`.
- [ ] Add `clear_nonestimable_state` to
  `InferenceMixinKKGLMMShared`'s required private methods.
- [ ] Add `y`, `w`, `X`, and `dead` to
  `InferenceMixinKKPassThrough`'s required private state.
- [ ] Add `supports_likelihood_tests` to
  `InferenceMixinKKPassThrough`'s required private methods until capability
  metadata replaces that hook.
- [ ] Classify `compute_fast_bootstrap_distr` and
  `compute_weighted_estimate_ivwc` as optional private hooks.
- [ ] For every class that composes a mixin, validate that the resolved parent
  chain plus host public/private lists satisfy all required component methods
  and state fields.
- [ ] Add a narrow allowlist for intentional `super$foo` calls in components;
  every allowlist entry must name the parent contract that supplies `foo`.
- [ ] Fail CI if a component gains a new `private$`, `self$`, or `super$`
  reference without a matching contract update.

### State Ownership TODOs

- [ ] Create a root-owned state list for `Inference` fields and a test that no
  component redeclares those fields in `owns_state`.
- [ ] Remove root-owned state redeclarations from
  `InferenceMixinKKPassThrough`, including `m`, `y_temp`, `dead`, `w`, `X`,
  `any_censoring`, `optimization_alg`, and `cached_mod`.
- [ ] Replace component state declarations with `requires_state` unless the
  component is the sole owner of that field.
- [ ] Add a state-collision test for public/private/state names across parent,
  components, and host lists.
- [ ] Add a migration allowlist for existing duplicate state fields, with a
  test that counts and forces the allowlist downward over phases.

### Function-Body Copy TODOs

- [ ] Add a source-level test that fails on `eval(body(Inference...))`.
- [ ] Replace copied calls to
  `InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T`
  with a named reusable helper or component-provided public method.
- [ ] Replace copied mixin bodies in KK OLS, KK mean-diff, KK Wilcox, KK
  quantile, KK survival, KK incidence, and KK count paths one family at a
  time.
- [ ] Add characterization tests before each replacement to prove bootstrap
  distributions and randomization paths remain numerically compatible.

### Super-Dependency TODOs

- [ ] Inventory every `super$foo` reference inside `InferenceMixin*` and
  `InferenceExt*` bodies.
- [ ] Move each `super$foo` dependency into either a named host requirement or
  an explicit parent-contract requirement.
- [ ] Refactor component code so it calls host/private service hooks where
  possible instead of depending on a particular superclass implementation.
- [ ] Add a test that a component cannot be installed on a parent that lacks a
  required `super$foo` contract.

### Capability Detection TODOs

- [ ] Replace `InferenceRandCI`'s `is_a_*` private-method sniffing with
  metadata queries.
- [ ] Remove statistical-family inference from `has_private_method()` calls;
  reserve that helper for optional implementation hooks only.
- [ ] Add metadata fields for response type, design family, estimator family,
  likelihood tier, and capabilities.
- [ ] Add tests that no source uses `has_private_method("is_a_` or
  `object_has_private_method(..., "is_a_` for semantic classification.
- [ ] Add tests that capability decisions call `obj$supports()` or metadata
  helpers rather than probing private method names.

### Alias And Discovery TODOs

- [ ] Add alias metadata for `InferenceAllKKCompoundMeanDiff` and
  `InferenceIncidExactZhang`.
- [ ] Add a deprecation schedule for each alias and document the canonical
  replacement class.
- [ ] Add a test that aliases resolve to the canonical generator but never
  appear as independent `InferenceSuite` candidates.
- [ ] Replace `InferenceSuite`'s manual denylist with metadata-based discovery
  once every generator has `abstract`, `exported`, and compatibility metadata.
- [ ] Add tests that unavailable optional packages are reported separately from
  design incompatibility.

### Class Metadata TODOs

- [ ] Create `inference_class_registry.R` with
  `EDI_INFERENCE_CLASS_REGISTRY`.
- [ ] Add `register_inference_class()`,
  `get_inference_class_metadata()`, `get_direct_components()`,
  `get_effective_components()`, `get_effective_capabilities()`, and
  `get_inference_alias_target()`.
- [ ] Add mandatory metadata fields: `name`, `abstract`, `exported`,
  `alias_of`, `response_types`, `design_families`, `compatibility`,
  `likelihood_tier`, `direct_components`, `required_packages`, and
  `migration_status`.
- [ ] Add metadata validation for allowed `likelihood_tier` values:
  `none`, `quasi`, `partial`, and `full`.
- [ ] Add tests that every metadata entry points to a real generator or real
  alias and every generator has exactly one canonical metadata entry.
- [ ] Add tests that metadata names match generator `classname` values.

### Component System TODOs

- [ ] Implement `InferenceComponent()` or an equivalent constructor for
  component records.
- [ ] Implement `resolve_inference_components()` with inherited components,
  deterministic dependency ordering, cycle detection, and duplicate detection.
- [ ] Implement `validate_component_contract()` for required hooks, required
  state, optional hooks, collisions, declared overrides, required packages,
  and likelihood-tier compatibility.
- [ ] Implement `assemble_public()` and `assemble_private()` so list splicing
  happens in one place.
- [ ] Implement `define_inference_class()` as a thin checked wrapper over
  `R6::R6Class()`.
- [ ] Add `migration_allow_legacy = TRUE` support only for temporary bridge
  classes; each use must include a deletion phase.
- [ ] Add an inspection helper that prints parent, direct components,
  effective components, capabilities, owned state, required hooks, optional
  hooks, and explicit overrides for a class.

### Component Inheritance TODOs

- [ ] Add tests proving `effective_components(child)` includes parent
  components without re-listing them.
- [ ] Add tests proving direct re-listing of an inherited component fails.
- [ ] Add tests proving component dependencies are inherited and ordered once.
- [ ] Add tests proving capability inheritance follows effective components.
- [ ] Add temporary suppression support only for migration bridges and count
  remaining suppressions in tests.

### Capability Rule TODOs

- [ ] Add a `capability_requires` table. It must encode rules such as:
  `parametric_likelihood_bootstrap` requires tier `partial` or `full`,
  `get_likelihood_test_spec`, and `simulate_under_lik_null`.
- [ ] Add a `public_methods_for_capability` table mapping each capability to
  public methods that should exist when the capability is present.
- [ ] Add tests that public optional method presence exactly equals the
  capability table for every concrete generator.
- [ ] Remove public optional throwing stubs after their capability components
  are extracted.
- [ ] Remove duplicate public/private `supports_*` hooks after metadata-backed
  `obj$supports()` is in place.
- [ ] Add tests that `likelihood_tier = "none"` never exposes LR, score,
  gradient, Bartlett, or parametric likelihood-bootstrap APIs.
- [ ] Add tests that `likelihood_tier = "quasi"` cannot expose likelihood-ratio
  APIs unless the method is renamed as an estimating-equation test and
  explicitly justified.
- [ ] Add tests that `partial` and `full` tiers do not imply parametric
  bootstrap without a null simulator.

### Component Extraction TODOs

- [ ] Extract `RandomizationTest` from `InferenceRand` as a component while
  leaving `InferenceRand` as an abstract bridge.
- [ ] Extract `RandomizationCI` from `InferenceRandCI` as a component.
- [ ] Extract `NonparametricBootstrap` from `InferenceNonParamBootstrap`.
- [ ] Extract `RandomizationBootstrap` from `InferenceRandBootstrap` and
  `InferenceRandBootstrapCI`.
- [ ] Extract `BayesianBootstrap` from `InferenceBayesianBootstrap`.
- [ ] Extract `Jackknife` from `InferenceJackknife`.
- [ ] Extract `Wald` and summary-table behavior from `InferenceAsymp` and
  `InferenceMLEorKMSummaryTable`.
- [ ] Extract `LikelihoodTests` from `InferenceAsympLik`.
- [ ] Extract `ParametricLikelihoodBootstrap` from `InferenceParamBootstrap`.
- [ ] Extract `StandardModelCache` from `InferenceAsympLikStdModCache`.
- [ ] Extract `CountLikelihoodPlumbing` from `InferenceCountLikelihood`.
- [ ] Extract KK pass-through, KK compound, KK GEE, and KK GLMM behavior into
  separately registered components with explicit dependencies.

### Likelihood Tier TODOs

- [ ] Add reviewed `likelihood_tier` metadata for every concrete class.
- [ ] Classify each concrete class from its implemented objective, not from
  filename or `OneLik` naming.
- [ ] Review all GEE, quasi-Poisson, modified-Poisson, robust, and sandwich
  classes as `quasi` candidates.
- [ ] Review Cox, stratified Cox, conditional logistic, conditional ordinal,
  and LWA Cox classes as `partial` candidates.
- [ ] Review GLM, count, beta, ordinal-link, Weibull, frailty, and copula
  classes as `full` candidates only when the implemented objective is a true
  normalized likelihood.
- [ ] Review exact, rank, Wilcoxon, quantile, matching/IVWC, g-computation,
  KM, log-rank, RMST, and Bai-style classes as `none` candidates.
- [ ] Add tests that tier changes require an explicit manifest diff and
  reviewer note.

### Suite Discovery TODOs

- [ ] Replace suite-level inferred compatibility metadata with registry-backed
  `compatibility` predicates once the class registry exists.
- [ ] Normalize design metadata in one helper used by every compatibility
  predicate.
- [ ] Add tests for continuous, incidence, count, proportion, survival, and
  ordinal designs.
- [ ] Add tests for KK and non-KK compatibility.
- [ ] Add tests for blocking and uncensored-only compatibility.
- [ ] Add tests that infrastructure, abstract, scaffold, and alias generators
  never appear in suite discovery.
- [ ] Add tests that constructor failures cannot influence discovery.

### Automated Invariant TODOs

- [ ] Create a dedicated structural test file, for example
  `EDI/tests/testthat/test-inference-hierarchy-invariants.R`.
- [ ] Implement each invariant listed in `Automated Invariants` as a separate
  `test_that()` block with failure messages naming offending classes.
- [ ] Add fixture helpers for parsing R6 generator definitions and component
  source bodies.
- [ ] Add a manifest-diff test that compares expected metadata to source
  metadata and fails on unreviewed changes.
- [ ] Add finite smoke tests for all classes that declare
  `parametric_likelihood_bootstrap`.
- [ ] Add a CI target that runs structural tests on every change under
  `EDI/R/inference*`.

### Behavioral Regression TODOs

- [ ] For each component extraction, add a before/after golden test for at
  least estimate, standard error, confidence interval, p-value, and bootstrap
  distribution where applicable.
- [ ] Reuse existing comprehensive tests as broad gates, but add focused tests
  for each migrated family so failures are diagnosable.
- [ ] Add targeted tests for Poisson, negative-binomial, OLS, Cox `OneLik`,
  copula `OneLik`, frailty `OneLik`, and conditional-logit `OneLik`
  parametric-bootstrap LR paths.
- [ ] Add targeted count-likelihood tests for Poisson, negative-binomial,
  hurdle, zero-inflated, and zero-augmented families.
- [ ] Add targeted standard-model-cache tests for the GLM, ordinal, survival,
  and proportion families listed in Phase 9.

### Cleanup TODOs

- [ ] Delete compatibility bridges after their descendant count reaches zero.
- [ ] Delete scaffold mixins that never acquire behavior.
- [ ] Delete obsolete `InferenceExt*` splices after they are converted to
  components or inlined into their single owner.
- [ ] Delete redundant `supports_*` flags after metadata-backed support is
  authoritative.
- [ ] Delete public optional throwing methods after public API presence is
  capability-driven.
- [ ] Delete alias exports only according to the documented deprecation
  schedule.
- [ ] Regenerate docs with `Rscript fast_roxygenize.R` after any exported API,
  class name, inheritance, or roxygen change.

## Migration Plan

### Migration constraints retained from the superseded plan

The following execution rules from the 2026-07-28 plan remain correct and
should apply to every phase:

1. Every phase must leave `devtools::load_all("EDI")` and the relevant
   targeted tests passing before the next phase starts.
2. Numerical outputs must stay unchanged unless a specific confirmed defect in
   this report is being fixed and covered by a regression test.
3. Run `Rscript fast_roxygenize.R` from the repository root after class
   renames, `inherit =` changes, exported method changes, or roxygen changes.
   Do not call `roxygen2::roxygenize()` directly.
4. Preserve `lock_objects = FALSE` on every R6 generator touched.
5. Keep the local naming convention: reusable components remain in
   `inference_mixin_<snake_case>.R` as `InferenceMixin<PascalCase>`;
   single-host file-split extensions, while they still exist, live in
   `inference_ext_<snake_case>.R` as `InferenceExt<PascalCase>`.
6. Use thin deprecated compatibility generators only as temporary migration
   bridges. They must preserve behavior for unmigrated concrete classes and
   must be deleted once no concrete generator descends from them.

### Phase 0: Freeze an authoritative inventory

Generate a machine-readable manifest for all current generators containing:

```text
canonical class, aliases, parent, abstract/exported status, response type,
design compatibility, likelihood tier, effective public methods,
effective private hooks, components/extensions, state fields, required packages
```

The existing hierarchy diagrams and reports are useful inputs, but source plus
the reviewed manifest becomes the migration truth. Add golden behavioral tests
before moving methods.

The manifest should explicitly identify every current list splice, including
the single-host "extension" objects named in the old plan:

- `InferenceMixinMinimumVolatilitySelector`
- `InferenceMixinExchangeableResamplingUnits`
- `InferenceMixinBcaBootstrapCI`
- `InferenceMixinPRWSubsampling`
- `InferenceMixinMOutOfNBootstrap`
- `InferenceMixinCIInversion`
- `InferenceMixinInformationMatrix`
- `InferenceMixinLikelihoodTestMemoization`
- `InferenceMixinParamBootstrapEstimate`
- `InferenceMixinBartlettApprox`
- `InferenceMixinCustomRandomizationStatistic`
- `InferenceMixinSequentialMCPval`
- `InferenceMixinQuantileRandCI`

If the manifest confirms that any of these are truly single-host
file-splitting devices, either rename them to `InferenceExt*` as a temporary
documentation cleanup or inline them into their owner during component
extraction. Do not treat the rename itself as architectural progress.

### Phase 1: Repair live faults and add guardrails

1. Fix the two count private parametric-bootstrap support overrides.
2. Remove the false public likelihood API from `InferenceAllKKMeanDiffIVWC`.
3. Remove the dead Newcombe generator-private expression.
4. Add the immediate missing abstract classes to the Suite denylist as a
   temporary containment measure.
5. Complete current mixin contracts and register all `InferenceExt*` splices.
6. Add collision, public/private-hook, and duplicate-component tests.
7. Add temporary characterization only as explicit metadata. Do not add new
   `private$is_a_*` probes for semantic discovery.

Do not wait for the full refactor to fix observable wrong behavior.

### Phase 2: Introduce checked composition

Implement `InferenceComponent()` and `define_inference_class()` without moving
behavior. Convert the current registry to the new schema. Make construction
fail on undeclared dependencies or collisions.

Provide an inspection helper that prints the resolved component order,
capabilities, owned state, required hooks, and explicit overrides for any
generator. This replaces hand inspection of nested list splices.

The first implementation may continue emitting Pattern-1 `list(public=,
private=)` bodies internally so behavior stays byte-identical. The important
change is that the lists are assembled from declared component contracts, with
duplicate components and accidental method collisions rejected up front.

### Phase 3: Remove the algorithm inheritance ladder

Extract one behavior group at a time in this order:

1. Randomization and randomization CI.
2. Nonparametric and randomization bootstrap.
3. Bayesian bootstrap.
4. Jackknife.
5. Wald/asymptotic inference.

For each extraction, leave a temporary deprecated compatibility generator for
unmigrated classes. New or migrated concrete classes must compose the component
directly. Delete each compatibility generator once it has no descendants.

The superseded plan's bridge strategy is still useful here. Convert large
algorithmic bases into components while leaving the old generator name as a
thin pass-through for existing descendants:

- `InferenceRand` and `InferenceRandCI` become randomization components.
- `InferenceNonParamBootstrap`, `InferenceRandBootstrap`, and
  `InferenceRandBootstrapCI` become bootstrap components.
- `InferenceBayesianBootstrap` becomes a Bayesian-bootstrap component.
- `InferenceJackknife` becomes a jackknife component.
- `InferenceAsymp` and `InferenceMLEorKMSummaryTable` become asymptotic
  estimate, standard-error, confidence-interval, and summary-table components.

Each extracted bridge must be mechanically equivalent first: copy current
public/private bodies verbatim, wire the old class name to the component, run
the focused family tests, then migrate concrete classes in small batches.

### Phase 4: Split likelihood semantics from likelihood APIs

1. Add required `likelihood_tier` metadata.
2. Extract likelihood score/gradient/LR behavior into
   `InferenceMixinLikelihoodTests`.
3. Extract null simulation and calibration into
   `InferenceMixinParametricLikelihoodBootstrap`.
4. Extract standard fitted-model cache behavior and count-specific likelihood
   plumbing into separately owned components.
5. Do not create an API-bearing common superclass for all likelihood tiers.

This phase supersedes the earlier proposal that all four tiers inherit
`InferenceAsympLik`.

Retain these extraction targets from the old plan, but attach them to the
component system instead of to likelihood-tier subclasses:

- Extract `InferenceParamBootstrap` into a parametric likelihood-bootstrap
  component. A temporary bridge may be named
  `InferenceMixinParamBootstrapSimulate` while migrating, but its final
  capability should be `parametric_likelihood_bootstrap`.
- Extract `inference_count_likelihood_public` and
  `inference_count_likelihood_private` from
  `InferenceCountLikelihood` into count-specific likelihood plumbing.
- Extract `inference_asymp_lik_std_mod_cache_public` and
  `inference_asymp_lik_std_mod_cache_private` from
  `InferenceAsympLikStdModCache` into standard fitted-model cache plumbing.
- Keep `InferenceParamBootstrap`, `InferenceCountLikelihood`, and
  `InferenceAsympLikStdModCache` only as deprecated compatibility generators
  while descendants are being migrated.

Do not preserve the old plan's instruction to make the compatibility
generators inherit a likelihood tier class. Their temporary inheritance may
remain whatever is required for behavioral parity until the concrete classes
compose the new components directly.

### Phase 5: Consolidate KK composition

Separate KK responsibilities into explicit components:

```text
matching/pass-through data access
compound or IVWC point estimation
bootstrap weight adaptation
GEE shared behavior
GLMM shared behavior
```

Resolve dependencies once, remove child re-splicing, eliminate
`eval(body(...))`, and move all shared state ownership to the root or one named
state component.

### Phase 6: Migrate concrete families

Migrate in reviewed batches:

1. Simple and no-likelihood estimators.
2. GEE, robust, sandwich, and other quasi-likelihood estimators.
3. Cox, conditional, and other partial-likelihood estimators.
4. Full-likelihood GLM, count, incidence, survival, and ordinal estimators.
5. KK and IVWC variants after their shared components are stable.

For each batch, compare the manifest before and after. A capability may change
only with a specific finding, implementation, and behavioral test.

Use `package_metadata/ivwc_capability_lattice.html` and
`package_metadata/non_ivwc_capability_lattice.html` as inputs to the reviewed
classification, not as unchecked truth. The source-level manifest and the
actual likelihood objective determine the final tier and capability set.

### Phase 7: Replace Suite discovery

Delete constructor probing and the manual base denylist. Select candidates from
metadata and use pure compatibility predicates. Report unavailable optional
packages separately from model incompatibility.

### Phase 8: Delete compatibility structure

Remove empty deprecated bases, redundant `supports_*` flags, throwing optional
methods, duplicate aliases, obsolete `InferenceExt*` splices, and unused mixin
scaffolds. Confirm that no concrete generator descends from an algorithmic
compatibility base.

### Phase 9: Enforce in CI

Run the structural invariant suite on every change to `EDI/R/inference*`. Run
focused family tests for changed components and the comprehensive package suite
before merging a migration phase.

The old plan's test gates should be kept as concrete checkpoints:

- Parametric bootstrap extraction must exercise Poisson, negative-binomial,
  OLS, Cox `OneLik`, copula `OneLik`, frailty `OneLik`, and conditional-logit
  `OneLik` paths that call
  `compute_lik_ratio_bootstrap_two_sided_pval()`.
- Count likelihood extraction must exercise Poisson, negative-binomial,
  hurdle, zero-inflated, and zero-augmented families, including
  `mark_count_likelihood_block_asymp_nonestimable()`.
- Standard-model-cache extraction must exercise logistic, probit, cauchit,
  cloglog, ordered-probit, proportional-odds, stereotype-logit,
  adjacent-category-logit, log-binomial, beta, zero-one-inflated-beta, Cox,
  stratified Cox, Weibull, dependent-censoring transform, incidence risk
  difference, and fractional-logit families.

## Definition of Done

The hierarchy is complete only when all of the following are true:

- Structural inheritance represents substitutable estimator types only.
- Likelihood tier is explicit and valid for every concrete class.
- Optional public method presence exactly equals declared capability.
- No class relies on private support flags to disable inherited public APIs.
- No public and private capability hooks share a name.
- No concrete class inherits an unsupported inference algorithm.
- No component or transitive dependency is composed twice.
- No method or state collision is implicit.
- Every component dependency is declared and validated.
- Every mutable field has one owner.
- No mixin uses `super`, implementation-body evaluation, or generator-private
  access.
- `InferenceSuite` lists only compatible, exported, concrete estimators.
- Aliases never appear as separate estimator types.
- All confirmed defects in this report have regression tests.
- Numerical golden tests pass across every migrated family.
- The old algorithmic compatibility bases have no descendants and are removed.

Until these conditions hold, the hierarchy should be treated as a migration in
progress rather than as a reliable capability model.
