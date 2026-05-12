# Nominal Response Type Report

## Scope

This report evaluates how difficult it would be to add a new
`response_type = "nominal"` to the package, where "nominal" means:

- a categorical outcome
- with no ordering of the levels
- represented most naturally as an unordered factor

This is different from the existing `ordinal` response type, which already has:

- integer-coded outcome storage
- ordered-factor support
- level bookkeeping in the `Design` base class
- several ordinal regression inference paths

The report covers:

- the core package plumbing needed to admit a new response type
- how the existing `inference_all_*` abstractions would respond
- a pragmatic set of nominal-specific inference paths to implement first

## Short Answer

Adding the **label** `nominal` is easy.

Adding a **coherent nominal response type** that works across designs,
`SimulationFramework`, `InferenceSuite`, and a useful set of inference methods
is a **moderate-to-hard** project.

The main reason is that the current package architecture is strongly organized
around response families with clear numeric semantics:

- continuous
- incidence
- proportion
- count
- survival
- ordinal

Nominal breaks that pattern in two ways:

1. it is not naturally scalar on an ordered numeric scale
2. many generic "all-subject" inference paths currently assume a scalar
   treatment effect like a mean difference, rank shift, or single regression
   coefficient

So the real work is not adding the enum. It is deciding what the nominal
estimand is and making the generic inference infrastructure respect that choice.

## How Common Are Nominal Outcomes In Experimental Literatures?

The short answer is:

- **very common in principle**
- **less dominant than binary and continuous outcomes in many trial literatures**
- **common enough that a serious experimental-design package should have a plan
  for them**

### Medical science and clinical trials

In clinical-trial practice, binary outcomes are probably the most common
categorical endpoint family. A 2020 review of 200 randomized trial reports in
medicine notes that binary outcomes are widely used and cites earlier work
finding that about half of trials calculated sample size from a binary outcome;
see Rombach et al. (2020), BMC Medicine:
https://link.springer.com/article/10.1186/s12916-020-01598-7

That means nominal outcomes are **not** the modal endpoint family in the core
RCT literature. However, they are clearly present. Two useful signals are:

- the Trials protocol by Selman et al. explicitly treats nominal outcomes as a
  standard trial outcome scale distinct from ordinal and continuous outcomes;
  see:
  https://link.springer.com/article/10.1186/s13063-023-07262-8
- Schmid, Trikalinos, and Olkin (2014) develop Bayesian network meta-analysis
  specifically for unordered categorical outcomes and apply it to 17 trials with
  mutually exclusive nominal outcomes such as cardiovascular death,
  non-cardiovascular death, and no death; see PubMed:
  https://pubmed.ncbi.nlm.nih.gov/26052655/

So the right characterization for medicine is:

- binary and ordinal endpoints are more common than nominal endpoints as primary
  trial outcomes
- but nominal outcomes are established enough that there is active methodological
  work for trials and evidence synthesis

### Economics field experiments

Field experimentation in economics has expanded substantially in scope. Levitt
and List's overview emphasizes that the range of field experiments and the types
of questions studied have grown tremendously; see NBER Working Paper 14356:
https://www.nber.org/papers/w14356

That matters for nominal outcomes because many economic experiments naturally
generate unordered categorical endpoints such as:

- choice among alternatives
- program participation status categories
- employment-status categories
- market or political choice categories

The exact prevalence is harder to summarize with one headline percentage than in
clinical-trial reviews, but nominal outcomes are conceptually very natural in
this literature.

### Social and political randomized experiments

Nominal outcomes are especially natural in social and political experiments,
because many substantive endpoints are genuinely category choices rather than
ordered scales.

Examples:

- Coppock, Green, and Porter study **vote choice / vote share** in a randomized
  field experiment on digital political advertising; see:
  https://journals.sagepub.com/doi/10.1177/20531680221076901
- Gerber, Huber, and Washington study **party affiliation** and related
  political beliefs in a field experiment; see NBER Working Paper 15365:
  https://www.nber.org/papers/w15365

In these settings, unordered categories such as:

- party identification
- vote choice
- policy preference category
- employment or schooling status category

are not edge cases. They are central substantive outcomes.

### Practical interpretation for this package

So the package should not treat nominal outcomes as exotic. The better summary
is:

- nominal outcomes are **less central than binary outcomes** in mainstream
  clinical-trial methodology
- nominal outcomes are **clearly present** in medical and biostatistical
  methodology
- nominal outcomes are **very natural** in economics, political science, and
  social experiments whenever the outcome is a choice among more than two
  unordered categories

That strengthens the case for adding `response_type = "nominal"` eventually,
but it also argues for implementing it in a way that matches the actual applied
use cases:

- global category-distribution tests
- focal-category contrasts
- multinomial choice / multinomial logit style models

## What Already Exists

The package already has one categorical response family beyond binary:

- `ordinal`

That matters because the `Design` base class already has special response-type
state for:

- `ordinal_levels`
- `original_ordinal_levels`
- coercion of ordered factors to integer storage

See [design_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:42),
[design_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:63),
and [design_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:119).

So nominal can reuse the idea that a categorical response may carry labels that
are distinct from the internal numeric storage. But unordered categories still
need their own semantics.

## Difficulty By Layer

### 1. Design base classes: moderate

The `Design` base class currently validates:

- `continuous`, `incidence`, `proportion`, `count`, `survival`, `ordinal`

in [design_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:56).

To add `nominal`, the `Design` layer would need:

- `nominal` added to the allowed `response_type` set
- support for unordered factors in `add_one_subject_response()`
- support for unordered factors in `add_all_subject_responses()`
- a `nominal_levels` / `original_nominal_levels` analogue to the ordinal fields
- an update to `assert_y()` so it validates nominal category encoding
- an update to `transform_y()` if transformed nominal outcomes are ever allowed

This is not difficult, but it is real work because the current code treats
non-ordinal responses as numeric and `ordinal` as the only factor-backed
special case.

### 2. Most design classes: easy

Most design classes are response-type agnostic after initialization. They carry:

- assignments
- covariates
- outcomes

and usually do not model the response directly.

So once `Design` itself understands `nominal`, many concrete design classes
would work unchanged or nearly unchanged.

### 3. KK21 / response-adaptive weighting designs: moderate to hard

This is where nominal is not free.

`DesignSeqOneByOneKK21` and related paths explicitly branch on `response_type`
to compute covariate weights from the observed response family; see
[design_seq_one_by_one_KK21.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_seq_one_by_one_KK21.R:253)
through [design_seq_one_by_one_KK21.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_seq_one_by_one_KK21.R:281).

Today they have hand-written logic for:

- continuous
- incidence
- count
- proportion
- survival
- ordinal

Nominal would need:

- either a genuine nominal-response weighting model
- or a documented approximation / speedup path

So response-adaptive KK weighting is one of the harder design-side additions.

### 4. `SimulationFramework`: hard

`SimulationFramework` is a major separate work item.

It currently:

- validates the response-type enum in
  [simulations_framework.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:499)
- transforms latent continuous signals to response-family scale in
  [simulations_framework.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:126)
- defines response-type-specific treatment-effect semantics in
  [simulations_framework.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:2133)
  and
  [simulations_framework.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:2887)
- chooses curated default inference classes by response type in
  [simulations_framework.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/simulations_framework.R:3106)

For nominal, the framework would need new answers to all of these questions:

- how do we generate a latent nominal outcome from `y_cont`?
- what is the treatment effect `betaT` for a nominal response?
- what is the scalar "true effect" used for MSE, coverage, and power summaries?
- which default inference classes are the curated nominal set?

This is one of the biggest architectural costs.

## The Main Conceptual Problem: What Is The Nominal Estimand?

For current response types, the package usually has a natural scalar estimand:

- mean difference
- log-odds ratio
- log-risk ratio
- log-rate ratio
- log-hazard / log-time ratio
- ordinal-location-type coefficient

For nominal, there is no single obvious scalar effect.

Common possibilities are:

- a vector of category-specific log-odds ratios relative to a baseline category
- a contrast for one chosen category vs all others
- a treatment effect on expected utility under user-specified category scores
- a global null statistic: no treatment effect on the entire category
  distribution

This choice drives everything else:

- model output shape
- CI/p-value APIs
- `SimulationFramework` truth
- summary-table shape
- which generic inference classes can make sense

Without fixing this, `nominal` remains a storage type but not an inferential
family.

## How The Existing `inference_all_*` Paths Would Respond

### `InferenceAllSimpleMeanDiff`: should reject

This class currently computes:

- mean of treated outcomes minus mean of control outcomes

with no response-type assertion in the initializer; see
[inference_all_mean_diff.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_mean_diff.R:29).

That means if `nominal` were added naively and internally encoded as integers,
this class would likely run and produce a meaningless result based on arbitrary
level codes.

So for nominal, this class should be updated to:

- explicitly reject `response_type = "nominal"`

unless the package deliberately chooses a scored nominal utility framework,
which would be a different feature.

### `InferenceAllSimpleWilcox`: should reject

This class currently supports:

- `continuous`, `count`, `proportion`, `survival`, `ordinal`

and explicitly rejects `incidence`; see
[inference_all_simple_wilcox.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_simple_wilcox.R:74).

Nominal should also be explicitly rejected because:

- rank-based location-shift logic does not apply to unordered categories
- any integer encoding of category labels would be arbitrary

So this is another path that needs an explicit guard, not accidental acceptance.

### `InferenceMLEorKMSummaryTable`: neutral abstract layer

This is an abstract summary-table layer; see
[inference_all_abstract_mle_or_KM_summary_table.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_mle_or_KM_summary_table.R:8).

By itself, it is not nominal-specific and does not block nominal. But its
assumption is still:

- a coefficient vector
- a distinguished treatment coefficient
- a scalar standard error for that coefficient

So it can support nominal only if the concrete nominal inference class chooses a
scalar treatment estimand or a designated baseline-category treatment
coefficient.

This layer is reusable, but only for certain nominal-model designs.

### `InferenceAsympLikStdModCache`: reusable for scalar nominal models

This abstract class is the common path for many likelihood-based families; see
[inference_all_abstract_asymp_lik_std_mod_cache.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R:6).

It assumes:

- `generate_mod()` returns an object with `b[2]` as the treatment effect
- the treatment effect is scalar

So this is reusable for a nominal model only if we define the nominal treatment
effect as:

- one coefficient, usually for treatment in a chosen category contrast

It is **not** a natural fit for a fully vector-valued multinomial treatment
effect without further extension.

### Randomization / bootstrap base classes: mostly neutral

The generic randomization and bootstrap base classes are not inherently tied to
likelihoods. In principle they can support nominal outcomes if a concrete class
provides:

- a scalar statistic
- re-estimation under permutations or resamples

So from a software-architecture perspective, these are not the main blocker.
The blocker is still the choice of nominal statistic / estimand.

### `InferenceSuite`: will adapt reasonably once errors are explicit

`InferenceSuite` discovers compatibility by trying to instantiate each class and
looking for incompatibility messages; see
[inference_suite.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_suite.R:110).

That means `InferenceSuite` will work fine for nominal **if and only if**:

- non-nominal classes reject nominal clearly
- nominal classes accept it clearly

If generic classes like mean-diff are not fenced off, `InferenceSuite` may
incorrectly list them as applicable. So explicit rejection paths are important.

## Package Areas That Need Explicit `nominal` Fences

Even before implementing any nominal-specific inference, the following classes
should likely be updated to reject nominal explicitly:

- `InferenceAllSimpleMeanDiff`
- `InferenceAllSimpleMeanDiffPooledVar`
- `InferenceAllSimpleWilcox`
- any quantile-regression-based abstractions
- any regression path that assumes a scalar ordered/numeric outcome but does not
  already `assertResponseType(...)`

The package is currently safest when a concrete class says:

- "I work for `ordinal`"
- "I work for `count`"
- etc.

Nominal will expose any generic path that currently relies on numeric outcome
storage without sufficiently specific response-type gating.

## Common Nominal-Specific Inference Paths To Implement

The best first nominal paths are the ones that give:

- a standard applied interpretation
- one or more scalar treatment estimands
- a reasonable fit with the package's current asymptotic / bootstrap APIs

### 1. Multinomial logistic regression (baseline-category logit)

This should be the main nominal likelihood path.

Typical setup:

- choose one response category as the baseline
- fit a multinomial logit model with treatment and covariates
- define treatment effects as the category-specific treatment log-odds ratios

Difficulty: **moderate to hard**

Why:

- statistically standard and familiar
- nominal-exclusive in the right way
- but naturally vector-valued, which conflicts with the package's mostly scalar
  treatment-effect API

A pragmatic first version could reduce this to one scalar estimand by:

- targeting a user-chosen focal category vs baseline

That would make it fit the current asymptotic interface much more easily.

### 2. Category-specific binary contrast inference

This is the easiest nominal path conceptually.

Idea:

- pick one focal category `k`
- define a derived binary outcome `1{Y = k}`
- run existing incidence-style inference on that contrast

Difficulty: **easy to moderate**

Why:

- this reuses existing incidence machinery
- it gives an immediately interpretable estimand
- it avoids building a full multinomial stack on day one

Limit:

- it is not a full nominal model
- it needs multiple-testing or family-wise guidance if used across many
  categories

Still, this is an attractive first nominal path because it is useful and fits
the current software design.

### 3. Global two-sample nominal-distribution test

Examples in the applied world:

- Pearson chi-squared test for treatment vs response table
- likelihood-ratio `G^2` test
- Fisher-Freeman-Halton exact test for small samples

Difficulty: **moderate**

Why:

- this gives a clean global null: treatment does not change the category
  distribution
- it does not require choosing a baseline category
- it is naturally nominal

Limitation:

- this is mainly a test, not a scalar treatment-effect estimate with a CI

So this is a strong candidate for:

- an exact / asymptotic / randomization p-value path

but not necessarily for the package's full estimate-plus-CI interface.

### 4. Stratified nominal test / CMH-style extension

For blocked or stratified designs, a nominal analogue of CMH is relevant.

Difficulty: **moderate to hard**

Why:

- practically useful for blocked designs
- but there are multiple CMH-type generalizations for nominal data
- and the package would need to decide whether the estimand is:
  - a common association test,
  - a category-specific contrast,
  - or a model-based coefficient

This is a second-wave nominal path, not the first one.

## Recommended First-Wave Nominal Inference Set

If the goal is to make `nominal` useful without exploding scope, the best first
wave is:

1. one global nominal test
2. one focal-category incidence-style contrast path
3. one multinomial-regression path if a scalar focal-category coefficient is
   acceptable

Concretely:

- `InferenceNominalPearsonChisq`
- `InferenceNominalCategoryContrast`
- `InferenceNominalMultinomLogit`

These three cover:

- global testing
- simple effect estimation
- covariate-adjusted model-based inference

without forcing the package to solve a full vector-valued nominal-effect API on
day one.

## `SimulationFramework` Implications

Nominal is especially expensive in `SimulationFramework`.

### Data generation

The current generator transforms a latent continuous signal into:

- Bernoulli probability
- proportion mean
- count rate
- survival scale
- ordinal category by cutpoints

Nominal would need something else, for example:

- latent utilities for each category with softmax probabilities
- or a reference-category logit generative model

This is not a one-line extension of the existing `switch(...)`.

### Treatment effect semantics

`betaT` is currently scalar and response-family specific.

For nominal, possibilities include:

- shifting the log-odds of one focal category vs baseline
- shifting a whole vector of category logits
- shifting latent utilities for all categories

This must be specified before simulation output is meaningful.

### Truth and summary metrics

The current framework computes scalar truth for:

- MSE
- coverage
- power

A fully multinomial nominal effect would be vector-valued, so either:

- `SimulationFramework` must be generalized to vector estimands
- or nominal simulation must target a scalar estimand such as one focal
  category's log-odds contrast

This is why nominal support in `SimulationFramework` is one of the hardest
parts.

## Recommended Implementation Plan

### Stage 1: Admit the type safely

Add `nominal` to the core response-type plumbing:

- `Design` validation
- response storage / factor coercion
- documentation

At the same time, add explicit rejection guards to generic inference classes
that would otherwise produce nonsense from integer-coded nominal levels.

### Stage 2: Add one nominal global test

Implement an asymptotic nominal distribution-comparison path such as:

- Pearson chi-squared
- or likelihood-ratio `G^2`

This gives immediate nominal utility with minimal ambiguity.

### Stage 3: Add focal-category contrast inference

Implement a scalar nominal effect by selecting a category and modeling:

- `Y = category_k` vs `Y != category_k`

This reuses a large amount of existing incidence machinery and fits the
package's scalar-estimand architecture.

### Stage 4: Add multinomial logistic regression

Implement a nominal-specific regression path, ideally with:

- focal-category scalar treatment effect first
- full vector-valued treatment effects later if desired

This is the right place to decide whether the package wants to support
vector-valued inferential outputs more generally.

### Stage 5: Extend `SimulationFramework`

Only after the estimand is settled should `SimulationFramework` be updated for:

- nominal data generation
- scalar truth definition
- curated default nominal inference set

## Bottom Line

Adding `response_type = "nominal"` is not hard at the enum / storage level.

What is hard is making it statistically coherent across the package.

The main engineering and design conclusions are:

- the `Design` base class can be extended with moderate effort
- most design classes will inherit nominal support cheaply
- response-adaptive KK weighting designs are a harder special case
- many generic `inference_all_*` paths should explicitly reject nominal rather
  than accidentally accepting integer codes
- the best first nominal inference paths are:
  - a global nominal test,
  - a focal-category contrast path,
  - and later a multinomial logistic regression path
- `SimulationFramework` is the hardest integration point because nominal does
  not come with a natural scalar `betaT` semantics

So the pragmatic recommendation is:

1. add `nominal` safely at the design/container level
2. fence off generic numeric-only inference paths
3. implement one or two nominal-specific scalar/test paths first
4. postpone full multinomial and `SimulationFramework` support until the
   nominal estimand is explicitly chosen
