# Fix Documentation TODOs

Generated from the current `EDI/man/*.Rd` public API surface. This file is a work plan, not user documentation.

- Public Rd topics found: `252`
- Public R6 methods found: `627`
- Total API entries listed: `879`

## General Instructions

These instructions apply to every numbered TODO in the thin-description list below.

### Documentation Standard

For every public API entry below, replace short or generic prose with documentation that answers:

- What estimand, parameter, contrast, or design quantity is being computed.
- What model, estimating equation, randomization distribution, bootstrap law, or likelihood is used.
- The mathematical form of the key statistic, loss, likelihood, score, variance estimator, or CI inversion when relevant.
- What assumptions are required for validity, including design, response type, censoring, matching, blocking, clustering, or asymptotics.
- What input conventions the method assumes. Document required dimensions, types, ordering, treatment coding, factor levels, matching or block identifiers, weights, offsets, censoring indicators, tie handling, missing-value behavior, zero-count behavior, and boundary cases.
- How function arguments map to mathematical symbols and model parameters. State distribution parameterizations explicitly, including link functions, baseline or reference categories, nuisance parameters, treatment-coefficient indices, and whether parameters are on the natural, log, logit, hazard, odds, probability, time, or transformed scale.
- How hypothesis tests and intervals are defined. State the null and alternative hypotheses, one-sided or two-sided tail conventions, confidence level convention, whether intervals are obtained by Wald formulas, test inversion, profile likelihood, bootstrap quantiles, or randomization inversion, and how exact, asymptotic, bootstrap, jackknife, Bayesian-bootstrap, and randomization paths differ.
- What happens when the quantity is not estimable or when the class does not support the requested inference path.
- Whether the method mutates R6 object state. Document cached fields, call-order requirements, stale-cache invalidation, cloning behavior, reproducibility state, and whether repeated calls are deterministic for fixed inputs and seeds.
- Numerical implementation details that affect reproducibility, accuracy, or interpretation. Document optimizers, starting values, constraints, tolerances, gradients, Hessians, convergence criteria, line-search or trust-region behavior, stabilization tricks, fallbacks, and warnings.
- Randomness and reproducibility behavior. Document which RNG source is used, how seeds are consumed, whether parallel execution changes draw order, how Monte Carlo error is reported or should be assessed, and whether randomization or bootstrap draws can be reused.
- Computational complexity and practical limits when relevant. State the dominant scaling in subjects, pairs, blocks, clusters, covariates, bootstrap replicates, randomization draws, or outcome levels, and identify known memory or runtime bottlenecks.
- Contracts for `fast_*` and C++ backend functions. Document indexing conventions, memory layout, unchecked assumptions, domain restrictions, overflow and underflow behavior, NA/NaN handling, convergence flags, and the relationship between low-level backends and safer R wrappers.
- Links to related parent methods instead of repeating shared bootstrap, randomization, jackknife, exact, or asymptotic contracts.
- Do not repeat yourself: place shared documentation at the highest possible level in the inheritance chain, then have subclasses link to that parent documentation and document only subclass-specific behavior.
- References with stable identifiers, preferably DOI, arXiv, or package/manual citations.
- References and links for specialized numerical or statistical ingredients, including approximations, optimizers, quadrature rules, transforms, corrections, and named algorithms. For example, document and cite Lanczos approximation, Stirling approximation, saddlepoint approximations, Gauss-Hermite quadrature, Cholesky decompositions, log-sum-exp stabilization, Bartlett corrections, and other esoteric implementation details when they affect formulas, accuracy, or interpretation.
- HTML links to analogous Python package documentation when it helps users compare APIs or verify model conventions.
- HTML links to Wikipedia pages only as secondary orientation aids; do not use Wikipedia as the primary source for statistical validity, formulas, or implementation details.
- Interpretation guidance and misuse warnings where needed. State causal, design-based, model-based, exchangeability, positivity, independent-censoring, proportional-hazards, proportional-odds, or large-sample assumptions clearly enough that users can tell when a result should not be used.
- Lifecycle status for APIs that are experimental, internal-but-exported, superseded, deprecated, or thin wrappers around lower-level routines.

### First-Rate Documentation Workstreams

These are cross-cutting documentation improvements that should be completed alongside the numbered API TODOs.

- Create theory vignettes by method family: `theory-incidence.Rmd`, `theory-count.Rmd`, `theory-survival.Rmd`, `theory-ordinal.Rmd`, `theory-bootstrap-randomization.Rmd`, and any additional family needed for continuous/proportion methods. Individual roxygen entries should link into these vignettes rather than repeating long derivations.
- Create a notation glossary defining package-wide symbols and conventions, including `W`, `Y`, `X`, `delta`, `beta_T`, matched pairs, reservoir subjects, blocks, clusters, censoring indicators, bootstrap weights, Bayesian-bootstrap weights, randomization permutations, treatment-effect scales, and transformed outcomes.
- Create support tables for all inference classes. Each row should include response type, allowed design types, estimand, effect scale, exact/asymptotic/bootstrap/Bayesian-bootstrap/jackknife/randomization support, likelihood-test support, required packages, and known fallback or unsupported behavior.
- Document return object schemas for every public method or `fast_*` backend that returns a list or structured object. Include field names, dimensions, parameter order, treatment-coefficient index, convergence fields, variance fields, and when fields may be `NA`, `NULL`, or omitted.
- Add explicit "scale of effect" notes throughout the inference docs. State whether the returned estimate is a mean difference, risk difference, log risk ratio, odds ratio, log odds ratio, log hazard ratio, log-time ratio, RMST difference, quantile shift, rank statistic, or other scale.
- Create a centralized validation and failure-semantics page. Cover `is_nonestimable()`, `get_nonestimable_reason()`, `get_nonestimable_stage()`, hardening, convergence failure, model separation, insufficient data, non-finite bootstrap draws, missing standard errors, unsupported exact/asymptotic/bootstrap/randomization methods, and how these states appear in return values.
- Create a method-selection guide. Include decision tables that help users choose an inference class from the design, outcome type, estimand, treatment-effect scale, sample size, censoring structure, clustering/blocking structure, and desired inference mode.
- Create an input-conventions page. Define package-wide conventions for `W`, `Y`, covariate matrices, model matrices, treatment coding, factor ordering, block and cluster identifiers, matched-pair order, weights, offsets, missingness, censoring indicators, ties, and boundary values.
- Create an R6 behavior page. Document which public methods are pure computations, which mutate object state, which cache results, which require a prior call to another method, how cloning behaves, and how users should reason about repeated calls.
- Create a reproducibility page. Document RNG usage for simulation, bootstrap, Bayesian-bootstrap, and randomization methods; explain seed handling, parallel execution, draw reuse, Monte Carlo error, and how to reproduce published examples.
- Create a backend-contracts page for `fast_*` and C++ utilities. Document low-level assumptions that are intentionally not checked, argument dimensions and storage order, numeric domains, overflow/underflow safeguards, convergence flags, and wrapper-to-backend equivalence.
- Add cross-language analog links for each model family. Use statsmodels, lifelines, scipy, scikit-survival, or other official package documentation where useful, and label these links as analogous APIs or implementation references rather than EDI dependencies.
- Add stable citation keys and a reference map. Prefer `\doi{...}`, arXiv links, package manuals, or canonical book references. Maintain a central `REFERENCES.md` or package citation map from method family/class to references.
- Add validation-evidence references. For important methods, link examples, tests, or vignettes to published numerical examples, simulation checks, package-to-package comparisons, or closed-form special cases that demonstrate the implementation matches the documented formulas.
- Improve examples for public classes. Each exported class should have one tiny runnable example, one realistic `\donttest{}` example where appropriate, and one example or note showing the returned estimand scale.
- Improve pkgdown organization. Group APIs by design family, outcome family, inference mode, and backend utility; add search-friendly aliases for common statistical terms; ensure equations render correctly; and verify that reference pages link to the relevant theory vignettes, support tables, and external references.
- Add documentation tests. The test should parse generated Rd files and fail on placeholder descriptions such as "TODO", "Compute pval", "Initialize object", missing references for theoretical methods, broken links, repeated boilerplate that should live in a parent, public methods without argument/return documentation, undocumented defaults, unresolved mathematical symbols, malformed equations, misspellings, and exported objects missing from the reference index.

### External HTML Link Policy

Use three classes of links in the roxygen pass:

- Primary statistical references: papers, arXiv pages, DOI landing pages, books, package vignettes, or package manuals that define the method.
- Analogous software documentation: official Python/R/SAS documentation showing a comparable model API, parameterization, likelihood, optimizer, or returned object.
- Orientation links: Wikipedia pages for common concepts such as likelihood functions, score tests, Cox models, bootstrap, or special functions. These are useful for readers but should not replace primary references.

When adding analogous Python links, explicitly say they are analogs rather than dependencies. For example, a `fast_poisson_regression_cpp()` doc can link to statsmodels' Poisson/GLM documentation as a comparable Python API while still documenting EDI's own parameterization and C++ return fields.

### Analogous Python Documentation Links

Use these links as a starter map while fixing TODOs:

- Generalized linear models and exponential-family conventions: [statsmodels GLM](https://www.statsmodels.org/stable/glm.html).
- Discrete, binary, count, hurdle, zero-inflated, conditional-logit, and conditional-Poisson models: [statsmodels discrete models](https://www.statsmodels.org/stable/discretemod.html).
- Conditional Poisson likelihood with grouped intercepts conditioned out: [statsmodels ConditionalPoisson](https://www.statsmodels.org/dev/generated/statsmodels.discrete.conditional_models.ConditionalPoisson.html).
- Cox proportional-hazards models and survival APIs: [lifelines CoxPHFitter](https://lifelines.readthedocs.io/en/latest/fitters/regression/CoxPHFitter.html), [statsmodels duration models](https://www.statsmodels.org/dev/duration.html), and [scikit-survival documentation](https://scikit-survival.readthedocs.io/).
- Survival weighting and robust/sandwich examples: [lifelines examples](https://lifelines.readthedocs.io/en/latest/Examples.html).
- Special functions used by `fast_*` math utilities: [SciPy special functions](https://scipy.github.io/devdocs/reference/special.html) and [SciPy digamma](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.digamma.html).
- Distribution functions and density parameterizations used by fast likelihood code: [SciPy stats beta](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.beta.html), [SciPy stats distributions index](https://docs.scipy.org/doc/scipy/reference/stats.html).

### Wikipedia Orientation Links

Use these sparingly as "See also" orientation links, paired with primary references:

- [Likelihood function](https://en.wikipedia.org/wiki/Likelihood_function)
- [Maximum likelihood estimation](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation)
- [Score test](https://en.wikipedia.org/wiki/Score_test)
- [Likelihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test)
- [Wald test](https://en.wikipedia.org/wiki/Wald_test)
- [Bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
- [Randomization test](https://en.wikipedia.org/wiki/Randomization_test)
- [Generalized linear model](https://en.wikipedia.org/wiki/Generalized_linear_model)
- [Logistic regression](https://en.wikipedia.org/wiki/Logistic_regression)
- [Poisson regression](https://en.wikipedia.org/wiki/Poisson_regression)
- [Negative binomial distribution](https://en.wikipedia.org/wiki/Negative_binomial_distribution)
- [Beta distribution](https://en.wikipedia.org/wiki/Beta_distribution)
- [Ordinal regression](https://en.wikipedia.org/wiki/Ordinal_regression)
- [Quantile regression](https://en.wikipedia.org/wiki/Quantile_regression)
- [Cox proportional hazards model](https://en.wikipedia.org/wiki/Proportional_hazards_model)
- [Kaplan-Meier estimator](https://en.wikipedia.org/wiki/Kaplan%E2%80%93Meier_estimator)
- [Log-rank test](https://en.wikipedia.org/wiki/Logrank_test)
- [Gamma function](https://en.wikipedia.org/wiki/Gamma_function)
- [Digamma function](https://en.wikipedia.org/wiki/Digamma_function)
- [Lanczos approximation](https://en.wikipedia.org/wiki/Lanczos_approximation)
- [Stirling's approximation](https://en.wikipedia.org/wiki/Stirling%27s_approximation)
- [Gauss-Hermite quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature)
- [LogSumExp](https://en.wikipedia.org/wiki/LogSumExp)

### Reference Backlog

Add exact citations while editing the method docs. At minimum, verify and cite the following where relevant:

- Zhang exact/randomization incidence inference paper on arXiv.
- Azriel et al. paper on CMH / covariate-adaptive incidence inference on arXiv.
- Bai adjusted t-test / matched-pair plus reservoir estimator references.
- Lin (2013) covariate-adjusted estimator reference.
- Kaplan-Meier, log-rank, Gehan-Wilcoxon, RMST, Cox PH, stratified Cox, Weibull AFT, frailty, and Clayton-copula survival references.
- GLM and likelihood references for Wald, score, likelihood-ratio, gradient, Bartlett correction, quasi-likelihood, GEE, GLMM, CLMM, conditional logit, hurdle, zero-inflated, negative-binomial, beta, fractional-logit, and quantile-regression models.
- Bootstrap references for percentile/basic, bootstrap-t/studentized, BCa, Bayesian bootstrap, m-out-of-n bootstrap, PRW subsampling, double bootstrap, prepivoting, and calibrated/smoothed variants.
- Randomization-inference references for sharp-null permutation tests, custom statistics, and confidence-interval inversion.

### Primary Reference Hierarchy

Use the strongest primary reference available for each method:

1. Method-defining papers for package-specific or named procedures. Use these for Zhang exact/randomization incidence methods, Azriel CMH/covariate-adaptive incidence methods, Bai adjusted-t methods, Lin covariate adjustment, Newcombe intervals, Miettinen-Nurminen intervals, Bartlett corrections, BCa bootstrap, PRW subsampling, and any KK matching-on-the-fly estimators.
2. Journal papers or arXiv preprints for recent methods that do not yet have canonical textbook treatment. Prefer the published journal version when available; cite arXiv when it is the only stable public version or when the implementation follows the preprint exactly.
3. Textbooks and monographs for standard model families and asymptotic theory. Use these for GLMs, likelihood theory, Wald/score/LR tests, survival analysis, categorical/ordinal models, count regression, randomization inference, bootstrap theory, and estimating equations.
4. Software manuals only for implementation analogs, API comparisons, parameterization checks, and returned-object conventions. Do not cite software manuals as the sole justification for the statistical method unless the documented object is itself a software algorithm.
5. Wikipedia only as a "See also" orientation link, never as the primary statistical reference.

Good textbook or monograph candidates to map into TODOs:

- GLM and likelihood theory: McCullagh and Nelder, *Generalized Linear Models*; Casella and Berger, *Statistical Inference*; Pawitan, *In All Likelihood*.
- Count regression: Cameron and Trivedi, *Regression Analysis of Count Data*.
- Categorical and ordinal data: Agresti, *Categorical Data Analysis* and *Analysis of Ordinal Categorical Data*.
- Survival analysis: Cox and Oakes, *Analysis of Survival Data*; Kalbfleisch and Prentice, *The Statistical Analysis of Failure Time Data*; Klein and Moeschberger, *Survival Analysis*.
- Bootstrap and resampling: Efron and Tibshirani, *An Introduction to the Bootstrap*; Davison and Hinkley, *Bootstrap Methods and Their Application*.
- Randomization and causal inference: Fisher, *The Design of Experiments*; Rosenbaum, *Observational Studies*; Imbens and Rubin, *Causal Inference*.
- Quantile regression: Koenker, *Quantile Regression*.
- GEE and mixed models: Liang and Zeger papers for GEE; standard GLMM references for mixed-model likelihood and quadrature.

## Thinly Described API TODOs

Each numbered item below is an API entry whose generated description is thin or needs review. Edit the roxygen source, apply the General Instructions above, then regenerate Rd. The current description length is included as a triage signal; consult the generated Rd or source file for the existing wording. The `TODO #NNN` identifiers are stable range targets for instructions such as "do 41-50".


### `build_cox_data_cache_cpp.Rd`

- [ ] TODO #001: Topic `build_cox_data_cache_cpp` (current description `95` chars).

### `build_stratified_cox_data_cache_cpp.Rd`

- [ ] TODO #002: Topic `build_stratified_cox_data_cache_cpp` (current description `117` chars).

### `check_package_installed.Rd`

- [ ] TODO #003: Topic `check_package_installed` (current description `105` chars).

### `compute_coxph_rand_bootstrap_cpp.Rd`

- [ ] TODO #004: Topic `compute_coxph_rand_bootstrap_cpp` (current description `392` chars).

### `create_model_matrix_from_features.Rd`

- [ ] TODO #005: Topic `create_model_matrix_from_features` (current description `203` chars).

### `DesignFixedAOptimal.Rd`

- [ ] TODO #006: Method `DesignFixedAOptimal$clone()` (current description `57` chars).
- [ ] TODO #007: Method `DesignFixedAOptimal$new()` (current description `56` chars).
- [ ] TODO #008: Topic `DesignFixedAOptimal` (current description `302` chars).

### `DesignFixedBernoulli.Rd`

- [ ] TODO #009: Method `DesignFixedBernoulli$clone()` (current description `57` chars).
- [ ] TODO #010: Method `DesignFixedBernoulli$is_a_bernoulli_capable()` (current description `56` chars).
- [ ] TODO #011: Method `DesignFixedBernoulli$new()` (current description `48` chars).
- [ ] TODO #012: Topic `DesignFixedBernoulli` (current description `144` chars).

### `DesignFixedBinaryMatch.Rd`

- [ ] TODO #013: Method `DesignFixedBinaryMatch$assign_w_to_all_subjects()` (current description `72` chars).
- [ ] TODO #014: Method `DesignFixedBinaryMatch$clone()` (current description `57` chars).
- [ ] TODO #015: Method `DesignFixedBinaryMatch$is_a_kk_matching_capable()` (current description `66` chars).
- [ ] TODO #016: Method `DesignFixedBinaryMatch$new()` (current description `51` chars).
- [ ] TODO #017: Method `DesignFixedBinaryMatch$supports_batch_w_pregeneration()` (current description `57` chars).
- [ ] TODO #018: Topic `DesignFixedBinaryMatch` (current description `281` chars).

### `DesignFixedBlockedCluster.Rd`

- [ ] TODO #019: Method `DesignFixedBlockedCluster$clone()` (current description `57` chars).
- [ ] TODO #020: Method `DesignFixedBlockedCluster$is_a_cluster_capable()` (current description `54` chars).
- [ ] TODO #021: Method `DesignFixedBlockedCluster$new()` (current description `69` chars).
- [ ] TODO #022: Topic `DesignFixedBlockedCluster` (current description `247` chars).

### `DesignFixedBlocking.Rd`

- [ ] TODO #023: Method `DesignFixedBlocking$clone()` (current description `57` chars).
- [ ] TODO #024: Method `DesignFixedBlocking$new()` (current description `58` chars).
- [ ] TODO #025: Topic `DesignFixedBlocking` (current description `140` chars).

### `DesignFixedCluster.Rd`

- [ ] TODO #026: Method `DesignFixedCluster$clone()` (current description `57` chars).
- [ ] TODO #027: Method `DesignFixedCluster$is_a_cluster_capable()` (current description `54` chars).
- [ ] TODO #028: Method `DesignFixedCluster$new()` (current description `57` chars).
- [ ] TODO #029: Topic `DesignFixedCluster` (current description `236` chars).

### `DesignFixedDOptimal.Rd`

- [ ] TODO #030: Method `DesignFixedDOptimal$clone()` (current description `57` chars).
- [ ] TODO #031: Method `DesignFixedDOptimal$new()` (current description `55` chars).
- [ ] TODO #032: Topic `DesignFixedDOptimal` (current description `291` chars).

### `DesignFixedFactorial.Rd`

- [ ] TODO #033: Method `DesignFixedFactorial$assign_w_to_all_subjects()` (current description `54` chars).
- [ ] TODO #034: Method `DesignFixedFactorial$clone()` (current description `57` chars).
- [ ] TODO #035: Method `DesignFixedFactorial$draw_ws_according_to_design()` (current description `89` chars).
- [ ] TODO #036: Method `DesignFixedFactorial$get_w_factorial()` (current description `58` chars).
- [ ] TODO #037: Method `DesignFixedFactorial$get_w()` (current description `58` chars).
- [ ] TODO #038: Method `DesignFixedFactorial$new()` (current description `48` chars).
- [ ] TODO #039: Topic `DesignFixedFactorial` (current description `224` chars).

### `DesignFixedGreedy.Rd`

- [ ] TODO #040: Method `DesignFixedGreedy$clone()` (current description `57` chars).
- [ ] TODO #041: Method `DesignFixedGreedy$new()` (current description `52` chars).
- [ ] TODO #042: Method `DesignFixedGreedy$supports_batch_w_pregeneration()` (current description `70` chars).
- [ ] TODO #043: Topic `DesignFixedGreedy` (current description `196` chars).

### `DesignFixediBCRD.Rd`

- [ ] TODO #044: Method `DesignFixediBCRD$clone()` (current description `57` chars).
- [ ] TODO #045: Method `DesignFixediBCRD$new()` (current description `69` chars).
- [ ] TODO #046: Topic `DesignFixediBCRD` (current description `162` chars).

### `DesignFixedMatchingGreedyPairSwitching.Rd`

- [ ] TODO #047: Method `DesignFixedMatchingGreedyPairSwitching$clone()` (current description `57` chars).
- [ ] TODO #048: Method `DesignFixedMatchingGreedyPairSwitching$new()` (current description `90` chars).
- [ ] TODO #049: Method `DesignFixedMatchingGreedyPairSwitching$supports_batch_w_pregeneration()` (current description `70` chars).
- [ ] TODO #050: Topic `DesignFixedMatchingGreedyPairSwitching` (current description `254` chars).

### `DesignFixedNaiveMatch.Rd`

- [ ] TODO #051: Method `DesignFixedNaiveMatch$clone()` (current description `57` chars).
- [ ] TODO #052: Topic `DesignFixedNaiveMatch` (current description `623` chars).

### `DesignFixedOptimalBlocks.Rd`

- [ ] TODO #053: Method `DesignFixedOptimalBlocks$clone()` (current description `57` chars).
- [ ] TODO #054: Method `DesignFixedOptimalBlocks$new()` (current description `41` chars).
- [ ] TODO #055: Method `DesignFixedOptimalBlocks$supports_batch_w_pregeneration()` (current description `57` chars).
- [ ] TODO #056: Topic `DesignFixedOptimalBlocks` (current description `395` chars).

### `DesignFixedRerandomization.Rd`

- [ ] TODO #057: Method `DesignFixedRerandomization$clone()` (current description `57` chars).
- [ ] TODO #058: Method `DesignFixedRerandomization$new()` (current description `54` chars).
- [ ] TODO #059: Topic `DesignFixedRerandomization` (current description `327` chars).

### `DesignSeqOneByOneAtkinson.Rd`

- [ ] TODO #060: Method `DesignSeqOneByOneAtkinson$assign_wt()` (current description `52` chars).
- [ ] TODO #061: Method `DesignSeqOneByOneAtkinson$clone()` (current description `57` chars).
- [ ] TODO #062: Method `DesignSeqOneByOneAtkinson$new()` (current description `53` chars).
- [ ] TODO #063: Topic `DesignSeqOneByOneAtkinson` (current description `230` chars).

### `DesignSeqOneByOneBernoulli.Rd`

- [ ] TODO #064: Method `DesignSeqOneByOneBernoulli$assign_wt()` (current description `52` chars).
- [ ] TODO #065: Method `DesignSeqOneByOneBernoulli$clone()` (current description `57` chars).
- [ ] TODO #066: Method `DesignSeqOneByOneBernoulli$is_a_bernoulli_capable()` (current description `56` chars).
- [ ] TODO #067: Method `DesignSeqOneByOneBernoulli$new()` (current description `53` chars).
- [ ] TODO #068: Topic `DesignSeqOneByOneBernoulli` (current description `144` chars).

### `DesignSeqOneByOneEfron.Rd`

- [ ] TODO #069: Method `DesignSeqOneByOneEfron$assign_wt()` (current description `52` chars).
- [ ] TODO #070: Method `DesignSeqOneByOneEfron$clone()` (current description `57` chars).
- [ ] TODO #071: Method `DesignSeqOneByOneEfron$new()` (current description `62` chars).
- [ ] TODO #072: Topic `DesignSeqOneByOneEfron` (current description `153` chars).

### `DesignSeqOneByOneiBCRD.Rd`

- [ ] TODO #073: Method `DesignSeqOneByOneiBCRD$add_one_subject_to_experiment_and_assign()` (current description `37` chars).
- [ ] TODO #074: Method `DesignSeqOneByOneiBCRD$assign_wt()` (current description `52` chars).
- [ ] TODO #075: Method `DesignSeqOneByOneiBCRD$clone()` (current description `57` chars).
- [ ] TODO #076: Method `DesignSeqOneByOneiBCRD$new()` (current description `52` chars).
- [ ] TODO #077: Topic `DesignSeqOneByOneiBCRD` (current description `228` chars).

### `DesignSeqOneByOneKK14.Rd`

- [ ] TODO #078: Method `DesignSeqOneByOneKK14$assign_wt()` (current description `52` chars).
- [ ] TODO #079: Method `DesignSeqOneByOneKK14$clone()` (current description `57` chars).
- [ ] TODO #080: Method `DesignSeqOneByOneKK14$is_a_kk_matching_capable()` (current description `66` chars).
- [ ] TODO #081: Method `DesignSeqOneByOneKK14$new()` (current description `48` chars).
- [ ] TODO #082: Topic `DesignSeqOneByOneKK14` (current description `282` chars).

### `DesignSeqOneByOneKK21.Rd`

- [ ] TODO #083: Method `DesignSeqOneByOneKK21$assign_wt()` (current description `70` chars).
- [ ] TODO #084: Method `DesignSeqOneByOneKK21$clone()` (current description `57` chars).
- [ ] TODO #085: Method `DesignSeqOneByOneKK21$get_covariate_weights()` (current description `62` chars).
- [ ] TODO #086: Method `DesignSeqOneByOneKK21$get_iteration_weights()` (current description `49` chars).
- [ ] TODO #087: Method `DesignSeqOneByOneKK21$new()` (current description `119` chars).
- [ ] TODO #088: Topic `DesignSeqOneByOneKK21` (current description `282` chars).

### `DesignSeqOneByOneKK21stepwise.Rd`

- [ ] TODO #089: Method `DesignSeqOneByOneKK21stepwise$clone()` (current description `57` chars).
- [ ] TODO #090: Method `DesignSeqOneByOneKK21stepwise$new()` (current description `198` chars).
- [ ] TODO #091: Topic `DesignSeqOneByOneKK21stepwise` (current description `282` chars).

### `DesignSeqOneByOnePocockSimon.Rd`

- [ ] TODO #092: Method `DesignSeqOneByOnePocockSimon$assign_wt()` (current description `64` chars).
- [ ] TODO #093: Method `DesignSeqOneByOnePocockSimon$clone()` (current description `57` chars).
- [ ] TODO #094: Method `DesignSeqOneByOnePocockSimon$new()` (current description `58` chars).
- [ ] TODO #095: Topic `DesignSeqOneByOnePocockSimon` (current description `232` chars).

### `DesignSeqOneByOneRandomBlockSize.Rd`

- [ ] TODO #096: Method `DesignSeqOneByOneRandomBlockSize$assign_wt()` (current description `52` chars).
- [ ] TODO #097: Method `DesignSeqOneByOneRandomBlockSize$clone()` (current description `57` chars).
- [ ] TODO #098: Method `DesignSeqOneByOneRandomBlockSize$new()` (current description `61` chars).
- [ ] TODO #099: Topic `DesignSeqOneByOneRandomBlockSize` (current description `339` chars).

### `DesignSeqOneByOneSPBR.Rd`

- [ ] TODO #100: Method `DesignSeqOneByOneSPBR$assign_wt()` (current description `52` chars).
- [ ] TODO #101: Method `DesignSeqOneByOneSPBR$clone()` (current description `57` chars).
- [ ] TODO #102: Method `DesignSeqOneByOneSPBR$new()` (current description `69` chars).
- [ ] TODO #103: Topic `DesignSeqOneByOneSPBR` (current description `251` chars).

### `DesignSeqOneByOneUrn.Rd`

- [ ] TODO #104: Method `DesignSeqOneByOneUrn$assign_wt()` (current description `52` chars).
- [ ] TODO #105: Method `DesignSeqOneByOneUrn$clone()` (current description `57` chars).
- [ ] TODO #106: Method `DesignSeqOneByOneUrn$new()` (current description `48` chars).
- [ ] TODO #107: Topic `DesignSeqOneByOneUrn` (current description `273` chars).

### `dot-normalize_optimizer_algorithm.Rd`

- [ ] TODO #108: Topic `.normalize_optimizer_algorithm` (current description `156` chars).

### `edi_build_info_cpp.Rd`

- [ ] TODO #109: Topic `edi_build_info_cpp` (current description `247` chars).

### `exact_jonckheere_terpstra_pval_cpp.Rd`

- [ ] TODO #110: Topic `exact_jonckheere_terpstra_pval_cpp` (current description `115` chars).

### `expand_adjacent_category_data_cpp.Rd`

- [ ] TODO #111: Topic `expand_adjacent_category_data_cpp` (current description `106` chars).

### `expand_continuation_ratio_data_cpp.Rd`

- [ ] TODO #112: Topic `expand_continuation_ratio_data_cpp` (current description `102` chars).

### `fast_adjacent_category_logit_cpp.Rd`

- [ ] TODO #113: Topic `fast_adjacent_category_logit_cpp` (current description `90` chars).

### `fast_adjacent_category_logit_with_var_cpp.Rd`

- [ ] TODO #114: Topic `fast_adjacent_category_logit_with_var_cpp` (current description `124` chars).

### `fast_beta_regression_cpp.Rd`

- [ ] TODO #115: Topic `fast_beta_regression_cpp` (current description `107` chars).

### `fast_beta_regression_weighted_cpp.Rd`

- [ ] TODO #116: Topic `fast_beta_regression_weighted_cpp` (current description `106` chars).

### `fast_beta_regression_with_var_cpp.Rd`

- [ ] TODO #117: Topic `fast_beta_regression_with_var_cpp` (current description `144` chars).

### `fast_beta_regression_with_var.Rd`

- [ ] TODO #118: Topic `fast_beta_regression_with_var` (current description `369` chars).

### `fast_beta_regression.Rd`

- [ ] TODO #119: Topic `fast_beta_regression` (current description `277` chars).

### `fast_continuation_ratio_regression_cpp.Rd`

- [ ] TODO #120: Topic `fast_continuation_ratio_regression_cpp` (current description `97` chars).

### `fast_continuation_ratio_regression_with_var_cpp.Rd`

- [ ] TODO #121: Topic `fast_continuation_ratio_regression_with_var_cpp` (current description `141` chars).

### `fast_coxph_regression_cpp.Rd`

- [ ] TODO #122: Topic `fast_coxph_regression_cpp` (current description `109` chars).

### `fast_coxph_regression_prebuilt_cpp.Rd`

- [ ] TODO #123: Topic `fast_coxph_regression_prebuilt_cpp` (current description `115` chars).

### `fast_coxph_regression.Rd`

- [ ] TODO #124: Topic `fast_coxph_regression` (current description `264` chars).

### `fast_cpoisson_combined_with_var_cpp.Rd`

- [ ] TODO #125: Topic `fast_cpoisson_combined_with_var_cpp` (current description `183` chars).

### `fast_digamma_vec_cpp.Rd`

- [ ] TODO #126: Topic `fast_digamma_vec_cpp` (current description `144` chars).

### `fast_dnbinom_mu_vec_cpp.Rd`

- [ ] TODO #127: Topic `fast_dnbinom_mu_vec_cpp` (current description `295` chars).

### `fast_hurdle_negbin_with_var_cpp.Rd`

- [ ] TODO #128: Topic `fast_hurdle_negbin_with_var_cpp` (current description `119` chars).

### `fast_identity_binomial_regression_cpp.Rd`

- [ ] TODO #129: Topic `fast_identity_binomial_regression_cpp` (current description `117` chars).

### `fast_identity_binomial_regression_weighted_cpp.Rd`

- [ ] TODO #130: Topic `fast_identity_binomial_regression_weighted_cpp` (current description `135` chars).

### `fast_identity_binomial_regression_with_var_cpp.Rd`

- [ ] TODO #131: Topic `fast_identity_binomial_regression_with_var_cpp` (current description `151` chars).

### `fast_lbeta_vec_cpp.Rd`

- [ ] TODO #132: Topic `fast_lbeta_vec_cpp` (current description `153` chars).

### `fast_lgamma_vec_cpp.Rd`

- [ ] TODO #133: Topic `fast_lgamma_vec_cpp` (current description `128` chars).

### `fast_log_binomial_regression_cpp.Rd`

- [ ] TODO #134: Topic `fast_log_binomial_regression_cpp` (current description `105` chars).

### `fast_log_binomial_regression_weighted_cpp.Rd`

- [ ] TODO #135: Topic `fast_log_binomial_regression_weighted_cpp` (current description `123` chars).

### `fast_log_binomial_regression_with_var_cpp.Rd`

- [ ] TODO #136: Topic `fast_log_binomial_regression_with_var_cpp` (current description `125` chars).

### `fast_log_dnorm_vec_cpp.Rd`

- [ ] TODO #137: Topic `fast_log_dnorm_vec_cpp` (current description `123` chars).

### `fast_log_pnorm_vec_cpp.Rd`

- [ ] TODO #138: Topic `fast_log_pnorm_vec_cpp` (current description `162` chars).

### `fast_logistic_regression_cpp.Rd`

- [ ] TODO #139: Topic `fast_logistic_regression_cpp` (current description `77` chars).

### `fast_logistic_regression_weighted_cpp.Rd`

- [ ] TODO #140: Topic `fast_logistic_regression_weighted_cpp` (current description `95` chars).

### `fast_logistic_regression_with_var_cpp.Rd`

- [ ] TODO #141: Topic `fast_logistic_regression_with_var_cpp` (current description `121` chars).

### `fast_logistic_regression_with_var.Rd`

- [ ] TODO #142: Topic `fast_logistic_regression_with_var` (current description `315` chars).

### `fast_logistic_regression.Rd`

- [ ] TODO #143: Topic `fast_logistic_regression` (current description `262` chars).

### `fast_neg_bin_cpp.Rd`

- [ ] TODO #144: Topic `fast_neg_bin_cpp` (current description `102` chars).

### `fast_neg_bin_weighted_cpp.Rd`

- [ ] TODO #145: Topic `fast_neg_bin_weighted_cpp` (current description `132` chars).

### `fast_neg_bin_with_var_cpp.Rd`

- [ ] TODO #146: Topic `fast_neg_bin_with_var_cpp` (current description `148` chars).

### `fast_negbin_regression_with_var.Rd`

- [ ] TODO #147: Topic `fast_negbin_regression_with_var` (current description `406` chars).

### `fast_negbin_regression.Rd`

- [ ] TODO #148: Topic `fast_negbin_regression` (current description `302` chars).

### `fast_ols_cpp.Rd`

- [ ] TODO #149: Topic `fast_ols_cpp` (current description `369` chars).

### `fast_ols_with_var_cpp.Rd`

- [ ] TODO #150: Topic `fast_ols_with_var_cpp` (current description `282` chars).

### `fast_ordinal_cauchit_regression_cpp.Rd`

- [ ] TODO #151: Topic `fast_ordinal_cauchit_regression_cpp` (current description `111` chars).

### `fast_ordinal_cauchit_regression_with_var_cpp.Rd`

- [ ] TODO #152: Topic `fast_ordinal_cauchit_regression_with_var_cpp` (current description `124` chars).

### `fast_ordinal_cloglog_regression_cpp.Rd`

- [ ] TODO #153: Topic `fast_ordinal_cloglog_regression_cpp` (current description `111` chars).

### `fast_ordinal_cloglog_regression_with_var_cpp.Rd`

- [ ] TODO #154: Topic `fast_ordinal_cloglog_regression_with_var_cpp` (current description `124` chars).

### `fast_ordinal_glmm_cpp.Rd`

- [ ] TODO #155: Topic `fast_ordinal_glmm_cpp` (current description `121` chars).

### `fast_ordinal_probit_regression_cpp.Rd`

- [ ] TODO #156: Topic `fast_ordinal_probit_regression_cpp` (current description `109` chars).

### `fast_ordinal_probit_regression_with_var_cpp.Rd`

- [ ] TODO #157: Topic `fast_ordinal_probit_regression_with_var_cpp` (current description `122` chars).

### `fast_ordinal_regression_cpp.Rd`

- [ ] TODO #158: Topic `fast_ordinal_regression_cpp` (current description `95` chars).

### `fast_ordinal_regression_weighted_cpp.Rd`

- [ ] TODO #159: Topic `fast_ordinal_regression_weighted_cpp` (current description `113` chars).

### `fast_ordinal_regression_with_var_cpp.Rd`

- [ ] TODO #160: Topic `fast_ordinal_regression_with_var_cpp` (current description `108` chars).

### `fast_poisson_regression_cpp.Rd`

- [ ] TODO #161: Topic `fast_poisson_regression_cpp` (current description `208` chars).

### `fast_poisson_regression_weighted_cpp.Rd`

- [ ] TODO #162: Topic `fast_poisson_regression_weighted_cpp` (current description `93` chars).

### `fast_poisson_regression_with_var_cpp.Rd`

- [ ] TODO #163: Topic `fast_poisson_regression_with_var_cpp` (current description `127` chars).

### `fast_probit_regression_cpp.Rd`

- [ ] TODO #164: Topic `fast_probit_regression_cpp` (current description `84` chars).

### `fast_probit_regression_with_var_cpp.Rd`

- [ ] TODO #165: Topic `fast_probit_regression_with_var_cpp` (current description `117` chars).

### `fast_qnorm_vec_cpp.Rd`

- [ ] TODO #166: Topic `fast_qnorm_vec_cpp` (current description `200` chars).

### `fast_quasipoisson_regression_with_var_cpp.Rd`

- [ ] TODO #167: Topic `fast_quasipoisson_regression_with_var_cpp` (current description `139` chars).

### `fast_robust_regression_cpp.Rd`

- [ ] TODO #168: Topic `fast_robust_regression_cpp` (current description `83` chars).

### `fast_stereotype_logit_cpp.Rd`

- [ ] TODO #169: Topic `fast_stereotype_logit_cpp` (current description `113` chars).

### `fast_stereotype_logit_with_var_cpp.Rd`

- [ ] TODO #170: Topic `fast_stereotype_logit_with_var_cpp` (current description `126` chars).

### `fast_stereotype_profile_loglik_cpp.Rd`

- [ ] TODO #171: Topic `fast_stereotype_profile_loglik_cpp` (current description `131` chars).

### `fast_trigamma_vec_cpp.Rd`

- [ ] TODO #172: Topic `fast_trigamma_vec_cpp` (current description `148` chars).

### `fast_weibull_regression_cpp.Rd`

- [ ] TODO #173: Topic `fast_weibull_regression_cpp` (current description `85` chars).

### `fast_weibull_regression.Rd`

- [ ] TODO #174: Topic `fast_weibull_regression` (current description `167` chars).

### `fast_zero_augmented_poisson_cpp.Rd`

- [ ] TODO #175: Topic `fast_zero_augmented_poisson_cpp` (current description `134` chars).

### `fast_zero_one_inflated_beta_cpp.Rd`

- [ ] TODO #176: Topic `fast_zero_one_inflated_beta_cpp` (current description `135` chars).

### `gcomp_fractional_logit_point_estimate_cpp.Rd`

- [ ] TODO #177: Topic `gcomp_fractional_logit_point_estimate_cpp` (current description `165` chars).

### `gcomp_logistic_point_estimate_cpp.Rd`

- [ ] TODO #178: Topic `gcomp_logistic_point_estimate_cpp` (current description `158` chars).

### `gcomp_logistic_post_fit_cpp.Rd`

- [ ] TODO #179: Topic `gcomp_logistic_post_fit_cpp` (current description `101` chars).

### `gcomp_ordinal_proportional_odds_post_fit_cpp.Rd`

- [ ] TODO #180: Topic `gcomp_ordinal_proportional_odds_post_fit_cpp` (current description `135` chars).

### `generate_covariate_dataset.Rd`

- [ ] TODO #181: Topic `generate_covariate_dataset` (current description `330` chars).

### `get_beta_regression_hessian_cpp.Rd`

- [ ] TODO #182: Topic `get_beta_regression_hessian_cpp` (current description `133` chars).

### `get_bootstrap_dispatch_policy.Rd`

- [ ] TODO #183: Topic `get_bootstrap_dispatch_policy` (current description `205` chars).

### `get_cold_start_dispatch_policy.Rd`

- [ ] TODO #184: Topic `get_cold_start_dispatch_policy` (current description `462` chars).

### `get_cpoisson_combined_hessian_cpp.Rd`

- [ ] TODO #185: Topic `get_cpoisson_combined_hessian_cpp` (current description `130` chars).

### `get_cpoisson_combined_score_cpp.Rd`

- [ ] TODO #186: Topic `get_cpoisson_combined_score_cpp` (current description `135` chars).

### `get_identity_binomial_regression_hessian_cpp.Rd`

- [ ] TODO #187: Topic `get_identity_binomial_regression_hessian_cpp` (current description `129` chars).

### `get_identity_binomial_regression_score_cpp.Rd`

- [ ] TODO #188: Topic `get_identity_binomial_regression_score_cpp` (current description `125` chars).

### `get_identity_binomial_regression_weighted_hessian_cpp.Rd`

- [ ] TODO #189: Topic `get_identity_binomial_regression_weighted_hessian_cpp` (current description `147` chars).

### `get_identity_binomial_regression_weighted_score_cpp.Rd`

- [ ] TODO #190: Topic `get_identity_binomial_regression_weighted_score_cpp` (current description `143` chars).

### `get_log_binomial_regression_hessian_cpp.Rd`

- [ ] TODO #191: Topic `get_log_binomial_regression_hessian_cpp` (current description `149` chars).

### `get_log_binomial_regression_score_cpp.Rd`

- [ ] TODO #192: Topic `get_log_binomial_regression_score_cpp` (current description `135` chars).

### `get_log_binomial_regression_weighted_hessian_cpp.Rd`

- [ ] TODO #193: Topic `get_log_binomial_regression_weighted_hessian_cpp` (current description `124` chars).

### `get_log_binomial_regression_weighted_score_cpp.Rd`

- [ ] TODO #194: Topic `get_log_binomial_regression_weighted_score_cpp` (current description `120` chars).

### `get_negbin_regression_hessian_cpp.Rd`

- [ ] TODO #195: Topic `get_negbin_regression_hessian_cpp` (current description `159` chars).

### `get_optimization_dispatch_policy.Rd`

- [ ] TODO #196: Topic `get_optimization_dispatch_policy` (current description `121` chars).

### `get_ordinal_regression_hessian_cpp.Rd`

- [ ] TODO #197: Topic `get_ordinal_regression_hessian_cpp` (current description `140` chars).

### `get_ordinal_regression_score_cpp.Rd`

- [ ] TODO #198: Topic `get_ordinal_regression_score_cpp` (current description `126` chars).

### `get_parallel_dispatch_policy.Rd`

- [ ] TODO #199: Topic `get_parallel_dispatch_policy` (current description `115` chars).

### `get_stereotype_logit_hessian_cpp.Rd`

- [ ] TODO #200: Topic `get_stereotype_logit_hessian_cpp` (current description `135` chars).

### `get_warm_start_dispatch_policy.Rd`

- [ ] TODO #201: Topic `get_warm_start_dispatch_policy` (current description `171` chars).

### `get_weibull_regression_hessian_cpp.Rd`

- [ ] TODO #202: Topic `get_weibull_regression_hessian_cpp` (current description `143` chars).

### `get_weibull_regression_score_cpp.Rd`

- [ ] TODO #203: Topic `get_weibull_regression_score_cpp` (current description `129` chars).

### `InferenceAllKKMeanDiffIVWC.Rd`

- [ ] TODO #204: Method `InferenceAllKKMeanDiffIVWC$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #205: Method `InferenceAllKKMeanDiffIVWC$clone()` (current description `57` chars).
- [ ] TODO #206: Method `InferenceAllKKMeanDiffIVWC$compute_asymp_confidence_interval()` (current description `198` chars).
- [ ] TODO #207: Method `InferenceAllKKMeanDiffIVWC$compute_asymp_two_sided_pval()` (current description `144` chars).
- [ ] TODO #208: Method `InferenceAllKKMeanDiffIVWC$compute_estimate()` (current description `69` chars).
- [ ] TODO #209: Method `InferenceAllKKMeanDiffIVWC$compute_rand_confidence_interval()` (current description `63` chars).
- [ ] TODO #210: Method `InferenceAllKKMeanDiffIVWC$get_likelihood_test_spec()` (current description `155` chars).
- [ ] TODO #211: Method `InferenceAllKKMeanDiffIVWC$simulate_under_lik_null()` (current description `49` chars).
- [ ] TODO #212: Method `InferenceAllKKMeanDiffIVWC$supports_lik_ratio_param_bootstrap()` (current description `59` chars).
- [ ] TODO #213: Method `InferenceAllKKMeanDiffIVWC$supports_likelihood_tests()` (current description `45` chars).
- [ ] TODO #214: Topic `InferenceAllKKMeanDiffIVWC` (current description `385` chars).

### `InferenceAllKKWilcoxIVWC.Rd`

- [ ] TODO #215: Method `InferenceAllKKWilcoxIVWC$clone()` (current description `57` chars).
- [ ] TODO #216: Method `InferenceAllKKWilcoxIVWC$compute_asymp_confidence_interval()` (current description `48` chars).
- [ ] TODO #217: Method `InferenceAllKKWilcoxIVWC$compute_asymp_two_sided_pval()` (current description `125` chars).
- [ ] TODO #218: Method `InferenceAllKKWilcoxIVWC$compute_estimate()` (current description `69` chars).
- [ ] TODO #219: Method `InferenceAllKKWilcoxIVWC$compute_jackknife_bias_estimate()` (current description `130` chars).
- [ ] TODO #220: Method `InferenceAllKKWilcoxIVWC$compute_jackknife_estimate()` (current description `45` chars).
- [ ] TODO #221: Method `InferenceAllKKWilcoxIVWC$compute_jackknife_std_error()` (current description `131` chars).
- [ ] TODO #222: Method `InferenceAllKKWilcoxIVWC$compute_jackknife_wald_confidence_interval()` (current description `63` chars).
- [ ] TODO #223: Method `InferenceAllKKWilcoxIVWC$compute_jackknife_wald_two_sided_pval()` (current description `62` chars).
- [ ] TODO #224: Method `InferenceAllKKWilcoxIVWC$new()` (current description `118` chars).
- [ ] TODO #225: Topic `InferenceAllKKWilcoxIVWC` (current description `495` chars).

### `InferenceAllSimpleMeanDiff.Rd`

- [ ] TODO #226: Method `InferenceAllSimpleMeanDiff$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #227: Method `InferenceAllSimpleMeanDiff$clone()` (current description `57` chars).
- [ ] TODO #228: Method `InferenceAllSimpleMeanDiff$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #229: Method `InferenceAllSimpleMeanDiff$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #230: Method `InferenceAllSimpleMeanDiff$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #231: Method `InferenceAllSimpleMeanDiff$compute_estimate()` (current description `58` chars).
- [ ] TODO #232: Method `InferenceAllSimpleMeanDiff$new()` (current description `53` chars).
- [ ] TODO #233: Topic `InferenceAllSimpleMeanDiff` (current description `285` chars).

### `InferenceAllSimpleMeanDiffPooledVar.Rd`

- [ ] TODO #234: Method `InferenceAllSimpleMeanDiffPooledVar$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #235: Method `InferenceAllSimpleMeanDiffPooledVar$clone()` (current description `57` chars).
- [ ] TODO #236: Method `InferenceAllSimpleMeanDiffPooledVar$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #237: Method `InferenceAllSimpleMeanDiffPooledVar$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #238: Method `InferenceAllSimpleMeanDiffPooledVar$compute_estimate()` (current description `58` chars).
- [ ] TODO #239: Method `InferenceAllSimpleMeanDiffPooledVar$new()` (current description `142` chars).
- [ ] TODO #240: Topic `InferenceAllSimpleMeanDiffPooledVar` (current description `329` chars).

### `InferenceAllSimpleWilcox.Rd`

- [ ] TODO #241: Method `InferenceAllSimpleWilcox$clone()` (current description `57` chars).
- [ ] TODO #242: Method `InferenceAllSimpleWilcox$compute_asymp_confidence_interval()` (current description `56` chars).
- [ ] TODO #243: Method `InferenceAllSimpleWilcox$compute_asymp_two_sided_pval()` (current description `64` chars).
- [ ] TODO #244: Method `InferenceAllSimpleWilcox$compute_estimate_with_bootstrap_weights()` (current description `62` chars).
- [ ] TODO #245: Method `InferenceAllSimpleWilcox$compute_estimate()` (current description `80` chars).
- [ ] TODO #246: Method `InferenceAllSimpleWilcox$compute_jackknife_bias_estimate()` (current description `128` chars).
- [ ] TODO #247: Method `InferenceAllSimpleWilcox$compute_jackknife_estimate()` (current description `45` chars).
- [ ] TODO #248: Method `InferenceAllSimpleWilcox$compute_jackknife_std_error()` (current description `129` chars).
- [ ] TODO #249: Method `InferenceAllSimpleWilcox$compute_jackknife_wald_confidence_interval()` (current description `63` chars).
- [ ] TODO #250: Method `InferenceAllSimpleWilcox$compute_jackknife_wald_two_sided_pval()` (current description `62` chars).
- [ ] TODO #251: Method `InferenceAllSimpleWilcox$new()` (current description `91` chars).
- [ ] TODO #252: Topic `InferenceAllSimpleWilcox` (current description `135` chars).

### `InferenceBaiAdjustedTKK14.Rd`

- [ ] TODO #253: Method `InferenceBaiAdjustedTKK14$clone()` (current description `57` chars).
- [ ] TODO #254: Topic `InferenceBaiAdjustedTKK14` (current description `229` chars).

### `InferenceBaiAdjustedTKK21.Rd`

- [ ] TODO #255: Method `InferenceBaiAdjustedTKK21$clone()` (current description `57` chars).
- [ ] TODO #256: Topic `InferenceBaiAdjustedTKK21` (current description `229` chars).

### `InferenceContinKKGLMM.Rd`

- [ ] TODO #257: Method `InferenceContinKKGLMM$clone()` (current description `57` chars).
- [ ] TODO #258: Method `InferenceContinKKGLMM$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #259: Method `InferenceContinKKGLMM$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #260: Method `InferenceContinKKGLMM$compute_estimate_with_bootstrap_weights()` (current description `62` chars).
- [ ] TODO #261: Method `InferenceContinKKGLMM$compute_estimate()` (current description `58` chars).
- [ ] TODO #262: Method `InferenceContinKKGLMM$new()` (current description `38` chars).
- [ ] TODO #263: Topic `InferenceContinKKGLMM` (current description `456` chars).

### `InferenceContinKKOLSIVWC.Rd`

- [ ] TODO #264: Method `InferenceContinKKOLSIVWC$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #265: Method `InferenceContinKKOLSIVWC$clone()` (current description `57` chars).
- [ ] TODO #266: Method `InferenceContinKKOLSIVWC$new()` (current description `110` chars).
- [ ] TODO #267: Topic `InferenceContinKKOLSIVWC` (current description `432` chars).

### `InferenceContinKKOLSOneLik.Rd`

- [ ] TODO #268: Method `InferenceContinKKOLSOneLik$clone()` (current description `57` chars).
- [ ] TODO #269: Method `InferenceContinKKOLSOneLik$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #270: Method `InferenceContinKKOLSOneLik$compute_asymp_two_sided_pval()` (current description `107` chars).
- [ ] TODO #271: Method `InferenceContinKKOLSOneLik$compute_estimate_with_bootstrap_weights()` (current description `62` chars).
- [ ] TODO #272: Method `InferenceContinKKOLSOneLik$compute_estimate()` (current description `58` chars).
- [ ] TODO #273: Method `InferenceContinKKOLSOneLik$get_likelihood_components()` (current description `74` chars).
- [ ] TODO #274: Method `InferenceContinKKOLSOneLik$new()` (current description `111` chars).
- [ ] TODO #275: Topic `InferenceContinKKOLSOneLik` (current description `416` chars).

### `InferenceContinKKQuantileRegrIVWC.Rd`

- [ ] TODO #276: Method `InferenceContinKKQuantileRegrIVWC$clone()` (current description `57` chars).
- [ ] TODO #277: Method `InferenceContinKKQuantileRegrIVWC$new()` (current description `73` chars).
- [ ] TODO #278: Topic `InferenceContinKKQuantileRegrIVWC` (current description `1586` chars).

### `InferenceContinKKQuantileRegrOneLik.Rd`

- [ ] TODO #279: Method `InferenceContinKKQuantileRegrOneLik$clone()` (current description `57` chars).
- [ ] TODO #280: Method `InferenceContinKKQuantileRegrOneLik$compute_estimate()` (current description `57` chars).
- [ ] TODO #281: Method `InferenceContinKKQuantileRegrOneLik$new()` (current description `144` chars).
- [ ] TODO #282: Topic `InferenceContinKKQuantileRegrOneLik` (current description `413` chars).

### `InferenceContinKKRobustRegrIVWC.Rd`

- [ ] TODO #283: Method `InferenceContinKKRobustRegrIVWC$clone()` (current description `57` chars).
- [ ] TODO #284: Method `InferenceContinKKRobustRegrIVWC$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #285: Method `InferenceContinKKRobustRegrIVWC$compute_asymp_two_sided_pval()` (current description `106` chars).
- [ ] TODO #286: Method `InferenceContinKKRobustRegrIVWC$compute_estimate()` (current description `57` chars).
- [ ] TODO #287: Method `InferenceContinKKRobustRegrIVWC$duplicate()` (current description `164` chars).
- [ ] TODO #288: Method `InferenceContinKKRobustRegrIVWC$new()` (current description `120` chars).
- [ ] TODO #289: Topic `InferenceContinKKRobustRegrIVWC` (current description `321` chars).

### `InferenceContinKKRobustRegrOneLik.Rd`

- [ ] TODO #290: Method `InferenceContinKKRobustRegrOneLik$clone()` (current description `57` chars).
- [ ] TODO #291: Method `InferenceContinKKRobustRegrOneLik$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #292: Method `InferenceContinKKRobustRegrOneLik$compute_asymp_two_sided_pval()` (current description `121` chars).
- [ ] TODO #293: Method `InferenceContinKKRobustRegrOneLik$compute_estimate_with_bootstrap_weights()` (current description `56` chars).
- [ ] TODO #294: Method `InferenceContinKKRobustRegrOneLik$compute_estimate()` (current description `72` chars).
- [ ] TODO #295: Method `InferenceContinKKRobustRegrOneLik$compute_wald_confidence_interval()` (current description `38` chars).
- [ ] TODO #296: Method `InferenceContinKKRobustRegrOneLik$compute_wald_two_sided_pval()` (current description `36` chars).
- [ ] TODO #297: Method `InferenceContinKKRobustRegrOneLik$duplicate()` (current description `172` chars).
- [ ] TODO #298: Method `InferenceContinKKRobustRegrOneLik$new()` (current description `114` chars).
- [ ] TODO #299: Topic `InferenceContinKKRobustRegrOneLik` (current description `221` chars).

### `InferenceContinLin.Rd`

- [ ] TODO #300: Method `InferenceContinLin$clone()` (current description `57` chars).
- [ ] TODO #301: Method `InferenceContinLin$compute_asymp_confidence_interval()` (current description `73` chars).
- [ ] TODO #302: Method `InferenceContinLin$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #303: Method `InferenceContinLin$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #304: Method `InferenceContinLin$compute_estimate()` (current description `48` chars).
- [ ] TODO #305: Method `InferenceContinLin$new()` (current description `41` chars).
- [ ] TODO #306: Topic `InferenceContinLin` (current description `360` chars).

### `InferenceContinOLS.Rd`

- [ ] TODO #307: Method `InferenceContinOLS$clone()` (current description `57` chars).
- [ ] TODO #308: Method `InferenceContinOLS$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #309: Method `InferenceContinOLS$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #310: Method `InferenceContinOLS$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #311: Method `InferenceContinOLS$compute_estimate()` (current description `50` chars).
- [ ] TODO #312: Method `InferenceContinOLS$new()` (current description `35` chars).
- [ ] TODO #313: Topic `InferenceContinOLS` (current description `317` chars).

### `InferenceContinQuantileRegr.Rd`

- [ ] TODO #314: Method `InferenceContinQuantileRegr$clone()` (current description `57` chars).
- [ ] TODO #315: Method `InferenceContinQuantileRegr$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #316: Method `InferenceContinQuantileRegr$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #317: Method `InferenceContinQuantileRegr$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #318: Method `InferenceContinQuantileRegr$compute_estimate()` (current description `66` chars).
- [ ] TODO #319: Method `InferenceContinQuantileRegr$new()` (current description `72` chars).
- [ ] TODO #320: Topic `InferenceContinQuantileRegr` (current description `651` chars).

### `InferenceContinRobustRegr.Rd`

- [ ] TODO #321: Method `InferenceContinRobustRegr$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #322: Method `InferenceContinRobustRegr$clone()` (current description `57` chars).
- [ ] TODO #323: Method `InferenceContinRobustRegr$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #324: Method `InferenceContinRobustRegr$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #325: Method `InferenceContinRobustRegr$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #326: Method `InferenceContinRobustRegr$compute_estimate()` (current description `64` chars).
- [ ] TODO #327: Method `InferenceContinRobustRegr$new()` (current description `70` chars).
- [ ] TODO #328: Topic `InferenceContinRobustRegr` (current description `650` chars).

### `InferenceCountHurdleNegBin.Rd`

- [ ] TODO #329: Method `InferenceCountHurdleNegBin$clone()` (current description `57` chars).
- [ ] TODO #330: Method `InferenceCountHurdleNegBin$compute_asymp_confidence_interval()` (current description `156` chars).
- [ ] TODO #331: Method `InferenceCountHurdleNegBin$compute_asymp_two_sided_pval()` (current description `172` chars).
- [ ] TODO #332: Method `InferenceCountHurdleNegBin$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #333: Method `InferenceCountHurdleNegBin$compute_gradient_confidence_interval()` (current description `145` chars).
- [ ] TODO #334: Method `InferenceCountHurdleNegBin$compute_gradient_two_sided_pval()` (current description `48` chars).
- [ ] TODO #335: Method `InferenceCountHurdleNegBin$compute_jackknife_bias_estimate()` (current description `166` chars).
- [ ] TODO #336: Method `InferenceCountHurdleNegBin$compute_jackknife_estimate()` (current description `59` chars).
- [ ] TODO #337: Method `InferenceCountHurdleNegBin$compute_jackknife_std_error()` (current description `167` chars).
- [ ] TODO #338: Method `InferenceCountHurdleNegBin$compute_jackknife_wald_confidence_interval()` (current description `63` chars).
- [ ] TODO #339: Method `InferenceCountHurdleNegBin$compute_jackknife_wald_two_sided_pval()` (current description `62` chars).
- [ ] TODO #340: Method `InferenceCountHurdleNegBin$new()` (current description `154` chars).
- [ ] TODO #341: Topic `InferenceCountHurdleNegBin` (current description `214` chars).

### `InferenceCountHurdlePoisson.Rd`

- [ ] TODO #342: Method `InferenceCountHurdlePoisson$clone()` (current description `57` chars).
- [ ] TODO #343: Method `InferenceCountHurdlePoisson$new()` (current description `45` chars).
- [ ] TODO #344: Topic `InferenceCountHurdlePoisson` (current description `194` chars).

### `InferenceCountKKCondPoissonOneLik.Rd`

- [ ] TODO #345: Method `InferenceCountKKCondPoissonOneLik$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #346: Method `InferenceCountKKCondPoissonOneLik$clone()` (current description `57` chars).
- [ ] TODO #347: Method `InferenceCountKKCondPoissonOneLik$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #348: Method `InferenceCountKKCondPoissonOneLik$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #349: Method `InferenceCountKKCondPoissonOneLik$compute_estimate_with_bootstrap_weights()` (current description `58` chars).
- [ ] TODO #350: Method `InferenceCountKKCondPoissonOneLik$compute_estimate()` (current description `130` chars).
- [ ] TODO #351: Method `InferenceCountKKCondPoissonOneLik$compute_gradient_confidence_interval()` (current description `56` chars).
- [ ] TODO #352: Method `InferenceCountKKCondPoissonOneLik$compute_gradient_two_sided_pval()` (current description `44` chars).
- [ ] TODO #353: Method `InferenceCountKKCondPoissonOneLik$compute_lik_ratio_confidence_interval()` (current description `64` chars).
- [ ] TODO #354: Method `InferenceCountKKCondPoissonOneLik$compute_lik_ratio_two_sided_pval()` (current description `52` chars).
- [ ] TODO #355: Method `InferenceCountKKCondPoissonOneLik$compute_score_confidence_interval()` (current description `53` chars).
- [ ] TODO #356: Method `InferenceCountKKCondPoissonOneLik$compute_score_two_sided_pval()` (current description `41` chars).
- [ ] TODO #357: Method `InferenceCountKKCondPoissonOneLik$compute_wald_confidence_interval()` (current description `54` chars).
- [ ] TODO #358: Method `InferenceCountKKCondPoissonOneLik$compute_wald_two_sided_pval()` (current description `34` chars).
- [ ] TODO #359: Method `InferenceCountKKCondPoissonOneLik$new()` (current description `120` chars).
- [ ] TODO #360: Method `InferenceCountKKCondPoissonOneLik$supports_lik_ratio_param_bootstrap()` (current description `59` chars).
- [ ] TODO #361: Topic `InferenceCountKKCondPoissonOneLik` (current description `139` chars).

### `InferenceCountKKGLMM.Rd`

- [ ] TODO #362: Method `InferenceCountKKGLMM$clone()` (current description `57` chars).
- [ ] TODO #363: Method `InferenceCountKKGLMM$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #364: Method `InferenceCountKKGLMM$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #365: Method `InferenceCountKKGLMM$compute_estimate_with_bootstrap_weights()` (current description `55` chars).
- [ ] TODO #366: Method `InferenceCountKKGLMM$compute_estimate()` (current description `58` chars).
- [ ] TODO #367: Method `InferenceCountKKGLMM$compute_lik_ratio_confidence_interval()` (current description `63` chars).
- [ ] TODO #368: Method `InferenceCountKKGLMM$compute_lik_ratio_two_sided_pval()` (current description `61` chars).
- [ ] TODO #369: Method `InferenceCountKKGLMM$compute_wald_confidence_interval()` (current description `54` chars).
- [ ] TODO #370: Method `InferenceCountKKGLMM$compute_wald_two_sided_pval()` (current description `34` chars).
- [ ] TODO #371: Method `InferenceCountKKGLMM$new()` (current description `46` chars).
- [ ] TODO #372: Topic `InferenceCountKKGLMM` (current description `348` chars).

### `InferenceCountKKHurdlePoissonOneLik.Rd`

- [ ] TODO #373: Method `InferenceCountKKHurdlePoissonOneLik$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #374: Method `InferenceCountKKHurdlePoissonOneLik$clone()` (current description `57` chars).
- [ ] TODO #375: Method `InferenceCountKKHurdlePoissonOneLik$compute_asymp_confidence_interval()` (current description `194` chars).
- [ ] TODO #376: Method `InferenceCountKKHurdlePoissonOneLik$compute_asymp_two_sided_pval()` (current description `192` chars).
- [ ] TODO #377: Method `InferenceCountKKHurdlePoissonOneLik$compute_estimate_with_bootstrap_weights()` (current description `53` chars).
- [ ] TODO #378: Method `InferenceCountKKHurdlePoissonOneLik$compute_estimate()` (current description `124` chars).
- [ ] TODO #379: Method `InferenceCountKKHurdlePoissonOneLik$compute_gradient_confidence_interval()` (current description `60` chars).
- [ ] TODO #380: Method `InferenceCountKKHurdlePoissonOneLik$compute_gradient_two_sided_pval()` (current description `48` chars).
- [ ] TODO #381: Method `InferenceCountKKHurdlePoissonOneLik$compute_lik_ratio_confidence_interval()` (current description `68` chars).
- [ ] TODO #382: Method `InferenceCountKKHurdlePoissonOneLik$compute_lik_ratio_two_sided_pval()` (current description `56` chars).
- [ ] TODO #383: Method `InferenceCountKKHurdlePoissonOneLik$compute_score_confidence_interval()` (current description `57` chars).
- [ ] TODO #384: Method `InferenceCountKKHurdlePoissonOneLik$compute_score_two_sided_pval()` (current description `45` chars).
- [ ] TODO #385: Method `InferenceCountKKHurdlePoissonOneLik$compute_wald_confidence_interval()` (current description `173` chars).
- [ ] TODO #386: Method `InferenceCountKKHurdlePoissonOneLik$compute_wald_two_sided_pval()` (current description `27` chars).
- [ ] TODO #387: Method `InferenceCountKKHurdlePoissonOneLik$new()` (current description `196` chars).
- [ ] TODO #388: Method `InferenceCountKKHurdlePoissonOneLik$supports_lik_ratio_param_bootstrap()` (current description `59` chars).
- [ ] TODO #389: Topic `InferenceCountKKHurdlePoissonOneLik` (current description `135` chars).

### `InferenceCountNegBin.Rd`

- [ ] TODO #390: Method `InferenceCountNegBin$clone()` (current description `57` chars).
- [ ] TODO #391: Method `InferenceCountNegBin$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #392: Method `InferenceCountNegBin$compute_jackknife_bias_estimate()` (current description `158` chars).
- [ ] TODO #393: Method `InferenceCountNegBin$compute_jackknife_estimate()` (current description `52` chars).
- [ ] TODO #394: Method `InferenceCountNegBin$compute_jackknife_std_error()` (current description `159` chars).
- [ ] TODO #395: Method `InferenceCountNegBin$compute_jackknife_wald_confidence_interval()` (current description `63` chars).
- [ ] TODO #396: Method `InferenceCountNegBin$compute_jackknife_wald_two_sided_pval()` (current description `62` chars).
- [ ] TODO #397: Method `InferenceCountNegBin$new()` (current description `59` chars).
- [ ] TODO #398: Topic `InferenceCountNegBin` (current description `200` chars).

### `InferenceCountPoisson.Rd`

- [ ] TODO #399: Method `InferenceCountPoisson$clone()` (current description `57` chars).
- [ ] TODO #400: Method `InferenceCountPoisson$compute_asymp_confidence_interval()` (current description `69` chars).
- [ ] TODO #401: Method `InferenceCountPoisson$compute_asymp_two_sided_pval()` (current description `67` chars).
- [ ] TODO #402: Method `InferenceCountPoisson$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #403: Method `InferenceCountPoisson$compute_gradient_confidence_interval()` (current description `60` chars).
- [ ] TODO #404: Method `InferenceCountPoisson$compute_gradient_two_sided_pval()` (current description `58` chars).
- [ ] TODO #405: Method `InferenceCountPoisson$compute_lik_ratio_bootstrap_confidence_interval()` (current description `58` chars).
- [ ] TODO #406: Method `InferenceCountPoisson$compute_lik_ratio_bootstrap_two_sided_pval()` (current description `63` chars).
- [ ] TODO #407: Method `InferenceCountPoisson$compute_lik_ratio_confidence_interval()` (current description `68` chars).
- [ ] TODO #408: Method `InferenceCountPoisson$compute_lik_ratio_two_sided_pval()` (current description `66` chars).
- [ ] TODO #409: Method `InferenceCountPoisson$compute_score_confidence_interval()` (current description `57` chars).
- [ ] TODO #410: Method `InferenceCountPoisson$compute_score_two_sided_pval()` (current description `55` chars).
- [ ] TODO #411: Method `InferenceCountPoisson$compute_wald_confidence_interval()` (current description `56` chars).
- [ ] TODO #412: Method `InferenceCountPoisson$compute_wald_two_sided_pval()` (current description `54` chars).
- [ ] TODO #413: Method `InferenceCountPoisson$new()` (current description `49` chars).
- [ ] TODO #414: Topic `InferenceCountPoisson` (current description `189` chars).

### `InferenceCountPoissonKKGEE.Rd`

- [ ] TODO #415: Method `InferenceCountPoissonKKGEE$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #416: Method `InferenceCountPoissonKKGEE$clone()` (current description `57` chars).
- [ ] TODO #417: Method `InferenceCountPoissonKKGEE$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #418: Method `InferenceCountPoissonKKGEE$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #419: Method `InferenceCountPoissonKKGEE$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #420: Method `InferenceCountPoissonKKGEE$compute_estimate()` (current description `170` chars).
- [ ] TODO #421: Method `InferenceCountPoissonKKGEE$new()` (current description `128` chars).
- [ ] TODO #422: Topic `InferenceCountPoissonKKGEE` (current description `297` chars).

### `InferenceCountQuasiPoisson.Rd`

- [ ] TODO #423: Method `InferenceCountQuasiPoisson$clone()` (current description `57` chars).
- [ ] TODO #424: Method `InferenceCountQuasiPoisson$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #425: Method `InferenceCountQuasiPoisson$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #426: Method `InferenceCountQuasiPoisson$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #427: Method `InferenceCountQuasiPoisson$compute_estimate()` (current description `58` chars).
- [ ] TODO #428: Method `InferenceCountQuasiPoisson$new()` (current description `55` chars).
- [ ] TODO #429: Topic `InferenceCountQuasiPoisson` (current description `201` chars).

### `InferenceCountRobustPoisson.Rd`

- [ ] TODO #430: Method `InferenceCountRobustPoisson$clone()` (current description `57` chars).
- [ ] TODO #431: Method `InferenceCountRobustPoisson$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #432: Method `InferenceCountRobustPoisson$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #433: Method `InferenceCountRobustPoisson$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #434: Method `InferenceCountRobustPoisson$compute_estimate()` (current description `58` chars).
- [ ] TODO #435: Method `InferenceCountRobustPoisson$new()` (current description `56` chars).
- [ ] TODO #436: Topic `InferenceCountRobustPoisson` (current description `244` chars).

### `InferenceCountZeroInflatedNegBin.Rd`

- [ ] TODO #437: Method `InferenceCountZeroInflatedNegBin$clone()` (current description `57` chars).
- [ ] TODO #438: Method `InferenceCountZeroInflatedNegBin$new()` (current description `62` chars).
- [ ] TODO #439: Topic `InferenceCountZeroInflatedNegBin` (current description `228` chars).

### `InferenceCountZeroInflatedPoisson.Rd`

- [ ] TODO #440: Method `InferenceCountZeroInflatedPoisson$clone()` (current description `57` chars).
- [ ] TODO #441: Method `InferenceCountZeroInflatedPoisson$new()` (current description `52` chars).
- [ ] TODO #442: Topic `InferenceCountZeroInflatedPoisson` (current description `208` chars).

### `InferenceIncidBinomialIdentityRiskDiff.Rd`

- [ ] TODO #443: Method `InferenceIncidBinomialIdentityRiskDiff$clone()` (current description `57` chars).
- [ ] TODO #444: Method `InferenceIncidBinomialIdentityRiskDiff$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #445: Method `InferenceIncidBinomialIdentityRiskDiff$compute_lik_ratio_confidence_interval()` (current description `68` chars).
- [ ] TODO #446: Method `InferenceIncidBinomialIdentityRiskDiff$new()` (current description `69` chars).
- [ ] TODO #447: Topic `InferenceIncidBinomialIdentityRiskDiff` (current description `282` chars).

### `InferenceIncidCMH.Rd`

- [ ] TODO #448: Method `InferenceIncidCMH$clone()` (current description `57` chars).
- [ ] TODO #449: Method `InferenceIncidCMH$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #450: Method `InferenceIncidCMH$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #451: Method `InferenceIncidCMH$new()` (current description `149` chars).
- [ ] TODO #452: Topic `InferenceIncidCMH` (current description `1121` chars).

### `InferenceIncidenceExactZhang.Rd`

- [ ] TODO #453: Method `InferenceIncidenceExactZhang$clone()` (current description `57` chars).
- [ ] TODO #454: Method `InferenceIncidenceExactZhang$compute_estimate()` (current description `47` chars).
- [ ] TODO #455: Method `InferenceIncidenceExactZhang$compute_exact_confidence_interval()` (current description `43` chars).
- [ ] TODO #456: Method `InferenceIncidenceExactZhang$new()` (current description `43` chars).
- [ ] TODO #457: Topic `InferenceIncidenceExactZhang` (current description `63` chars).

### `InferenceIncidExactBinomial.Rd`

- [ ] TODO #458: Method `InferenceIncidExactBinomial$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #459: Method `InferenceIncidExactBinomial$clone()` (current description `57` chars).
- [ ] TODO #460: Method `InferenceIncidExactBinomial$compute_estimate()` (current description `66` chars).
- [ ] TODO #461: Method `InferenceIncidExactBinomial$new()` (current description `72` chars).
- [ ] TODO #462: Topic `InferenceIncidExactBinomial` (current description `332` chars).

### `InferenceIncidExactFisher.Rd`

- [ ] TODO #463: Method `InferenceIncidExactFisher$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #464: Method `InferenceIncidExactFisher$clone()` (current description `57` chars).
- [ ] TODO #465: Method `InferenceIncidExactFisher$compute_estimate()` (current description `66` chars).
- [ ] TODO #466: Method `InferenceIncidExactFisher$new()` (current description `57` chars).
- [ ] TODO #467: Topic `InferenceIncidExactFisher` (current description `65` chars).

### `InferenceIncidExtendedRobins.Rd`

- [ ] TODO #468: Method `InferenceIncidExtendedRobins$clone()` (current description `57` chars).
- [ ] TODO #469: Method `InferenceIncidExtendedRobins$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #470: Method `InferenceIncidExtendedRobins$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #471: Method `InferenceIncidExtendedRobins$new()` (current description `62` chars).
- [ ] TODO #472: Topic `InferenceIncidExtendedRobins` (current description `329` chars).

### `InferenceIncidGCompRiskDiff.Rd`

- [ ] TODO #473: Method `InferenceIncidGCompRiskDiff$clone()` (current description `57` chars).
- [ ] TODO #474: Topic `InferenceIncidGCompRiskDiff` (current description `339` chars).

### `InferenceIncidGCompRiskRatio.Rd`

- [ ] TODO #475: Method `InferenceIncidGCompRiskRatio$clone()` (current description `57` chars).
- [ ] TODO #476: Topic `InferenceIncidGCompRiskRatio` (current description `329` chars).

### `InferenceIncidKKCondLogitIVWC.Rd`

- [ ] TODO #477: Method `InferenceIncidKKCondLogitIVWC$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #478: Method `InferenceIncidKKCondLogitIVWC$clone()` (current description `57` chars).
- [ ] TODO #479: Method `InferenceIncidKKCondLogitIVWC$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #480: Method `InferenceIncidKKCondLogitIVWC$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #481: Method `InferenceIncidKKCondLogitIVWC$compute_estimate()` (current description `58` chars).
- [ ] TODO #482: Method `InferenceIncidKKCondLogitIVWC$new()` (current description `129` chars).
- [ ] TODO #483: Topic `InferenceIncidKKCondLogitIVWC` (current description `129` chars).

### `InferenceIncidKKCondLogitOneLik.Rd`

- [ ] TODO #484: Method `InferenceIncidKKCondLogitOneLik$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #485: Method `InferenceIncidKKCondLogitOneLik$clone()` (current description `57` chars).
- [ ] TODO #486: Method `InferenceIncidKKCondLogitOneLik$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #487: Method `InferenceIncidKKCondLogitOneLik$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #488: Method `InferenceIncidKKCondLogitOneLik$compute_estimate_with_bootstrap_weights()` (current description `59` chars).
- [ ] TODO #489: Method `InferenceIncidKKCondLogitOneLik$compute_estimate()` (current description `58` chars).
- [ ] TODO #490: Method `InferenceIncidKKCondLogitOneLik$new()` (current description `209` chars).
- [ ] TODO #491: Topic `InferenceIncidKKCondLogitOneLik` (current description `308` chars).

### `InferenceIncidKKCondLogitPlusGLMMIVWC.Rd`

- [ ] TODO #492: Method `InferenceIncidKKCondLogitPlusGLMMIVWC$clone()` (current description `57` chars).
- [ ] TODO #493: Topic `InferenceIncidKKCondLogitPlusGLMMIVWC` (current description `368` chars).

### `InferenceIncidKKCondLogitPlusGLMMOneLik.Rd`

- [ ] TODO #494: Method `InferenceIncidKKCondLogitPlusGLMMOneLik$clone()` (current description `57` chars).
- [ ] TODO #495: Method `InferenceIncidKKCondLogitPlusGLMMOneLik$new()` (current description `122` chars).
- [ ] TODO #496: Topic `InferenceIncidKKCondLogitPlusGLMMOneLik` (current description `383` chars).

### `InferenceIncidKKGCompRiskDiff.Rd`

- [ ] TODO #497: Method `InferenceIncidKKGCompRiskDiff$clone()` (current description `57` chars).
- [ ] TODO #498: Topic `InferenceIncidKKGCompRiskDiff` (current description `497` chars).

### `InferenceIncidKKGCompRiskRatio.Rd`

- [ ] TODO #499: Method `InferenceIncidKKGCompRiskRatio$clone()` (current description `57` chars).
- [ ] TODO #500: Topic `InferenceIncidKKGCompRiskRatio` (current description `486` chars).

### `InferenceIncidKKGEE.Rd`

- [ ] TODO #501: Method `InferenceIncidKKGEE$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #502: Method `InferenceIncidKKGEE$clone()` (current description `57` chars).
- [ ] TODO #503: Method `InferenceIncidKKGEE$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #504: Method `InferenceIncidKKGEE$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #505: Method `InferenceIncidKKGEE$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #506: Method `InferenceIncidKKGEE$compute_estimate()` (current description `173` chars).
- [ ] TODO #507: Method `InferenceIncidKKGEE$new()` (current description `129` chars).
- [ ] TODO #508: Topic `InferenceIncidKKGEE` (current description `274` chars).

### `InferenceIncidKKModifiedPoisson.Rd`

- [ ] TODO #509: Method `InferenceIncidKKModifiedPoisson$clone()` (current description `57` chars).
- [ ] TODO #510: Topic `InferenceIncidKKModifiedPoisson` (current description `382` chars).

### `InferenceIncidKKNewcombeRiskDiff.Rd`

- [ ] TODO #511: Method `InferenceIncidKKNewcombeRiskDiff$clone()` (current description `57` chars).
- [ ] TODO #512: Method `InferenceIncidKKNewcombeRiskDiff$compute_estimate()` (current description `58` chars).
- [ ] TODO #513: Method `InferenceIncidKKNewcombeRiskDiff$new()` (current description `118` chars).
- [ ] TODO #514: Topic `InferenceIncidKKNewcombeRiskDiff` (current description `586` chars).

### `InferenceIncidLogBinomial.Rd`

- [ ] TODO #515: Method `InferenceIncidLogBinomial$clone()` (current description `57` chars).
- [ ] TODO #516: Method `InferenceIncidLogBinomial$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #517: Method `InferenceIncidLogBinomial$compute_gradient_confidence_interval()` (current description `60` chars).
- [ ] TODO #518: Method `InferenceIncidLogBinomial$compute_score_confidence_interval()` (current description `57` chars).
- [ ] TODO #519: Method `InferenceIncidLogBinomial$new()` (current description `54` chars).
- [ ] TODO #520: Topic `InferenceIncidLogBinomial` (current description `207` chars).

### `InferenceIncidLogRegr.Rd`

- [ ] TODO #521: Method `InferenceIncidLogRegr$clone()` (current description `57` chars).
- [ ] TODO #522: Method `InferenceIncidLogRegr$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #523: Method `InferenceIncidLogRegr$new()` (current description `50` chars).
- [ ] TODO #524: Topic `InferenceIncidLogRegr` (current description `199` chars).

### `InferenceIncidMiettinenNurminenRiskDiff.Rd`

- [ ] TODO #525: Method `InferenceIncidMiettinenNurminenRiskDiff$clone()` (current description `57` chars).
- [ ] TODO #526: Method `InferenceIncidMiettinenNurminenRiskDiff$compute_asymp_confidence_interval()` (current description `66` chars).
- [ ] TODO #527: Method `InferenceIncidMiettinenNurminenRiskDiff$compute_asymp_two_sided_pval()` (current description `54` chars).
- [ ] TODO #528: Method `InferenceIncidMiettinenNurminenRiskDiff$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #529: Method `InferenceIncidMiettinenNurminenRiskDiff$compute_estimate()` (current description `47` chars).
- [ ] TODO #530: Method `InferenceIncidMiettinenNurminenRiskDiff$new()` (current description `70` chars).
- [ ] TODO #531: Topic `InferenceIncidMiettinenNurminenRiskDiff` (current description `598` chars).

### `InferenceIncidModifiedPoisson.Rd`

- [ ] TODO #532: Method `InferenceIncidModifiedPoisson$clone()` (current description `57` chars).
- [ ] TODO #533: Method `InferenceIncidModifiedPoisson$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #534: Method `InferenceIncidModifiedPoisson$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #535: Method `InferenceIncidModifiedPoisson$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #536: Method `InferenceIncidModifiedPoisson$compute_estimate()` (current description `58` chars).
- [ ] TODO #537: Method `InferenceIncidModifiedPoisson$new()` (current description `58` chars).
- [ ] TODO #538: Topic `InferenceIncidModifiedPoisson` (current description `316` chars).

### `InferenceIncidNewcombeRiskDiff.Rd`

- [ ] TODO #539: Method `InferenceIncidNewcombeRiskDiff$clone()` (current description `57` chars).
- [ ] TODO #540: Method `InferenceIncidNewcombeRiskDiff$compute_asymp_confidence_interval()` (current description `50` chars).
- [ ] TODO #541: Method `InferenceIncidNewcombeRiskDiff$compute_asymp_two_sided_pval()` (current description `64` chars).
- [ ] TODO #542: Method `InferenceIncidNewcombeRiskDiff$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #543: Method `InferenceIncidNewcombeRiskDiff$compute_estimate()` (current description `47` chars).
- [ ] TODO #544: Method `InferenceIncidNewcombeRiskDiff$new()` (current description `55` chars).
- [ ] TODO #545: Topic `InferenceIncidNewcombeRiskDiff` (current description `496` chars).

### `InferenceIncidProbitRegr.Rd`

- [ ] TODO #546: Method `InferenceIncidProbitRegr$clone()` (current description `57` chars).
- [ ] TODO #547: Method `InferenceIncidProbitRegr$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #548: Method `InferenceIncidProbitRegr$new()` (current description `48` chars).
- [ ] TODO #549: Topic `InferenceIncidProbitRegr` (current description `195` chars).

### `InferenceIncidRiskDiff.Rd`

- [ ] TODO #550: Method `InferenceIncidRiskDiff$clone()` (current description `57` chars).
- [ ] TODO #551: Method `InferenceIncidRiskDiff$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #552: Method `InferenceIncidRiskDiff$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #553: Method `InferenceIncidRiskDiff$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #554: Method `InferenceIncidRiskDiff$compute_estimate()` (current description `58` chars).
- [ ] TODO #555: Method `InferenceIncidRiskDiff$new()` (current description `46` chars).
- [ ] TODO #556: Topic `InferenceIncidRiskDiff` (current description `263` chars).

### `InferenceIncidWald.Rd`

- [ ] TODO #557: Method `InferenceIncidWald$clone()` (current description `57` chars).
- [ ] TODO #558: Method `InferenceIncidWald$new()` (current description `166` chars).
- [ ] TODO #559: Topic `InferenceIncidWald` (current description `211` chars).

### `InferenceOrdinalAdjCatLogitRegr.Rd`

- [ ] TODO #560: Method `InferenceOrdinalAdjCatLogitRegr$clone()` (current description `57` chars).
- [ ] TODO #561: Method `InferenceOrdinalAdjCatLogitRegr$compute_estimate_with_bootstrap_weights()` (current description `55` chars).
- [ ] TODO #562: Method `InferenceOrdinalAdjCatLogitRegr$new()` (current description `55` chars).
- [ ] TODO #563: Topic `InferenceOrdinalAdjCatLogitRegr` (current description `217` chars).

### `InferenceOrdinalCauchitRegr.Rd`

- [ ] TODO #564: Method `InferenceOrdinalCauchitRegr$clone()` (current description `57` chars).
- [ ] TODO #565: Method `InferenceOrdinalCauchitRegr$compute_estimate_with_bootstrap_weights()` (current description `55` chars).
- [ ] TODO #566: Method `InferenceOrdinalCauchitRegr$new()` (current description `46` chars).
- [ ] TODO #567: Topic `InferenceOrdinalCauchitRegr` (current description `193` chars).

### `InferenceOrdinalCloglogRegr.Rd`

- [ ] TODO #568: Method `InferenceOrdinalCloglogRegr$clone()` (current description `57` chars).
- [ ] TODO #569: Method `InferenceOrdinalCloglogRegr$compute_estimate_with_bootstrap_weights()` (current description `55` chars).
- [ ] TODO #570: Method `InferenceOrdinalCloglogRegr$new()` (current description `49` chars).
- [ ] TODO #571: Topic `InferenceOrdinalCloglogRegr` (current description `193` chars).

### `InferenceOrdinalContRatioRegr.Rd`

- [ ] TODO #572: Method `InferenceOrdinalContRatioRegr$clone()` (current description `57` chars).
- [ ] TODO #573: Method `InferenceOrdinalContRatioRegr$compute_estimate_with_bootstrap_weights()` (current description `58` chars).
- [ ] TODO #574: Method `InferenceOrdinalContRatioRegr$new()` (current description `49` chars).
- [ ] TODO #575: Topic `InferenceOrdinalContRatioRegr` (current description `182` chars).

### `InferenceOrdinalGCompMeanDiff.Rd`

- [ ] TODO #576: Method `InferenceOrdinalGCompMeanDiff$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #577: Method `InferenceOrdinalGCompMeanDiff$clone()` (current description `57` chars).
- [ ] TODO #578: Method `InferenceOrdinalGCompMeanDiff$compute_asymp_confidence_interval()` (current description `72` chars).
- [ ] TODO #579: Method `InferenceOrdinalGCompMeanDiff$compute_asymp_two_sided_pval()` (current description `65` chars).
- [ ] TODO #580: Method `InferenceOrdinalGCompMeanDiff$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #581: Method `InferenceOrdinalGCompMeanDiff$compute_estimate()` (current description `80` chars).
- [ ] TODO #582: Method `InferenceOrdinalGCompMeanDiff$compute_wald_confidence_interval()` (current description `67` chars).
- [ ] TODO #583: Method `InferenceOrdinalGCompMeanDiff$compute_wald_two_sided_pval()` (current description `65` chars).
- [ ] TODO #584: Method `InferenceOrdinalGCompMeanDiff$new()` (current description `55` chars).
- [ ] TODO #585: Topic `InferenceOrdinalGCompMeanDiff` (current description `352` chars).

### `InferenceOrdinalJonckheereTerpstraTest.Rd`

- [ ] TODO #586: Method `InferenceOrdinalJonckheereTerpstraTest$clone()` (current description `57` chars).
- [ ] TODO #587: Method `InferenceOrdinalJonckheereTerpstraTest$compute_asymp_confidence_interval()` (current description `65` chars).
- [ ] TODO #588: Method `InferenceOrdinalJonckheereTerpstraTest$compute_asymp_two_sided_pval()` (current description `63` chars).
- [ ] TODO #589: Method `InferenceOrdinalJonckheereTerpstraTest$compute_estimate_with_bootstrap_weights()` (current description `82` chars).
- [ ] TODO #590: Method `InferenceOrdinalJonckheereTerpstraTest$compute_estimate()` (current description `64` chars).
- [ ] TODO #591: Method `InferenceOrdinalJonckheereTerpstraTest$compute_exact_two_sided_pval_for_treatment_effect()` (current description `36` chars).
- [ ] TODO #592: Method `InferenceOrdinalJonckheereTerpstraTest$new()` (current description `30` chars).
- [ ] TODO #593: Topic `InferenceOrdinalJonckheereTerpstraTest` (current description `310` chars).

### `InferenceOrdinalKKCLMM.Rd`

- [ ] TODO #594: Method `InferenceOrdinalKKCLMM$clone()` (current description `57` chars).
- [ ] TODO #595: Method `InferenceOrdinalKKCLMM$new()` (current description `98` chars).
- [ ] TODO #596: Topic `InferenceOrdinalKKCLMM` (current description `186` chars).

### `InferenceOrdinalKKCLMMCauchit.Rd`

- [ ] TODO #597: Method `InferenceOrdinalKKCLMMCauchit$clone()` (current description `57` chars).
- [ ] TODO #598: Method `InferenceOrdinalKKCLMMCauchit$new()` (current description `100` chars).
- [ ] TODO #599: Topic `InferenceOrdinalKKCLMMCauchit` (current description `85` chars).

### `InferenceOrdinalKKCLMMCloglog.Rd`

- [ ] TODO #600: Method `InferenceOrdinalKKCLMMCloglog$clone()` (current description `57` chars).
- [ ] TODO #601: Method `InferenceOrdinalKKCLMMCloglog$new()` (current description `100` chars).
- [ ] TODO #602: Topic `InferenceOrdinalKKCLMMCloglog` (current description `113` chars).

### `InferenceOrdinalKKCLMMProbit.Rd`

- [ ] TODO #603: Method `InferenceOrdinalKKCLMMProbit$clone()` (current description `57` chars).
- [ ] TODO #604: Method `InferenceOrdinalKKCLMMProbit$new()` (current description `99` chars).
- [ ] TODO #605: Topic `InferenceOrdinalKKCLMMProbit` (current description `83` chars).

### `InferenceOrdinalKKCondAdjCatLogitRegr.Rd`

- [ ] TODO #606: Method `InferenceOrdinalKKCondAdjCatLogitRegr$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #607: Method `InferenceOrdinalKKCondAdjCatLogitRegr$clone()` (current description `57` chars).
- [ ] TODO #608: Method `InferenceOrdinalKKCondAdjCatLogitRegr$compute_asymp_confidence_interval()` (current description `43` chars).
- [ ] TODO #609: Method `InferenceOrdinalKKCondAdjCatLogitRegr$compute_asymp_two_sided_pval()` (current description `128` chars).
- [ ] TODO #610: Method `InferenceOrdinalKKCondAdjCatLogitRegr$compute_estimate_with_bootstrap_weights()` (current description `53` chars).
- [ ] TODO #611: Method `InferenceOrdinalKKCondAdjCatLogitRegr$compute_estimate()` (current description `38` chars).
- [ ] TODO #612: Method `InferenceOrdinalKKCondAdjCatLogitRegr$new()` (current description `131` chars).
- [ ] TODO #613: Topic `InferenceOrdinalKKCondAdjCatLogitRegr` (current description `260` chars).

### `InferenceOrdinalKKGEE.Rd`

- [ ] TODO #614: Method `InferenceOrdinalKKGEE$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #615: Method `InferenceOrdinalKKGEE$clone()` (current description `57` chars).
- [ ] TODO #616: Method `InferenceOrdinalKKGEE$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #617: Method `InferenceOrdinalKKGEE$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #618: Method `InferenceOrdinalKKGEE$compute_estimate_with_bootstrap_weights()` (current description `54` chars).
- [ ] TODO #619: Method `InferenceOrdinalKKGEE$compute_estimate()` (current description `166` chars).
- [ ] TODO #620: Method `InferenceOrdinalKKGEE$new()` (current description `127` chars).
- [ ] TODO #621: Topic `InferenceOrdinalKKGEE` (current description `264` chars).

### `InferenceOrdinalKKGLMM.Rd`

- [ ] TODO #622: Method `InferenceOrdinalKKGLMM$clone()` (current description `57` chars).
- [ ] TODO #623: Method `InferenceOrdinalKKGLMM$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #624: Method `InferenceOrdinalKKGLMM$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #625: Method `InferenceOrdinalKKGLMM$compute_estimate_with_bootstrap_weights()` (current description `55` chars).
- [ ] TODO #626: Method `InferenceOrdinalKKGLMM$compute_estimate()` (current description `58` chars).
- [ ] TODO #627: Method `InferenceOrdinalKKGLMM$new()` (current description `131` chars).
- [ ] TODO #628: Topic `InferenceOrdinalKKGLMM` (current description `430` chars).

### `InferenceOrdinalOrderedProbitRegr.Rd`

- [ ] TODO #629: Method `InferenceOrdinalOrderedProbitRegr$clone()` (current description `57` chars).
- [ ] TODO #630: Method `InferenceOrdinalOrderedProbitRegr$compute_estimate_with_bootstrap_weights()` (current description `55` chars).
- [ ] TODO #631: Method `InferenceOrdinalOrderedProbitRegr$new()` (current description `46` chars).
- [ ] TODO #632: Topic `InferenceOrdinalOrderedProbitRegr` (current description `199` chars).

### `InferenceOrdinalPairedSignTest.Rd`

- [ ] TODO #633: Method `InferenceOrdinalPairedSignTest$approximate_bootstrap_distribution_beta_hat_T()` (current description `182` chars).
- [ ] TODO #634: Method `InferenceOrdinalPairedSignTest$approximate_jackknife_distribution_beta_hat_T()` (current description `162` chars).
- [ ] TODO #635: Method `InferenceOrdinalPairedSignTest$clone()` (current description `57` chars).
- [ ] TODO #636: Method `InferenceOrdinalPairedSignTest$compute_asymp_confidence_interval()` (current description `62` chars).
- [ ] TODO #637: Method `InferenceOrdinalPairedSignTest$compute_asymp_two_sided_pval()` (current description `39` chars).
- [ ] TODO #638: Method `InferenceOrdinalPairedSignTest$compute_estimate_with_bootstrap_weights()` (current description `51` chars).
- [ ] TODO #639: Method `InferenceOrdinalPairedSignTest$compute_estimate()` (current description `73` chars).
- [ ] TODO #640: Method `InferenceOrdinalPairedSignTest$new()` (current description `122` chars).
- [ ] TODO #641: Topic `InferenceOrdinalPairedSignTest` (current description `320` chars).

### `InferenceOrdinalPartialProportionalOddsRegr.Rd`

- [ ] TODO #642: Method `InferenceOrdinalPartialProportionalOddsRegr$benchmark_asymp_two_sided_pval_breakdown()` (current description `134` chars).
- [ ] TODO #643: Method `InferenceOrdinalPartialProportionalOddsRegr$clone()` (current description `57` chars).
- [ ] TODO #644: Method `InferenceOrdinalPartialProportionalOddsRegr$compute_asymp_confidence_interval()` (current description `143` chars).
- [ ] TODO #645: Method `InferenceOrdinalPartialProportionalOddsRegr$compute_asymp_two_sided_pval()` (current description `140` chars).
- [ ] TODO #646: Method `InferenceOrdinalPartialProportionalOddsRegr$compute_estimate_with_bootstrap_weights()` (current description `59` chars).
- [ ] TODO #647: Method `InferenceOrdinalPartialProportionalOddsRegr$compute_estimate()` (current description `48` chars).
- [ ] TODO #648: Method `InferenceOrdinalPartialProportionalOddsRegr$compute_wald_confidence_interval()` (current description `143` chars).
- [ ] TODO #649: Method `InferenceOrdinalPartialProportionalOddsRegr$compute_wald_two_sided_pval()` (current description `140` chars).
- [ ] TODO #650: Method `InferenceOrdinalPartialProportionalOddsRegr$new()` (current description `40` chars).
- [ ] TODO #651: Topic `InferenceOrdinalPartialProportionalOddsRegr` (current description `306` chars).

### `InferenceOrdinalPropOddsRegr.Rd`

- [ ] TODO #652: Method `InferenceOrdinalPropOddsRegr$clone()` (current description `57` chars).
- [ ] TODO #653: Method `InferenceOrdinalPropOddsRegr$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #654: Method `InferenceOrdinalPropOddsRegr$new()` (current description `48` chars).
- [ ] TODO #655: Topic `InferenceOrdinalPropOddsRegr` (current description `204` chars).

### `InferenceOrdinalRidit.Rd`

- [ ] TODO #656: Method `InferenceOrdinalRidit$clone()` (current description `57` chars).
- [ ] TODO #657: Method `InferenceOrdinalRidit$compute_asymp_confidence_interval()` (current description `69` chars).
- [ ] TODO #658: Method `InferenceOrdinalRidit$compute_asymp_two_sided_pval()` (current description `67` chars).
- [ ] TODO #659: Method `InferenceOrdinalRidit$compute_estimate_with_bootstrap_weights()` (current description `45` chars).
- [ ] TODO #660: Method `InferenceOrdinalRidit$compute_estimate()` (current description `58` chars).
- [ ] TODO #661: Method `InferenceOrdinalRidit$get_mean_ridit_treatment()` (current description `47` chars).
- [ ] TODO #662: Method `InferenceOrdinalRidit$get_ridit_scores()` (current description `42` chars).
- [ ] TODO #663: Method `InferenceOrdinalRidit$new()` (current description `45` chars).
- [ ] TODO #664: Topic `InferenceOrdinalRidit` (current description `358` chars).

### `InferenceOrdinalStereotypeLogitRegr.Rd`

- [ ] TODO #665: Method `InferenceOrdinalStereotypeLogitRegr$clone()` (current description `57` chars).
- [ ] TODO #666: Method `InferenceOrdinalStereotypeLogitRegr$compute_estimate_with_bootstrap_weights()` (current description `56` chars).
- [ ] TODO #667: Method `InferenceOrdinalStereotypeLogitRegr$new()` (current description `47` chars).
- [ ] TODO #668: Topic `InferenceOrdinalStereotypeLogitRegr` (current description `202` chars).

### `InferencePropBetaRegr.Rd`

- [ ] TODO #669: Method `InferencePropBetaRegr$clone()` (current description `57` chars).
- [ ] TODO #670: Method `InferencePropBetaRegr$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #671: Method `InferencePropBetaRegr$compute_estimate()` (current description `58` chars).
- [ ] TODO #672: Method `InferencePropBetaRegr$new()` (current description `46` chars).
- [ ] TODO #673: Topic `InferencePropBetaRegr` (current description `208` chars).

### `InferencePropFractionalLogit.Rd`

- [ ] TODO #674: Method `InferencePropFractionalLogit$clone()` (current description `57` chars).
- [ ] TODO #675: Method `InferencePropFractionalLogit$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #676: Method `InferencePropFractionalLogit$compute_estimate()` (current description `58` chars).
- [ ] TODO #677: Method `InferencePropFractionalLogit$new()` (current description `47` chars).
- [ ] TODO #678: Topic `InferencePropFractionalLogit` (current description `241` chars).

### `InferencePropGCompMeanDiff.Rd`

- [ ] TODO #679: Method `InferencePropGCompMeanDiff$clone()` (current description `57` chars).
- [ ] TODO #680: Topic `InferencePropGCompMeanDiff` (current description `362` chars).

### `InferencePropKKGEE.Rd`

- [ ] TODO #681: Method `InferencePropKKGEE$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #682: Method `InferencePropKKGEE$clone()` (current description `57` chars).
- [ ] TODO #683: Method `InferencePropKKGEE$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #684: Method `InferencePropKKGEE$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #685: Method `InferencePropKKGEE$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #686: Method `InferencePropKKGEE$compute_estimate()` (current description `171` chars).
- [ ] TODO #687: Method `InferencePropKKGEE$new()` (current description `131` chars).
- [ ] TODO #688: Topic `InferencePropKKGEE` (current description `300` chars).

### `InferencePropKKGLMM.Rd`

- [ ] TODO #689: Method `InferencePropKKGLMM$clone()` (current description `57` chars).
- [ ] TODO #690: Method `InferencePropKKGLMM$new()` (current description `95` chars).
- [ ] TODO #691: Topic `InferencePropKKGLMM` (current description `194` chars).

### `InferencePropKKQuantileRegrIVWC.Rd`

- [ ] TODO #692: Method `InferencePropKKQuantileRegrIVWC$clone()` (current description `57` chars).
- [ ] TODO #693: Method `InferencePropKKQuantileRegrIVWC$new()` (current description `101` chars).
- [ ] TODO #694: Topic `InferencePropKKQuantileRegrIVWC` (current description `1458` chars).

### `InferencePropKKQuantileRegrOneLik.Rd`

- [ ] TODO #695: Method `InferencePropKKQuantileRegrOneLik$clone()` (current description `57` chars).
- [ ] TODO #696: Method `InferencePropKKQuantileRegrOneLik$compute_estimate()` (current description `57` chars).
- [ ] TODO #697: Method `InferencePropKKQuantileRegrOneLik$new()` (current description `185` chars).
- [ ] TODO #698: Topic `InferencePropKKQuantileRegrOneLik` (current description `542` chars).

### `InferencePropQuantileRegr.Rd`

- [ ] TODO #699: Method `InferencePropQuantileRegr$clone()` (current description `57` chars).
- [ ] TODO #700: Method `InferencePropQuantileRegr$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #701: Method `InferencePropQuantileRegr$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #702: Method `InferencePropQuantileRegr$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #703: Method `InferencePropQuantileRegr$compute_estimate()` (current description `66` chars).
- [ ] TODO #704: Method `InferencePropQuantileRegr$new()` (current description `72` chars).
- [ ] TODO #705: Topic `InferencePropQuantileRegr` (current description `834` chars).

### `InferencePropZeroOneInflatedBetaRegr.Rd`

- [ ] TODO #706: Method `InferencePropZeroOneInflatedBetaRegr$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #707: Method `InferencePropZeroOneInflatedBetaRegr$clone()` (current description `57` chars).
- [ ] TODO #708: Method `InferencePropZeroOneInflatedBetaRegr$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #709: Method `InferencePropZeroOneInflatedBetaRegr$compute_estimate()` (current description `58` chars).
- [ ] TODO #710: Method `InferencePropZeroOneInflatedBetaRegr$new()` (current description `64` chars).
- [ ] TODO #711: Topic `InferencePropZeroOneInflatedBetaRegr` (current description `655` chars).

### `InferenceSuite.Rd`

- [ ] TODO #712: Method `InferenceSuite$clone()` (current description `57` chars).
- [ ] TODO #713: Method `InferenceSuite$new()` (current description `164` chars).
- [ ] TODO #714: Topic `InferenceSuite` (current description `706` chars).

### `InferenceSurvivalCoxPHRegr.Rd`

- [ ] TODO #715: Method `InferenceSurvivalCoxPHRegr$clone()` (current description `57` chars).
- [ ] TODO #716: Method `InferenceSurvivalCoxPHRegr$compute_estimate_with_bootstrap_weights()` (current description `46` chars).
- [ ] TODO #717: Method `InferenceSurvivalCoxPHRegr$compute_estimate()` (current description `53` chars).
- [ ] TODO #718: Method `InferenceSurvivalCoxPHRegr$new()` (current description `37` chars).
- [ ] TODO #719: Topic `InferenceSurvivalCoxPHRegr` (current description `220` chars).

### `InferenceSurvivalDepCensTransformRegr.Rd`

- [ ] TODO #720: Method `InferenceSurvivalDepCensTransformRegr$approximate_randomization_distribution_beta_hat_T()` (current description `69` chars).
- [ ] TODO #721: Method `InferenceSurvivalDepCensTransformRegr$clone()` (current description `57` chars).
- [ ] TODO #722: Method `InferenceSurvivalDepCensTransformRegr$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #723: Method `InferenceSurvivalDepCensTransformRegr$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #724: Method `InferenceSurvivalDepCensTransformRegr$compute_bootstrap_confidence_interval_basic()` (current description `62` chars).
- [ ] TODO #725: Method `InferenceSurvivalDepCensTransformRegr$compute_bootstrap_confidence_interval_bca()` (current description `60` chars).
- [ ] TODO #726: Method `InferenceSurvivalDepCensTransformRegr$compute_bootstrap_confidence_interval_studentized()` (current description `68` chars).
- [ ] TODO #727: Method `InferenceSurvivalDepCensTransformRegr$compute_bootstrap_confidence_interval()` (current description `56` chars).
- [ ] TODO #728: Method `InferenceSurvivalDepCensTransformRegr$compute_estimate_with_bootstrap_weights()` (current description `58` chars).
- [ ] TODO #729: Method `InferenceSurvivalDepCensTransformRegr$compute_estimate()` (current description `58` chars).
- [ ] TODO #730: Method `InferenceSurvivalDepCensTransformRegr$compute_jackknife_bias_estimate()` (current description `62` chars).
- [ ] TODO #731: Method `InferenceSurvivalDepCensTransformRegr$compute_jackknife_estimate()` (current description `114` chars).
- [ ] TODO #732: Method `InferenceSurvivalDepCensTransformRegr$compute_jackknife_std_error()` (current description `63` chars).
- [ ] TODO #733: Method `InferenceSurvivalDepCensTransformRegr$compute_jackknife_wald_confidence_interval()` (current description `73` chars).
- [ ] TODO #734: Method `InferenceSurvivalDepCensTransformRegr$compute_jackknife_wald_two_sided_pval()` (current description `71` chars).
- [ ] TODO #735: Method `InferenceSurvivalDepCensTransformRegr$compute_lik_ratio_confidence_interval()` (current description `58` chars).
- [ ] TODO #736: Method `InferenceSurvivalDepCensTransformRegr$compute_rand_confidence_interval()` (current description `72` chars).
- [ ] TODO #737: Method `InferenceSurvivalDepCensTransformRegr$compute_rand_two_sided_pval()` (current description `122` chars).
- [ ] TODO #738: Method `InferenceSurvivalDepCensTransformRegr$compute_score_two_sided_pval()` (current description `89` chars).
- [ ] TODO #739: Method `InferenceSurvivalDepCensTransformRegr$new()` (current description `65` chars).
- [ ] TODO #740: Topic `InferenceSurvivalDepCensTransformRegr` (current description `240` chars).

### `InferenceSurvivalGehanWilcox.Rd`

- [ ] TODO #741: Method `InferenceSurvivalGehanWilcox$clone()` (current description `57` chars).
- [ ] TODO #742: Method `InferenceSurvivalGehanWilcox$compute_asymp_confidence_interval()` (current description `163` chars).
- [ ] TODO #743: Method `InferenceSurvivalGehanWilcox$compute_asymp_two_sided_pval()` (current description `186` chars).
- [ ] TODO #744: Method `InferenceSurvivalGehanWilcox$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #745: Method `InferenceSurvivalGehanWilcox$compute_estimate()` (current description `156` chars).
- [ ] TODO #746: Method `InferenceSurvivalGehanWilcox$compute_rand_confidence_interval()` (current description `154` chars).
- [ ] TODO #747: Method `InferenceSurvivalGehanWilcox$new()` (current description `99` chars).
- [ ] TODO #748: Topic `InferenceSurvivalGehanWilcox` (current description `931` chars).

### `InferenceSurvivalKKClaytonCopulaIVWC.Rd`

- [ ] TODO #749: Method `InferenceSurvivalKKClaytonCopulaIVWC$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #750: Method `InferenceSurvivalKKClaytonCopulaIVWC$clone()` (current description `57` chars).
- [ ] TODO #751: Method `InferenceSurvivalKKClaytonCopulaIVWC$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #752: Method `InferenceSurvivalKKClaytonCopulaIVWC$compute_asymp_two_sided_pval()` (current description `160` chars).
- [ ] TODO #753: Method `InferenceSurvivalKKClaytonCopulaIVWC$compute_estimate_with_bootstrap_weights()` (current description `59` chars).
- [ ] TODO #754: Method `InferenceSurvivalKKClaytonCopulaIVWC$compute_estimate()` (current description `65` chars).
- [ ] TODO #755: Method `InferenceSurvivalKKClaytonCopulaIVWC$duplicate()` (current description `57` chars).
- [ ] TODO #756: Method `InferenceSurvivalKKClaytonCopulaIVWC$new()` (current description `111` chars).
- [ ] TODO #757: Topic `InferenceSurvivalKKClaytonCopulaIVWC` (current description `453` chars).

### `InferenceSurvivalKKClaytonCopulaOneLik.Rd`

- [ ] TODO #758: Method `InferenceSurvivalKKClaytonCopulaOneLik$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #759: Method `InferenceSurvivalKKClaytonCopulaOneLik$clone()` (current description `57` chars).
- [ ] TODO #760: Method `InferenceSurvivalKKClaytonCopulaOneLik$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #761: Method `InferenceSurvivalKKClaytonCopulaOneLik$compute_asymp_two_sided_pval()` (current description `122` chars).
- [ ] TODO #762: Method `InferenceSurvivalKKClaytonCopulaOneLik$compute_estimate()` (current description `65` chars).
- [ ] TODO #763: Method `InferenceSurvivalKKClaytonCopulaOneLik$duplicate()` (current description `57` chars).
- [ ] TODO #764: Method `InferenceSurvivalKKClaytonCopulaOneLik$new()` (current description `106` chars).
- [ ] TODO #765: Topic `InferenceSurvivalKKClaytonCopulaOneLik` (current description `119` chars).

### `InferenceSurvivalKKLWACoxPHIVWC.Rd`

- [ ] TODO #766: Method `InferenceSurvivalKKLWACoxPHIVWC$clone()` (current description `57` chars).
- [ ] TODO #767: Method `InferenceSurvivalKKLWACoxPHIVWC$new()` (current description `81` chars).
- [ ] TODO #768: Topic `InferenceSurvivalKKLWACoxPHIVWC` (current description `292` chars).

### `InferenceSurvivalKKLWACoxPHOneLik.Rd`

- [ ] TODO #769: Method `InferenceSurvivalKKLWACoxPHOneLik$clone()` (current description `57` chars).
- [ ] TODO #770: Method `InferenceSurvivalKKLWACoxPHOneLik$new()` (current description `91` chars).
- [ ] TODO #771: Topic `InferenceSurvivalKKLWACoxPHOneLik` (current description `208` chars).

### `InferenceSurvivalKKRankRegrIVWC.Rd`

- [ ] TODO #772: Method `InferenceSurvivalKKRankRegrIVWC$clone()` (current description `57` chars).
- [ ] TODO #773: Method `InferenceSurvivalKKRankRegrIVWC$new()` (current description `127` chars).
- [ ] TODO #774: Topic `InferenceSurvivalKKRankRegrIVWC` (current description `265` chars).

### `InferenceSurvivalKKStratCoxPHIVWC.Rd`

- [ ] TODO #775: Method `InferenceSurvivalKKStratCoxPHIVWC$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #776: Method `InferenceSurvivalKKStratCoxPHIVWC$clone()` (current description `57` chars).
- [ ] TODO #777: Method `InferenceSurvivalKKStratCoxPHIVWC$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #778: Method `InferenceSurvivalKKStratCoxPHIVWC$compute_asymp_two_sided_pval()` (current description `125` chars).
- [ ] TODO #779: Method `InferenceSurvivalKKStratCoxPHIVWC$compute_estimate()` (current description `58` chars).
- [ ] TODO #780: Method `InferenceSurvivalKKStratCoxPHIVWC$new()` (current description `115` chars).
- [ ] TODO #781: Topic `InferenceSurvivalKKStratCoxPHIVWC` (current description `746` chars).

### `InferenceSurvivalKKStratCoxPHOneLik.Rd`

- [ ] TODO #782: Method `InferenceSurvivalKKStratCoxPHOneLik$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #783: Method `InferenceSurvivalKKStratCoxPHOneLik$clone()` (current description `57` chars).
- [ ] TODO #784: Method `InferenceSurvivalKKStratCoxPHOneLik$compute_asymp_confidence_interval()` (current description `68` chars).
- [ ] TODO #785: Method `InferenceSurvivalKKStratCoxPHOneLik$compute_asymp_two_sided_pval()` (current description `49` chars).
- [ ] TODO #786: Method `InferenceSurvivalKKStratCoxPHOneLik$compute_estimate_with_bootstrap_weights()` (current description `54` chars).
- [ ] TODO #787: Method `InferenceSurvivalKKStratCoxPHOneLik$compute_estimate()` (current description `70` chars).
- [ ] TODO #788: Method `InferenceSurvivalKKStratCoxPHOneLik$new()` (current description `109` chars).
- [ ] TODO #789: Topic `InferenceSurvivalKKStratCoxPHOneLik` (current description `137` chars).

### `InferenceSurvivalKKWeibullFrailtyIVWC.Rd`

- [ ] TODO #790: Method `InferenceSurvivalKKWeibullFrailtyIVWC$clone()` (current description `57` chars).
- [ ] TODO #791: Method `InferenceSurvivalKKWeibullFrailtyIVWC$new()` (current description `53` chars).
- [ ] TODO #792: Topic `InferenceSurvivalKKWeibullFrailtyIVWC` (current description `102` chars).

### `InferenceSurvivalKKWeibullFrailtyOneLik.Rd`

- [ ] TODO #793: Method `InferenceSurvivalKKWeibullFrailtyOneLik$clone()` (current description `57` chars).
- [ ] TODO #794: Method `InferenceSurvivalKKWeibullFrailtyOneLik$new()` (current description `63` chars).
- [ ] TODO #795: Topic `InferenceSurvivalKKWeibullFrailtyOneLik` (current description `121` chars).

### `InferenceSurvivalKKWeibullMarginal.Rd`

- [ ] TODO #796: Method `InferenceSurvivalKKWeibullMarginal$approximate_bootstrap_distribution_beta_hat_T()` (current description `68` chars).
- [ ] TODO #797: Method `InferenceSurvivalKKWeibullMarginal$clone()` (current description `57` chars).
- [ ] TODO #798: Method `InferenceSurvivalKKWeibullMarginal$compute_asymp_confidence_interval()` (current description `61` chars).
- [ ] TODO #799: Method `InferenceSurvivalKKWeibullMarginal$compute_asymp_two_sided_pval()` (current description `59` chars).
- [ ] TODO #800: Method `InferenceSurvivalKKWeibullMarginal$compute_estimate_with_bootstrap_weights()` (current description `67` chars).
- [ ] TODO #801: Method `InferenceSurvivalKKWeibullMarginal$compute_estimate()` (current description `68` chars).
- [ ] TODO #802: Method `InferenceSurvivalKKWeibullMarginal$duplicate()` (current description `57` chars).
- [ ] TODO #803: Method `InferenceSurvivalKKWeibullMarginal$new()` (current description `66` chars).
- [ ] TODO #804: Topic `InferenceSurvivalKKWeibullMarginal` (current description `1180` chars).

### `InferenceSurvivalKMDiff.Rd`

- [ ] TODO #805: Method `InferenceSurvivalKMDiff$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #806: Method `InferenceSurvivalKMDiff$clone()` (current description `57` chars).
- [ ] TODO #807: Method `InferenceSurvivalKMDiff$compute_asymp_confidence_interval()` (current description `704` chars).
- [ ] TODO #808: Method `InferenceSurvivalKMDiff$compute_asymp_log_rank_two_sided_pval_for_treatment_effect()` (current description `48` chars).
- [ ] TODO #809: Method `InferenceSurvivalKMDiff$compute_asymp_two_sided_pval()` (current description `68` chars).
- [ ] TODO #810: Method `InferenceSurvivalKMDiff$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #811: Method `InferenceSurvivalKMDiff$compute_estimate()` (current description `58` chars).
- [ ] TODO #812: Method `InferenceSurvivalKMDiff$compute_rand_confidence_interval()` (current description `63` chars).
- [ ] TODO #813: Method `InferenceSurvivalKMDiff$new()` (current description `112` chars).
- [ ] TODO #814: Topic `InferenceSurvivalKMDiff` (current description `296` chars).

### `InferenceSurvivalLogRank.Rd`

- [ ] TODO #815: Method `InferenceSurvivalLogRank$clone()` (current description `57` chars).
- [ ] TODO #816: Method `InferenceSurvivalLogRank$compute_asymp_confidence_interval()` (current description `135` chars).
- [ ] TODO #817: Method `InferenceSurvivalLogRank$compute_asymp_log_rank_two_sided_pval_for_treatment_effect()` (current description `77` chars).
- [ ] TODO #818: Method `InferenceSurvivalLogRank$compute_asymp_two_sided_pval()` (current description `75` chars).
- [ ] TODO #819: Method `InferenceSurvivalLogRank$compute_estimate_with_bootstrap_weights()` (current description `77` chars).
- [ ] TODO #820: Method `InferenceSurvivalLogRank$compute_estimate()` (current description `88` chars).
- [ ] TODO #821: Method `InferenceSurvivalLogRank$compute_rand_confidence_interval()` (current description `152` chars).
- [ ] TODO #822: Method `InferenceSurvivalLogRank$new()` (current description `92` chars).
- [ ] TODO #823: Topic `InferenceSurvivalLogRank` (current description `435` chars).

### `InferenceSurvivalRestrictedMeanDiff.Rd`

- [ ] TODO #824: Method `InferenceSurvivalRestrictedMeanDiff$approximate_bootstrap_distribution_beta_hat_T()` (current description `66` chars).
- [ ] TODO #825: Method `InferenceSurvivalRestrictedMeanDiff$clone()` (current description `57` chars).
- [ ] TODO #826: Method `InferenceSurvivalRestrictedMeanDiff$compute_asymp_confidence_interval()` (current description `262` chars).
- [ ] TODO #827: Method `InferenceSurvivalRestrictedMeanDiff$compute_asymp_two_sided_pval()` (current description `48` chars).
- [ ] TODO #828: Method `InferenceSurvivalRestrictedMeanDiff$compute_estimate_with_bootstrap_weights()` (current description `76` chars).
- [ ] TODO #829: Method `InferenceSurvivalRestrictedMeanDiff$compute_estimate()` (current description `58` chars).
- [ ] TODO #830: Method `InferenceSurvivalRestrictedMeanDiff$compute_rand_confidence_interval()` (current description `63` chars).
- [ ] TODO #831: Method `InferenceSurvivalRestrictedMeanDiff$new()` (current description `116` chars).
- [ ] TODO #832: Topic `InferenceSurvivalRestrictedMeanDiff` (current description `296` chars).

### `InferenceSurvivalStratCoxPHRegr.Rd`

- [ ] TODO #833: Method `InferenceSurvivalStratCoxPHRegr$clone()` (current description `57` chars).
- [ ] TODO #834: Method `InferenceSurvivalStratCoxPHRegr$compute_estimate_with_bootstrap_weights()` (current description `57` chars).
- [ ] TODO #835: Method `InferenceSurvivalStratCoxPHRegr$compute_rand_confidence_interval()` (current description `148` chars).
- [ ] TODO #836: Method `InferenceSurvivalStratCoxPHRegr$new()` (current description `103` chars).
- [ ] TODO #837: Topic `InferenceSurvivalStratCoxPHRegr` (current description `304` chars).

### `InferenceSurvivalWeibullRegr.Rd`

- [ ] TODO #838: Method `InferenceSurvivalWeibullRegr$clone()` (current description `57` chars).
- [ ] TODO #839: Method `InferenceSurvivalWeibullRegr$compute_asymp_confidence_interval()` (current description `60` chars).
- [ ] TODO #840: Method `InferenceSurvivalWeibullRegr$compute_asymp_two_sided_pval()` (current description `58` chars).
- [ ] TODO #841: Method `InferenceSurvivalWeibullRegr$compute_estimate_with_bootstrap_weights()` (current description `51` chars).
- [ ] TODO #842: Method `InferenceSurvivalWeibullRegr$compute_estimate()` (current description `58` chars).
- [ ] TODO #843: Method `InferenceSurvivalWeibullRegr$new()` (current description `49` chars).
- [ ] TODO #844: Topic `InferenceSurvivalWeibullRegr` (current description `267` chars).

### `inv_logit.Rd`

- [ ] TODO #845: Topic `inv_logit` (current description `77` chars).

### `logit.Rd`

- [ ] TODO #846: Topic `logit` (current description `49` chars).

### `mn_pvalue_cpp.Rd`

- [ ] TODO #847: Topic `mn_pvalue_cpp` (current description `175` chars).

### `newcombe_independent_ci_cpp.Rd`

- [ ] TODO #848: Topic `newcombe_independent_ci_cpp` (current description `172` chars).

### `ols_hc2_post_fit_cpp.Rd`

- [ ] TODO #849: Topic `ols_hc2_post_fit_cpp` (current description `87` chars).

### `ordinal_gcomp_post_fit_cpp.Rd`

- [ ] TODO #850: Topic `ordinal_gcomp_post_fit_cpp` (current description `115` chars).

### `pocock_simon_assign_and_update_cpp.Rd`

- [ ] TODO #851: Topic `pocock_simon_assign_and_update_cpp` (current description `95` chars).

### `pocock_simon_assign_cpp.Rd`

- [ ] TODO #852: Topic `pocock_simon_assign_cpp` (current description `85` chars).

### `pocock_simon_redraw_w_cpp.Rd`

- [ ] TODO #853: Topic `pocock_simon_redraw_w_cpp` (current description `89` chars).

### `robust_negbinreg.Rd`

- [ ] TODO #854: Topic `robust_negbinreg` (current description `154` chars).

### `robust_survreg_with_surv_object.Rd`

- [ ] TODO #855: Topic `robust_survreg_with_surv_object` (current description `190` chars).

### `robust_survreg.Rd`

- [ ] TODO #856: Topic `robust_survreg` (current description `210` chars).

### `sample_mode.Rd`

- [ ] TODO #857: Topic `sample_mode` (current description `70` chars).

### `set_cold_start_dispatch_policy.Rd`

- [ ] TODO #858: Topic `set_cold_start_dispatch_policy` (current description `75` chars).

### `set_num_cores.Rd`

- [ ] TODO #859: Topic `set_num_cores` (current description `279` chars).

### `set_optimization_dispatch_policy.Rd`

- [ ] TODO #860: Topic `set_optimization_dispatch_policy` (current description `79` chars).

### `set_parallel_dispatch_policy.Rd`

- [ ] TODO #861: Topic `set_parallel_dispatch_policy` (current description `278` chars).

### `set_warm_start_dispatch_policy.Rd`

- [ ] TODO #862: Topic `set_warm_start_dispatch_policy` (current description `75` chars).

### `SimulationFramework.Rd`

- [ ] TODO #863: Method `SimulationFramework$clear_all_intermediate_data_and_gc()` (current description `199` chars).
- [ ] TODO #864: Method `SimulationFramework$clone()` (current description `57` chars).
- [ ] TODO #865: Method `SimulationFramework$get_all_intermediate_data()` (current description `143` chars).
- [ ] TODO #866: Method `SimulationFramework$new()` (current description `99` chars).
- [ ] TODO #867: Method `SimulationFramework$run()` (current description `151` chars).
- [ ] TODO #868: Topic `SimulationFramework` (current description `417` chars).

### `SimulationFrameworkReport.Rd`

- [ ] TODO #869: Method `SimulationFrameworkReport$clone()` (current description `57` chars).
- [ ] TODO #870: Method `SimulationFrameworkReport$get_errors()` (current description `53` chars).
- [ ] TODO #871: Method `SimulationFrameworkReport$get_results()` (current description `36` chars).
- [ ] TODO #872: Method `SimulationFrameworkReport$new()` (current description `106` chars).
- [ ] TODO #873: Method `SimulationFrameworkReport$print()` (current description `38` chars).
- [ ] TODO #874: Method `SimulationFrameworkReport$summarize()` (current description `43` chars).
- [ ] TODO #875: Topic `SimulationFrameworkReport` (current description `370` chars).

### `summary_glm_lean.Rd`

- [ ] TODO #876: Topic `summary_glm_lean` (current description `104` chars).

### `toggle_asserts.Rd`

- [ ] TODO #877: Topic `toggle_asserts` (current description `731` chars).

### `transform_cont_y_based_on_response_type.Rd`

- [ ] TODO #878: Topic `transform_cont_y_based_on_response_type` (current description `267` chars).

### `unset_num_cores.Rd`

- [ ] TODO #879: Topic `unset_num_cores` (current description `187` chars).
