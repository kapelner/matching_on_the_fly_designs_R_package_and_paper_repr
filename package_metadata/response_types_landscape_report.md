# Response Types Across Experimental Literatures

## Scope

This report summarizes how common the package's response-type families are
across broad branches of experimental research:

- `continuous`
- `incidence` (binary)
- `count`
- `proportion`
- `survival`
- `ordinal`
- `nominal`

The target literatures are intentionally broad:

- medical science and clinical trials
- economics field experiments
- political experiments
- sociology and other social-science field experiments
- chemistry and materials science
- physics and related measurement-heavy experimental sciences

## Important Limitation

There is no single authoritative census of outcome types across all
experimental literatures.

So this report should be read as a **qualitative landscape**:

- `dominant`
- `common`
- `present but secondary`
- `specialized`
- `rare`

The evidence is strongest for clinical trials, where review papers exist for
binary, ordinal, and time-to-event outcomes. It is weaker but still informative
for the social and physical sciences, where method overviews and domain examples
are more common than outcome-type frequency audits.

## Short Answer

Across the experimental sciences as a whole:

- `continuous` is the most universal response type
- `incidence` is extremely common in medicine and very common in social-policy
  and behavioral experiments
- `survival` is common and methodologically central in clinical research, but
  much less common outside it
- `ordinal` is common in some medical and social-science subfields, but not as
  universal as continuous or binary outcomes
- `count` is common in parts of biology, epidemiology, political behavior, and
  physics, but not usually the dominant general-purpose endpoint family
- `nominal` is less dominant than binary or continuous outcomes in medicine, but
  very natural in social, political, and choice-based experiments
- `proportion` is important but comparatively specialized; it is often used when
  the outcome is inherently bounded in `[0, 1]` or is a share/fraction

If the package wants to mirror broad experimental practice, the current
response-type menu is already well chosen. The biggest missing "general"
category is `nominal`, not because it dominates all fields, but because it is
important in several literatures and conceptually distinct from `ordinal`.

## High-Level Ranking

Very roughly, across all the literatures considered:

1. `continuous`
2. `incidence`
3. `ordinal` / `count` / `survival` depending on field
4. `nominal`
5. `proportion`

This ranking is deliberately rough. It is more accurate to think field by
field.

## Landscape Table

| Field | Continuous | Incidence | Count | Proportion | Survival | Ordinal | Nominal |
|---|---|---|---|---|---|---|---|
| Clinical trials / medical RCTs | Dominant | Dominant | Common | Specialized | Common | Common in some areas | Present but secondary |
| Epidemiology / public health experiments | Dominant | Dominant | Common | Common | Common | Common | Present |
| Economics field experiments | Dominant | Common | Common | Common | Rare | Common | Common |
| Political experiments | Common | Common | Common | Common | Rare | Common | Common |
| Sociology / social-policy experiments | Common | Common | Common | Common | Rare | Common | Common |
| Chemistry / materials science / engineering DoE | Dominant | Present | Present | Present | Rare | Present | Present |
| Physics / detector / counting sciences | Dominant | Present | Common | Rare | Rare | Rare | Present but secondary |

## By Response Type

### 1. Continuous

This is the most universal response family in experimental work.

Why:

- many experimental measurements are naturally quantitative
- classical design-of-experiments and regression workflows are built around
  continuous responses
- continuous outcomes are information-rich and statistically convenient

In clinical trials, continuous endpoints remain standard, especially for lab
values, symptom scales treated approximately continuously, physiological
measures, and repeated biomarkers. A clinical-trial methods paper notes that the
classical linear model is widely used for continuous outcomes:
https://link.springer.com/article/10.1186/s12874-022-01534-8

In chemistry and materials science, continuous responses are often the default.
The nanocrystal DoE guide lists size, yield, and polydispersity as common
continuous responses:
https://pubs.acs.org/doi/10.1021/acs.chemmater.2c02924

JMP's DoE guidance states that, wherever possible, continuous responses are
preferred because they contain more information than categorical responses:
https://www.jmp.com/en/statistics-knowledge-portal/design-of-experiments/key-design-of-experiments-concepts/response-variables

Verdict:

- `continuous` is the broadest and most universally used response type across
  experimental literatures.

### 2. Incidence (binary)

Binary outcomes are exceptionally important in medicine and common far beyond
it.

In clinical trials, binary outcomes are a major endpoint family. Rombach et al.
reviewed randomized trial reports and note that many trials use binary primary
outcomes; they cite earlier work finding that around half of trials calculated
sample size from a binary outcome:
https://link.springer.com/article/10.1186/s12916-020-01598-7

Binary endpoints are also common in:

- development and policy experiments: participation vs non-participation,
  employed vs unemployed, enrolled vs not enrolled
- political experiments: turnout vs no turnout, donation vs no donation,
  persuasion success vs no success
- social experiments more broadly: completion / compliance / take-up indicators

Verdict:

- `incidence` is one of the core endpoint families in medicine and is also very
  common in social and behavioral experimentation.

### 3. Count

Count outcomes are common but field-dependent.

They arise naturally when the outcome is:

- number of events
- number of visits
- number of crimes
- number of publications
- number of clicks / messages / contacts
- particle or photon counts

In medicine and public health, counts are common for admissions, episodes,
events, and utilization.

In political and social experiments, counts are common for:

- turnout contacts
- campaign interactions
- number of actions taken
- program participation episodes

In physics and detector sciences, count data are foundational. For example,
Lass et al. discuss Poisson and multinomial count-data analysis in neutron
scattering:
https://journals.sagepub.com/doi/10.3233/JNR-190145

Verdict:

- `count` is common and sometimes central, but it is not as universally
  dominant as `continuous` or as broadly standard in clinical trials as binary
  endpoints.

### 4. Proportion

Proportion outcomes are important but relatively specialized.

They are natural when the outcome is:

- a fraction of successes out of a denominator
- a share
- a bounded rate in `[0,1]`
- a composition component

Examples:

- adherence proportion
- market share
- vote share
- treated fraction in cluster summaries
- conversion rate

In clinical research, proportions are common but often get recoded into:

- binary outcomes
- counts with denominators
- or continuous approximations

In economics, political science, and social experiments, proportions often
appear at aggregate levels:

- vote share
- treatment take-up share
- participation rate

Verdict:

- `proportion` is a useful and important specialized family, but it is usually
  less dominant than continuous, binary, or sometimes count endpoints.

### 5. Survival

Survival or time-to-event outcomes are methodologically central in clinical
research and much less common elsewhere.

In medical trials, time-to-event outcomes are a major family. Jachno et al.
reviewed high-impact trial reports with time-to-event primary outcomes and found
Kaplan-Meier curves in 98%, log-rank tests in 88%, and Cox models in 97% of the
reviewed set:
https://link.springer.com/article/10.1186/s12874-019-0749-1

Another review found 373 RCTs with time-to-event outcomes versus 361 with
binary outcomes in a set of 734 eligible RCTs from major journals:
https://link.springer.com/article/10.1186/1471-2288-9-21

That does not mean survival is dominant everywhere. Outside medicine and some
event-history applications in the social sciences, it is much less central.

Verdict:

- `survival` is a major and indispensable family in clinical trials and some
  epidemiology/event-history research, but not a general-purpose response type
  across all experimental disciplines.

### 6. Ordinal

Ordinal outcomes are common in some literatures and increasingly important in
clinical trials where multiple patient states matter.

Selman et al. reviewed RCTs with ordinal outcomes in four major medical
journals and identified 144 eligible studies in 2017-2022, showing that
ordinal endpoints are clearly substantial in the medical literature:
https://link.springer.com/article/10.1186/s13063-024-08072-2

They also note high-profile examples such as:

- modified Rankin scale
- Glasgow Outcome Scale-Extended
- WHO COVID-19 Clinical Progression Scale

In social science, ordinal outcomes are also common:

- Likert-scale attitudes
- satisfaction ratings
- perceived fairness / trust / ideology scales

Verdict:

- `ordinal` is not as universal as continuous or binary, but it is common and
  often substantively important in medicine and survey-based experiments.

### 7. Nominal

Nominal outcomes are not usually the dominant endpoint family in clinical
trials, but they are clearly legitimate and often natural in other
experimental literatures.

The Trials protocol on ordinal outcomes explicitly distinguishes nominal
outcomes as categorical and unranked:
https://link.springer.com/article/10.1186/s13063-023-07262-8

In medicine, nominal outcomes appear in:

- mutually exclusive competing categories
- diagnostic classes
- adverse-event type categories
- treatment response categories without a defensible ordering

Schmid, Trikalinos, and Olkin develop network meta-analysis for unordered
categorical outcomes and apply it to 17 trials with mutually exclusive nominal
categories:
https://pubmed.ncbi.nlm.nih.gov/26052655/

In social and political experiments, nominal outcomes are especially natural:

- vote choice
- party identification
- occupational / schooling status category
- choice among alternatives

Examples:

- vote choice in a randomized political advertising field experiment:
  https://journals.sagepub.com/doi/10.1177/20531680221076901
- party affiliation in a field experiment:
  https://www.nber.org/papers/w15365

Verdict:

- `nominal` is less dominant than binary and continuous outcomes in medicine,
  but it is clearly important in social, political, and choice-based
  experiments and is the main missing general response family from the package.

## By Field

### Clinical trials and medical science

Best-supported conclusions:

- `incidence` is dominant
- `continuous` is dominant
- `survival` is common and methodologically central
- `ordinal` is common in some specialties and has grown in visibility
- `nominal` is present but secondary
- `count` and `proportion` are common but more specialized

Evidence:

- binary outcomes are a major trial family:
  https://link.springer.com/article/10.1186/s12916-020-01598-7
- time-to-event outcomes are common in major journals:
  https://link.springer.com/article/10.1186/1471-2288-9-21
  and
  https://link.springer.com/article/10.1186/s12874-019-0749-1
- ordinal outcomes are clearly substantial:
  https://link.springer.com/article/10.1186/s13063-024-08072-2

### Economics field experiments

Best-supported conclusions:

- `continuous` and `incidence` are very common
- `count`, `proportion`, and `nominal` are also common depending on the
  substantive domain
- `survival` is relatively rare

Why:

- field experiments study earnings, spending, hours worked, scores, and other
  continuous outcomes
- they also study take-up, employment, compliance, and participation as binary
  outcomes
- many applied endpoints are shares, counts, or category choices

Levitt and List's overview emphasizes the breadth of modern field experiments:
https://www.nber.org/papers/w14356

Verdict:

- economics uses a very mixed outcome palette, but continuous and binary
  outcomes probably dominate overall.

### Political experiments

Best-supported conclusions:

- `incidence`, `proportion`, `count`, `ordinal`, and `nominal` all matter
- `continuous` is also common when using scales or margins
- `survival` is rare

Examples:

- turnout or participation: binary
- number of actions taken: count
- vote share: proportion
- ideology / trust / feeling thermometers: ordinal or treated-as-continuous
- vote choice / party affiliation: nominal

Relevant examples and reviews:

- field experiments on political behavior:
  https://www.annualreviews.org/doi/pdf/10.1146/annurev.polisci.12.060107.154037
- vote choice field experiment:
  https://journals.sagepub.com/doi/10.1177/20531680221076901
- party affiliation field experiment:
  https://www.nber.org/papers/w15365

Verdict:

- political experimentation is one of the clearest cases for supporting
  `nominal`, alongside binary, proportion, and ordinal outcomes.

### Sociology and broader social-science field experiments

Best-supported conclusions:

- `incidence` and `continuous` are common
- `ordinal` and `nominal` are also very natural
- `count` and `proportion` are common in behavioral and policy contexts
- `survival` is relatively uncommon

Baldassarri and Abascal's review emphasizes the breadth of field experiments
across the social sciences:
https://www.annualreviews.org/doi/10.1146/annurev-soc-073014-112445

Verdict:

- social-science experiments use a broad mix of outcome types, especially
  binary, continuous, ordinal, and nominal.

### Chemistry, materials science, and engineering DoE

Best-supported conclusions:

- `continuous` is dominant
- discrete outcomes (`incidence`, `ordinal`, `nominal`) are present but
  secondary
- `count` and `proportion` occur in some settings
- `survival` is rare

The nanocrystal DoE guide is explicit: common responses include size,
polydispersity, yield, shape, and crystal phase, with continuous responses like
size/yield and discrete responses like shape/phase:
https://pubs.acs.org/doi/10.1021/acs.chemmater.2c02924

JMP's DOE guidance also advises using continuous responses where possible:
https://www.jmp.com/en/statistics-knowledge-portal/design-of-experiments/key-design-of-experiments-concepts/response-variables

Verdict:

- in chemistry and engineering DOE, continuous outcomes dominate far more than
  in medicine or social science.

### Physics and counting sciences

Best-supported conclusions:

- `continuous` is dominant for many measurement problems
- `count` is especially important in detector, nuclear, particle, and
  scattering settings
- `nominal`, `ordinal`, and `proportion` are much less central
- `survival` is rare outside specialized reliability/event-history contexts

Count-data methodology is fundamental in parts of physics:

- neutron scattering count-data analysis:
  https://journals.sagepub.com/doi/10.3233/JNR-190145
- experimental investigation of count distributions in nuclear measurement:
  https://www.sciencedirect.com/science/article/pii/016890029290209M

Verdict:

- physics is one of the strongest reasons to keep `count` as a first-class
  response family.

## Implications For The Package

If the package is meant to reflect broad experimental practice, the current
response families are well aligned:

- `continuous`: absolutely essential
- `incidence`: absolutely essential
- `count`: essential
- `proportion`: justified specialized family
- `survival`: essential for clinical-trial relevance
- `ordinal`: clearly justified

The biggest conceptual omission is:

- `nominal`

not because nominal outcomes are the single most common outcome type across all
fields, but because:

- they are common enough in important literatures
- they are qualitatively distinct from `ordinal`
- and they are not well represented by coercing category labels into numeric
  scales

## Additional Response Families Worth Noting

The package's current menu is already a sensible first-order summary of common
experimental outcomes. Still, several additional families are common enough in
practice that they are worth naming explicitly.

These should be read as **extra families or structures** that sit beyond the
current scalar response menu, not as tiny variants such as two-step survival or
other minor permutations.

### 1. Longitudinal / repeated-measures outcomes

These are outcomes observed repeatedly over time for the same unit:

- repeated biomarkers in clinical trials
- panel survey outcomes in political and social experiments
- repeated production / assay measurements in chemistry and engineering
- repeated detector or instrument readings in physics

This is one of the most common structures in experimental science.

Implementation difficulty in the package: `hard`

Why:

- the package is currently organized mostly around one realized outcome vector
  per experimental run
- repeated measures require a new data container, not just a new
  `response_type`
- many designs are still compatible, but the inference layer would need new
  estimands, covariance handling, and likely mixed-model / GEE / panel-summary
  style paths
- `SimulationFramework` would need a repeated-outcome data generator and a
  coherent truth object, which may no longer be a scalar `betaT`

Practical conclusion:

- this is important scientifically, but it is probably a second-generation
  architectural expansion, not a small additive response type.

### 2. Multivariate / vector outcomes

These occur when the experiment targets multiple outcomes jointly:

- co-primary clinical endpoints
- biomarker panels
- educational or social outcomes measured as a bundle
- joint physics or chemistry response surfaces

Implementation difficulty in the package: `very hard`

Why:

- the package would need to move from scalar estimands toward vector-valued
  estimands
- summary tables, power summaries, bootstrap logic, and inference-path outputs
  would need redesign
- many current abstractions assume a single target parameter or a single
  treatment-effect summary per path

Practical conclusion:

- this is conceptually important, but it is not just "one more response type";
  it is a broad generalization of the package's inferential core.

### 3. Time-to-recurrent-event / panel event-history outcomes

These are not ordinary one-event survival outcomes. They track repeated events
over time:

- repeated hospitalizations
- repeated infections
- repeated failures or breakdowns
- repeated criminal or behavioral events

Implementation difficulty in the package: `hard`

Why:

- current survival support is much closer to first-event survival
- recurrent-event analysis needs counting-process style data, gap-time or
  calendar-time choices, and more elaborate censoring logic
- inference paths would need Andersen-Gill, frailty, WLW/PWP, or simpler
  robust marginal alternatives

Practical conclusion:

- this is a meaningful extension of `survival`, but still substantial enough to
  be its own project.

### 4. Semi-continuous / zero-inflated continuous outcomes

These appear when there is a pile at zero plus a positive continuous tail:

- medical costs
- health-care utilization expenditures
- environmental or laboratory measurements below or near detection
- treatment intensity with many zeros

Implementation difficulty in the package: `moderate to hard`

Why:

- these outcomes are not well represented by plain `continuous`
- a useful implementation would need either hurdle / two-part inference paths
  or a consciously approximate strategy
- `SimulationFramework` would need a generator that separates the zero mass from
  the positive continuous component

Practical conclusion:

- feasible, but only if the package wants to support two-part estimands
  explicitly.

### 5. Compositional / multinomial-share outcomes

These are vector-valued shares that sum to one:

- market-share allocations
- vote-share vectors across parties
- compositional chemistry outputs
- budget-share or time-allocation outcomes

This is different from a single `proportion`, which covers one bounded scalar.

Implementation difficulty in the package: `hard`

Why:

- the response is constrained and multivariate
- inference needs either a multinomial / Dirichlet / log-ratio framework or a
  scalarization strategy
- design-side compatibility is mostly fine, but the inference and summary layer
  would need broad changes

Practical conclusion:

- scientifically useful, especially for political and compositional sciences,
  but not a small extension of the current `proportion` family.

### 6. Bounded discrete scores / scale totals

These are common in practice:

- Likert totals
- symptom indices
- educational test scores with limited support
- composite scales with floor or ceiling effects

Implementation difficulty in the package: `easy to moderate`

Why:

- many users already analyze these as `continuous` or `ordinal`
- a dedicated response family is not strictly necessary
- if formalized, the package would need a clear policy on whether these are
  treated as approximately continuous, ordinal, or something else

Practical conclusion:

- this is more of a documentation and user-guidance issue than a priority for a
  new package-level response type.

### 7. Rank / preference / choice-set outcomes

These arise when subjects rank alternatives or choose among many options:

- ranked candidate preferences
- product or policy rankings
- best-worst or discrete-choice style outcomes

Implementation difficulty in the package: `hard`

Why:

- plain `nominal` does not capture ranking structure
- inference would need rank-ordered logit, Plackett-Luce, choice models, or
  specialized score summaries
- simulation support would need a generator for full or partial rankings

Practical conclusion:

- common enough in political science, economics, and marketing-style
  experiments to matter, but too structurally different to be folded cleanly
  into the current scalar response menu.

### 8. Censored or truncated continuous outcomes

Examples:

- assay measurements below detection limit
- top-coded income-like outcomes
- machine limits in engineering experiments

Implementation difficulty in the package: `moderate`

Why:

- the response is still close to `continuous`, but inference paths should honor
  censoring/truncation
- this is often handled with Tobit-like or censored-regression machinery, or
  with analysis strategies tailored to detection limits

Practical conclusion:

- this is a plausible later extension of `continuous`, but probably not a
  first-priority standalone response family.

### 9. Interval-censored outcomes

These occur when the event is only known to happen between observation times.

Implementation difficulty in the package: `hard`

Why:

- it is a specialized survival/data-collection structure
- inference paths differ materially from ordinary right-censored survival

Practical conclusion:

- better treated as a later survival-subfamily project than as one of the core
  response types.

### 10. Functional / curve-valued outcomes

Examples:

- full spectra
- growth curves
- time-series traces
- dose-response curves treated as whole objects

Implementation difficulty in the package: `very hard`

Why:

- the target is no longer a scalar or small vector outcome
- inference would need basis expansions, functional regression, or carefully
  chosen scalar summaries

Practical conclusion:

- important in some physical sciences, but far outside the package's current
  architecture.

### 11. Network / relational outcomes

These include:

- ties formed in social experiments
- peer-link formation
- diffusion graphs
- group interaction patterns

Implementation difficulty in the package: `very hard`

Why:

- the outcome lives on pairs or graphs, not units
- interference and dependence are usually central rather than incidental
- the current design and inference abstractions are built around unit-level
  responses

Practical conclusion:

- this is a different methodological domain, not a near-term add-on.

## Bottom Line

There is no single global frequency table for response types across all
experimental literatures. But the broad picture is clear:

- `continuous` is the most universal response family
- `incidence` is one of the most common and important families, especially in
  medicine and policy experiments
- `survival` is indispensable in clinical research
- `count` is central in several scientific and behavioral domains
- `ordinal` is common and important in medicine and survey-based experiments
- `proportion` is useful but comparatively specialized
- `nominal` is less dominant than binary or continuous overall, but is
  important enough across social, political, medical, and choice-based research
  that it deserves explicit support if the package aims to cover experimental
  practice broadly

Beyond those core types, the most important extra families are not just more
labels. Many of them require a broader change in package architecture:

- `longitudinal / repeated-measures` and `multivariate` outcomes would push the
  package beyond a scalar-response design
- `recurrent-event`, `interval-censored`, and `choice/ranking` outcomes are
  substantial specialized expansions
- `semi-continuous`, `censored continuous`, and `bounded discrete scores` are
  more feasible medium-term extensions

So for practical package planning:

1. `nominal` is still the clearest next scalar response type to add.
2. `longitudinal` and `multivariate` are scientifically important, but they are
   architectural projects more than ordinary new response types.
3. Several other families are real and common, but they fit better as later
   domain-specific extensions than as immediate additions to the core menu.

## Sources

- Binary outcomes in clinical trials:
  https://link.springer.com/article/10.1186/s12916-020-01598-7
- Ordinal outcomes in clinical trials:
  https://link.springer.com/article/10.1186/s13063-024-08072-2
- Outcome-type taxonomy including nominal:
  https://link.springer.com/article/10.1186/s13063-023-07262-8
- Time-to-event outcomes in clinical trials:
  https://link.springer.com/article/10.1186/1471-2288-9-21
  and
  https://link.springer.com/article/10.1186/s12874-019-0749-1
- Nominal outcomes in trial/meta-analysis methodology:
  https://pubmed.ncbi.nlm.nih.gov/26052655/
- Field experiments in economics:
  https://www.nber.org/papers/w14356
- Field experiments across the social sciences:
  https://www.annualreviews.org/doi/10.1146/annurev-soc-073014-112445
- Political field experiments:
  https://www.annualreviews.org/doi/pdf/10.1146/annurev.polisci.12.060107.154037
- Vote-choice field experiment:
  https://journals.sagepub.com/doi/10.1177/20531680221076901
- Party-affiliation field experiment:
  https://www.nber.org/papers/w15365
- Chemistry/materials DOE response types:
  https://pubs.acs.org/doi/10.1021/acs.chemmater.2c02924
- DOE response-variable guidance:
  https://www.jmp.com/en/statistics-knowledge-portal/design-of-experiments/key-design-of-experiments-concepts/response-variables
- Physics / count-data methodology:
  https://journals.sagepub.com/doi/10.3233/JNR-190145
  and
  https://www.sciencedirect.com/science/article/pii/016890029290209M
