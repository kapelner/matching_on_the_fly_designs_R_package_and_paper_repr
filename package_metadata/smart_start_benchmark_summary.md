# Smart-Start Benchmark Summary

Benchmark outputs are in [/tmp/smart_start_reports](/tmp/smart_start_reports:1):

- [initialization_report.csv](/tmp/smart_start_reports/initialization_report.csv:1)
- [ci_inversion_score_report.csv](/tmp/smart_start_reports/ci_inversion_score_report.csv:1)
- [ci_inversion_lik_ratio_report.csv](/tmp/smart_start_reports/ci_inversion_lik_ratio_report.csv:1)
- [randomization_ci_reuse_report.csv](/tmp/smart_start_reports/randomization_ci_reuse_report.csv:1)
- [bootstrap_reuse_report.csv](/tmp/smart_start_reports/bootstrap_reuse_report.csv:1)

The benchmark configuration was:

- `n = 100`
- `p = 10`, interpreted as intercept + treatment + 8 covariates
- direct-fit benchmark compares smart default initialization against the legacy initialization branch
- CI, randomization, and bootstrap benchmarks compare warm-reuse `on` versus `off`

## Main Takeaways

The warm-reuse framework looks worthwhile.

- likelihood-ratio CI inversion improved for every benchmarked family
- score CI inversion improved for Poisson, negative binomial, Weibull, and ordinal logit
- randomization CI improved materially for Poisson
- bootstrap improved materially for Poisson

The direct initialization story is weaker and currently mixed.

- ordinal logit improved
- Poisson was effectively unchanged
- logistic, negative binomial, and Weibull were slower under the smart-start branch in this benchmark

Because of that, logistic, negative binomial, and Weibull smart default
initialization have since been removed. Those families now use their legacy
default starts unless the caller supplies an explicit start or opts into the
smart branch.

That split matters. The strongest performance case for this work is not the one-shot OLS-style initialization by itself. It is the combination of explicit-start plumbing plus warm reuse across sequences of nearby fits.

## Initialization Regressions

The surprising result is that the smart default start was not uniformly better even on plain estimation:

- negative binomial: `0.023s` smart vs `0.011s` legacy
- Weibull: `0.012s` smart vs `0.008s` legacy

Negative binomial is the clearest regression. It also took many more iterations in this run:

- negbin smart median iterations: `28`
- negbin legacy median iterations: `4`

That suggests the current smart NB initialization is not yet aligned with the optimizer as well as the legacy branch on this design. Since the smart NB start currently mixes OLS slopes with the legacy dispersion start, the most likely issue is that the `beta` and `theta` blocks are not jointly well calibrated.

Weibull also regressed modestly. Since the smart branch uses uncensored-only OLS while the legacy branch uses the older all-row start, this may be a case where the uncensored-only start is statistically more principled but not always cheaper in finite samples.

## Warm-Reuse Results

The stronger results came from reuse across repeated related fits.

Likelihood-ratio CI inversion improved across all benchmarked families:

- logistic: about `84%` faster
- Poisson: about `29%` faster
- negbin: about `51%` faster
- Weibull: about `44%` faster
- ordinal logit: about `74%` faster

Score CI inversion also improved in most families:

- Poisson: about `30%` faster
- negbin: about `63%` faster
- Weibull: about `46%` faster
- ordinal logit: about `50%` faster

The exception was logistic score inversion, which was slightly slower with reuse enabled in this run. That is small enough that it may be noise, but it should be rechecked.

Randomization CI and bootstrap also showed the expected benefit from reuse:

- randomization CI, Poisson: about `65%` faster
- bootstrap CI, Poisson: about `54%` faster

These are the most convincing validation points for the warm-start framework.

## Interpretation

At this point the engineering conclusion should be:

- keep the explicit-start and warm-reuse framework
- treat CI inversion, randomization CI, and bootstrap reuse as the main performance win
- do not assume the new smart default initialization is uniformly better family by family

The immediate follow-up work should focus on the initialization regressions:

- re-tune the negative binomial smart start first
- then re-check Weibull on a broader grid of sample sizes and signal strengths
- keep benchmarking initialization separately from repeated-fit reuse, because they are behaving differently
