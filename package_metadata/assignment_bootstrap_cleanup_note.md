# Cleanup Note On `bootstrap_type = "assignment"`

## Short Version

The current `bootstrap_type = "assignment"` branch in
[inference_all_abstract_non_param_boot.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract_non_param_boot.R:47)
is probably **misnamed**.

It does **not** currently generate fresh treatment assignments from the design’s
assignment mechanism. Instead, it calls `resample_assignment()` on the design,
and the design-level implementations in:

- [design_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_abstract.R:491)
- [design_blocking_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_blocking_abstract.R:209)
- [design_matching_abstract.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/design_matching_abstract.R:140)

simply resample the observed `(w, y, dead, m)` tuples by index.

So the current behavior is still much closer to **nonparametric case
resampling** than to a true **assignment / reassignment bootstrap**.

## Why The Name Is Misleading

“Assignment resampling” suggests something like:

1. hold subjects and observed outcomes fixed
2. redraw treatment assignments from the original randomization mechanism
3. recompute the estimator under those new assignments

That is not what the current code does.

What it actually does is more like:

1. sample indices with replacement
2. replace `w`, `y`, `dead`, and maybe `m` by those resampled rows

So the current branch is not:

- a true reassignment/bootstrap over the assignment mechanism
- a parametric bootstrap
- a randomization test

It is just a different path to empirical resampling.

## Recommended Cleanup

### Option 1: Rename It

If the current behavior is kept, rename:

```text
bootstrap_type = "assignment"
```

to something closer to what it actually does, e.g.

- `"resample_observed_tuples"`
- `"case_resample_assignments_and_outcomes"`
- `"empirical_assignment_resample"`

That would reduce conceptual confusion immediately.

### Option 2: Implement A Real Reassignment Bootstrap

If the package actually wants assignment-based resampling, add a separate mode
whose semantics are:

1. keep `X`, `y`, and censoring/outcome data fixed
2. redraw `w` from the design’s assignment rule
3. preserve design-specific constraints:
   - blocks
   - matched sets
   - cluster assignment
   - sequential rules where possible

This should be a separate API from the current empirical bootstrap.

For example:

- `bootstrap_type = "reassignment"`
- `bootstrap_type = "design_reassignment"`

## Recommendation

The cleanest short-term step is:

1. keep the current code behavior
2. rename `bootstrap_type = "assignment"`
3. reserve a new name like `"reassignment"` for a future true design-based
   assignment redraw

That keeps terminology aligned with what the code actually does and makes later
bootstrap-LR / design-based resampling work less confusing.
