# Inference Hierarchy Tier + Mixin Restructuring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ground the `EDI` package's `Inference*` R6 hierarchy in a 4-tier likelihood taxonomy (`InferenceNoLik` / `InferenceQuasiLik` / `InferencePartialLik` / `InferenceFullLik`), convert every cross-cutting concern that isn't actually about likelihood tier into an orthogonal mixin, and rename single-class-use "mixins" to "ext" (they're file-splitting, not mixing).

**Architecture:** Four new abstract tier classes sit under `InferenceAsympLik`, each fixing the tier-appropriate default capability flags. `InferenceParamBootstrap`, `InferenceCountLikelihood`, and `InferenceAsympLikStdModCache` — each of which today conflates "sets capability defaults" with "shares plumbing code" — are split into a genuine Pattern-1 mixin (the plumbing) plus a thin deprecated pass-through class (so every existing concrete kernel keeps working, unmigrated, until its own batch lands). Concrete kernels migrate tier-by-tier in later batches, each batch independently testable.

**Tech Stack:** R6 classes, R's own S3-free OOP via `R6::R6Class`, Pattern-1 mixins (plain `list(public=, private=)` spliced via `c()`/`modifyList()`), roxygen2 docs regenerated via `Rscript fast_roxygenize.R` (per project convention — never call `roxygen2::roxygenize()` directly).

## Global Constraints

- Every phase must leave `devtools::load_all("EDI")` and the existing test suite (`comprehensive_tests.R` / targeted tests) passing before moving to the next phase — no phase depends on a later phase's classes existing yet.
- No concrete kernel's numeric output (estimate, CI, p-value, bootstrap distribution) may change during Phases 0–3. Only `inherit=` lines and mixin names change; behavior is identical.
- Never delete or rename `InferenceMixinCordeiroFerrariApprox`, `InferenceMixinLemonteGradientApprox`, `InferenceMixinOffOptimumLikelihoodEval` — explicitly out of scope, leave these three dead-scaffold mixins untouched.
- Always run `Rscript fast_roxygenize.R` from the repo root after any class rename or `inherit=` change, never `roxygen2::roxygenize()` directly.
- `lock_objects = FALSE` must be preserved on every R6 class definition touched (project convention for test subclassing — see repo memory `feedback_r6_subclass_lock_objects`).
- Follow existing file/class naming convention exactly: mixin files are `inference_mixin_<snake_case>.R` defining `InferenceMixin<PascalCase>`; ext files (new) are `inference_ext_<snake_case>.R` defining `InferenceExt<PascalCase>`.

---

## Phase 0: Rename single-class-use mixins to "Ext"

**Why first:** these renames are pure mechanical find-and-replace with zero behavior risk, and Phase 2/3 need to reference the post-rename names, so doing this first avoids touching the same lines twice.

### Task 0.1: Rename the 5 `InferenceNonParamBootstrap`-only mixins

**Files:**
- Rename: `inference_mixin_minimum_volatility_selector.R` → `inference_ext_minimum_volatility_selector.R`
- Rename: `inference_mixin_exchangeable_resampling_units.R` → `inference_ext_exchangeable_resampling_units.R`
- Rename: `inference_mixin_bca_bootstrap_ci.R` → `inference_ext_bca_bootstrap_ci.R`
- Rename: `inference_mixin_prw_subsampling.R` → `inference_ext_prw_subsampling.R`
- Rename: `inference_mixin_m_out_of_n_bootstrap.R` → `inference_ext_m_out_of_n_bootstrap.R`
- Modify: `EDI/R/inference_all_abstract_non_param_boot.R`
- Modify: `EDI/R/mixin_contracts.R` (the `EDI_MIXIN_CONTRACTS` / `EDI_MIXIN_COMPOSITIONS` registry)

**Interfaces:**
- Consumes: nothing new — these are leaf files with no dependents besides `InferenceNonParamBootstrap`.
- Produces: `InferenceExtMinimumVolatilitySelector`, `InferenceExtExchangeableResamplingUnits`, `InferenceExtBcaBootstrapCI`, `InferenceExtPRWSubsampling`, `InferenceExtMOutOfNBootstrap` — same `list(public=, private=)` shape as before, only the top-level assignment name changes.

- [ ] **Step 1: Rename the files with git mv**

```bash
cd EDI/R
git mv inference_mixin_minimum_volatility_selector.R inference_ext_minimum_volatility_selector.R
git mv inference_mixin_exchangeable_resampling_units.R inference_ext_exchangeable_resampling_units.R
git mv inference_mixin_bca_bootstrap_ci.R inference_ext_bca_bootstrap_ci.R
git mv inference_mixin_prw_subsampling.R inference_ext_prw_subsampling.R
git mv inference_mixin_m_out_of_n_bootstrap.R inference_ext_m_out_of_n_bootstrap.R
```

- [ ] **Step 2: Rename the class assignment inside each renamed file**

In each file, change the top-level assignment and any self-referential doc comments, e.g. in `inference_ext_minimum_volatility_selector.R`:

```r
# before
InferenceMixinMinimumVolatilitySelector = list(
# after
InferenceExtMinimumVolatilitySelector = list(
```

Repeat for the other four: `InferenceMixinExchangeableResamplingUnits` → `InferenceExtExchangeableResamplingUnits`, `InferenceMixinBcaBootstrapCI` → `InferenceExtBcaBootstrapCI`, `InferenceMixinPRWSubsampling` → `InferenceExtPRWSubsampling`, `InferenceMixinMOutOfNBootstrap` → `InferenceExtMOutOfNBootstrap`. Also update each file's leading roxygen comment block (the `#' A Pattern-1 mixin...` description) to say "Pattern-1 file-split extension, spliced into exactly one class (`InferenceNonParamBootstrap`)" instead of describing it as shared.

- [ ] **Step 3: Update the two splice sites in `inference_all_abstract_non_param_boot.R`**

```r
# before (public splice, ~line 128)
InferenceMixinMOutOfNBootstrap$public,
InferenceMixinPRWSubsampling$public,

# after
InferenceExtMOutOfNBootstrap$public,
InferenceExtPRWSubsampling$public,
```

```r
# before (private splice, ~line 698)
InferenceMixinBcaBootstrapCI$private,
InferenceMixinExchangeableResamplingUnits$private,
InferenceMixinMinimumVolatilitySelector$private,
InferenceMixinMOutOfNBootstrap$private,
InferenceMixinPRWSubsampling$private,

# after
InferenceExtBcaBootstrapCI$private,
InferenceExtExchangeableResamplingUnits$private,
InferenceExtMinimumVolatilitySelector$private,
InferenceExtMOutOfNBootstrap$private,
InferenceExtPRWSubsampling$private,
```

- [ ] **Step 4: Update `mixin_contracts.R`**

Move the five renamed entries out of `EDI_MIXIN_CONTRACTS` (which is keyed by mixin name and used to guard against silent method-name overwrites across mixins) into a new, parallel `EDI_EXT_CONTRACTS` list with the same `file`/`private_methods`/`private_state` shape, renaming each key. Also remove their entries from `EDI_MIXIN_COMPOSITIONS` (that list is specifically "every base class that combines two or more mixins" — a single-use Ext no longer qualifies once renamed, though `InferenceNonParamBootstrap` still combines multiple *Ext* files, so add a parallel `EDI_EXT_COMPOSITIONS` entry for it instead of deleting the composition fact entirely).

- [ ] **Step 5: Verify no stale references remain**

Run: `grep -rn "InferenceMixinMinimumVolatilitySelector\|InferenceMixinExchangeableResamplingUnits\|InferenceMixinBcaBootstrapCI\|InferenceMixinPRWSubsampling\|InferenceMixinMOutOfNBootstrap" EDI/R/`
Expected: no output (zero matches).

- [ ] **Step 6: Reload and run the non-param-bootstrap tests**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
Expected: loads with no errors. Then run whatever test file covers `InferenceNonParamBootstrap` descendants (bootstrap CI tests) and confirm all pass, output unchanged.

- [ ] **Step 7: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_ext_minimum_volatility_selector.R EDI/R/inference_ext_exchangeable_resampling_units.R EDI/R/inference_ext_bca_bootstrap_ci.R EDI/R/inference_ext_prw_subsampling.R EDI/R/inference_ext_m_out_of_n_bootstrap.R EDI/R/inference_all_abstract_non_param_boot.R EDI/R/mixin_contracts.R EDI/NAMESPACE EDI/man/
git commit -m "Rename 5 single-class InferenceNonParamBootstrap mixins to Ext"
```

---

### Task 0.2: Rename the 3 `InferenceAsympLik`-only mixins

**Files:**
- Rename: `inference_mixin_ci_inversion.R` → `inference_ext_ci_inversion.R`
- Rename: `inference_mixin_information_matrix.R` → `inference_ext_information_matrix.R`
- Rename: `inference_mixin_likelihood_test_memoization.R` → `inference_ext_likelihood_test_memoization.R`
- Modify: `EDI/R/inference_all_abstract_asymp_lik.R`
- Modify: `EDI/R/mixin_contracts.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `InferenceExtCIInversion`, `InferenceExtInformationMatrix`, `InferenceExtLikelihoodTestMemoization`.

- [ ] **Step 1: git mv the three files**

```bash
cd EDI/R
git mv inference_mixin_ci_inversion.R inference_ext_ci_inversion.R
git mv inference_mixin_information_matrix.R inference_ext_information_matrix.R
git mv inference_mixin_likelihood_test_memoization.R inference_ext_likelihood_test_memoization.R
```

- [ ] **Step 2: Rename the class assignment in each file**

`InferenceMixinCIInversion` → `InferenceExtCIInversion`, `InferenceMixinInformationMatrix` → `InferenceExtInformationMatrix`, `InferenceMixinLikelihoodTestMemoization` → `InferenceExtLikelihoodTestMemoization`. Update each file's roxygen header the same way as Task 0.1 Step 2.

- [ ] **Step 3: Update the splice site in `inference_all_abstract_asymp_lik.R`**

```r
# before (~line 278)
private = c(InferenceMixinCIInversion$private, InferenceMixinInformationMatrix$private, InferenceMixinLikelihoodTestMemoization$private, list(

# after
private = c(InferenceExtCIInversion$private, InferenceExtInformationMatrix$private, InferenceExtLikelihoodTestMemoization$private, list(
```

Also update the three doc-comment cross-references inside this file (the `get_likelihood_test_spec` doc block that mentions `InferenceMixinInformationMatrix`, `InferenceMixinLikelihoodTestMemoization`, `InferenceMixinOffOptimumLikelihoodEval` by name — update the first two, leave `InferenceMixinOffOptimumLikelihoodEval` as-is per the Global Constraints).

- [ ] **Step 4: Update `mixin_contracts.R`** the same way as Task 0.1 Step 4, for these three names.

- [ ] **Step 5: Verify no stale references**

Run: `grep -rn "InferenceMixinCIInversion\|InferenceMixinInformationMatrix\|InferenceMixinLikelihoodTestMemoization" EDI/R/`
Expected: no output.

- [ ] **Step 6: Reload and run the asymptotic-likelihood test suite**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
Expected: loads cleanly. Run the score/gradient/LR/Bartlett test files and confirm all pass unchanged (these three Ext files back the entire `InferenceAsympLik` test-dispatch machinery, so this is the highest-traffic rename in Phase 0 — check carefully).

- [ ] **Step 7: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_ext_ci_inversion.R EDI/R/inference_ext_information_matrix.R EDI/R/inference_ext_likelihood_test_memoization.R EDI/R/inference_all_abstract_asymp_lik.R EDI/R/mixin_contracts.R EDI/NAMESPACE EDI/man/
git commit -m "Rename 3 single-class InferenceAsympLik mixins to Ext"
```

---

### Task 0.3: Rename the 2 `InferenceParamBootstrap`-only mixins

**Files:**
- Rename: `inference_mixin_param_bootstrap_estimate.R` → `inference_ext_param_bootstrap_estimate.R`
- Rename: `inference_mixin_bartlett_approx.R` → `inference_ext_bartlett_approx.R`
- Modify: `EDI/R/inference_all_abstract_param_boot.R`
- Modify: `EDI/R/mixin_contracts.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `InferenceExtParamBootstrapEstimate`, `InferenceExtBartlettApprox`. **Note for Task 2.1:** Phase 2 will further restructure the class that hosts these two — do this rename first so Phase 2 references the final names.

- [ ] **Step 1: git mv both files**

```bash
cd EDI/R
git mv inference_mixin_param_bootstrap_estimate.R inference_ext_param_bootstrap_estimate.R
git mv inference_mixin_bartlett_approx.R inference_ext_bartlett_approx.R
```

- [ ] **Step 2: Rename the class assignment in each file**

`InferenceMixinParamBootstrapEstimate` → `InferenceExtParamBootstrapEstimate`, `InferenceMixinBartlettApprox` → `InferenceExtBartlettApprox`. Update roxygen headers (note `inference_ext_bartlett_approx.R`'s header currently says "Splice into a daughter class (in practice, just `InferenceParamBootstrap` itself, once)" — keep that sentence, it's still true after the rename).

- [ ] **Step 3: Update the splice site in `inference_all_abstract_param_boot.R`**

```r
# before (~line 543)
private = c(InferenceMixinBartlettApprox$private, InferenceMixinParamBootstrapEstimate$private, list(

# after
private = c(InferenceExtBartlettApprox$private, InferenceExtParamBootstrapEstimate$private, list(
```

Also update the doc-comment cross-reference near line 554 (`InferenceMixinBartlettApprox, get_bartlett_factor_approx()`).

- [ ] **Step 4: Update `mixin_contracts.R`** the same way as prior tasks.

- [ ] **Step 5: Verify no stale references**

Run: `grep -rn "InferenceMixinParamBootstrapEstimate\|InferenceMixinBartlettApprox" EDI/R/`
Expected: no output.

- [ ] **Step 6: Reload and run parametric-bootstrap tests**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
Run the parametric-bootstrap LR / Bartlett-approx test files (these cover every `InferenceParamBootstrap` descendant — Poisson, NegBin, OLS, Cox OneLik, etc.) and confirm identical output.

- [ ] **Step 7: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_ext_param_bootstrap_estimate.R EDI/R/inference_ext_bartlett_approx.R EDI/R/inference_all_abstract_param_boot.R EDI/R/mixin_contracts.R EDI/NAMESPACE EDI/man/
git commit -m "Rename 2 single-class InferenceParamBootstrap mixins to Ext"
```

---

### Task 0.4: Rename the 2 `InferenceRand`-only mixins and the 1 `InferenceAbstractQuantileRandCI`-only mixin

**Files:**
- Rename: `inference_mixin_custom_randomization_statistic.R` → `inference_ext_custom_randomization_statistic.R`
- Rename: `inference_mixin_sequential_mc_pval.R` → `inference_ext_sequential_mc_pval.R`
- Rename: `inference_mixin_quantile_rand_ci.R` → `inference_ext_quantile_rand_ci.R`
- Modify: `EDI/R/inference_all_abstract_rand.R`
- Modify: `EDI/R/inference_all_abstract_quantile_rand_ci.R`
- Modify: `EDI/R/mixin_contracts.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `InferenceExtCustomRandomizationStatistic`, `InferenceExtSequentialMCPval`, `InferenceExtQuantileRandCI`.

- [ ] **Step 1: git mv all three files**

```bash
cd EDI/R
git mv inference_mixin_custom_randomization_statistic.R inference_ext_custom_randomization_statistic.R
git mv inference_mixin_sequential_mc_pval.R inference_ext_sequential_mc_pval.R
git mv inference_mixin_quantile_rand_ci.R inference_ext_quantile_rand_ci.R
```

- [ ] **Step 2: Rename the class assignment in each file** (same pattern as prior tasks).

- [ ] **Step 3: Update `inference_all_abstract_rand.R`**

```r
# before (~line 9)
public = c(InferenceMixinCustomRandomizationStatistic$public, list(
# after
public = c(InferenceExtCustomRandomizationStatistic$public, list(

# before (~line 353)
private = c(InferenceMixinCustomRandomizationStatistic$private, InferenceMixinSequentialMCPval$private, list(
# after
private = c(InferenceExtCustomRandomizationStatistic$private, InferenceExtSequentialMCPval$private, list(
```

- [ ] **Step 4: Update `inference_all_abstract_quantile_rand_ci.R`**

```r
# before (~line 12)
public = utils::modifyList(as.list(InferenceMixinQuantileRandCI$public), list(
# after
public = utils::modifyList(as.list(InferenceExtQuantileRandCI$public), list(

# before (~line 22)
private = utils::modifyList(as.list(InferenceMixinQuantileRandCI$private), list())
# after
private = utils::modifyList(as.list(InferenceExtQuantileRandCI$private), list())
```

- [ ] **Step 5: Update `mixin_contracts.R`** the same way as prior tasks.

- [ ] **Step 6: Verify no stale references**

Run: `grep -rn "InferenceMixinCustomRandomizationStatistic\|InferenceMixinSequentialMCPval\|InferenceMixinQuantileRandCI" EDI/R/`
Expected: no output.

- [ ] **Step 7: Reload and run randomization-inference and quantile-regression tests**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
Confirm identical output on the randomization p-value / sequential Monte Carlo tests and the quantile-regression rand-CI tests.

- [ ] **Step 8: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_ext_custom_randomization_statistic.R EDI/R/inference_ext_sequential_mc_pval.R EDI/R/inference_ext_quantile_rand_ci.R EDI/R/inference_all_abstract_rand.R EDI/R/inference_all_abstract_quantile_rand_ci.R EDI/R/mixin_contracts.R EDI/NAMESPACE EDI/man/
git commit -m "Rename remaining 3 single-class mixins (Rand x2, QuantileRandCI) to Ext"
```

---

## Phase 1: Add `is_a_*` characterization functions to the 4 true mixins

**Why:** these are the only mixins spliced into more than one class, so they're the only ones that legitimately represent a reusable characterization of "what kind of class is this" — worth a real predicate function.

### Task 1.1: Convert `InferenceMixinKKPassThrough`'s field flag to a function

**Files:**
- Modify: `EDI/R/inference_mixin_kk_passthrough.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `private$is_a_kk_passthrough_design()` returning `TRUE`, callable on any of the ~21 classes that splice this mixin.

- [ ] **Step 1: Add the function alongside the existing field**

```r
# in inference_mixin_kk_passthrough.R private list, near kk_passthrough = TRUE
kk_passthrough = TRUE,
is_a_kk_passthrough_design = function() TRUE,
```

Keep the existing `kk_passthrough` field too — do not delete it in this task, since other code may already read it directly (check with `grep -rn "private\$kk_passthrough\b" EDI/R/` and only remove the field in a later, separate cleanup task if nothing reads it besides the doc comment).

- [ ] **Step 2: Update the roxygen header**

Change `#' Capability flag: \code{private$kk_passthrough == TRUE}.` to also mention `\code{private$is_a_kk_passthrough_design()}` as the preferred check going forward.

- [ ] **Step 3: Reload and confirm no regression**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
This is purely additive (new function, nothing removed), so every existing test must still pass unchanged.

- [ ] **Step 4: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_mixin_kk_passthrough.R EDI/NAMESPACE EDI/man/
git commit -m "Add is_a_kk_passthrough_design() characterization to InferenceMixinKKPassThrough"
```

---

### Task 1.2: Add `is_a_kk_compound_estimator()` to `InferenceMixinKKPassThroughCompound`

**Files:**
- Modify: `EDI/R/inference_mixin_kk_passthrough_compound.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `private$is_a_kk_compound_estimator()` returning `TRUE`, callable on both classes built from this mixin (`InferenceKKPassThroughCompound`, `InferenceKKPassThroughCompoundNoParamBootstrap`) and, after Phase 3, on every concrete IVWC class.

- [ ] **Step 1: Add the function**

```r
# in inference_mixin_kk_passthrough_compound.R private list
is_a_kk_compound_estimator = function() TRUE,
```

- [ ] **Step 2: Reload and confirm no regression**

```bash
Rscript -e 'devtools::load_all("EDI")'
```

- [ ] **Step 3: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_mixin_kk_passthrough_compound.R EDI/NAMESPACE EDI/man/
git commit -m "Add is_a_kk_compound_estimator() characterization to InferenceMixinKKPassThroughCompound"
```

---

### Task 1.3: Add `is_a_gee_family()` to `InferenceMixinKKGEEShared`

**Files:**
- Modify: `EDI/R/inference_mixin_kk_gee_shared.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `private$is_a_gee_family()` returning `TRUE`, callable on the 4 classes that splice this mixin (count/incidence/ordinal/proportion GEE variants).

- [ ] **Step 1: Add the function**

```r
# in inference_mixin_kk_gee_shared.R private list
is_a_gee_family = function() TRUE,
```

- [ ] **Step 2: Reload and confirm no regression**

```bash
Rscript -e 'devtools::load_all("EDI")'
```

- [ ] **Step 3: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_mixin_kk_gee_shared.R EDI/NAMESPACE EDI/man/
git commit -m "Add is_a_gee_family() characterization to InferenceMixinKKGEEShared"
```

---

### Task 1.4: Add `is_a_glmm_family()` to `InferenceMixinKKGLMMShared`

**Files:**
- Modify: `EDI/R/inference_mixin_kk_glmm_shared.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `private$is_a_glmm_family()` returning `TRUE`, callable on the 3 classes that splice this mixin (continuous/count/ordinal GLMM variants).

- [ ] **Step 1: Add the function**

```r
# in inference_mixin_kk_glmm_shared.R private list
is_a_glmm_family = function() TRUE,
```

- [ ] **Step 2: Reload and confirm no regression**

```bash
Rscript -e 'devtools::load_all("EDI")'
```

- [ ] **Step 3: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_mixin_kk_glmm_shared.R EDI/NAMESPACE EDI/man/
git commit -m "Add is_a_glmm_family() characterization to InferenceMixinKKGLMMShared"
```

---

## Phase 2: Introduce the four tier classes (purely additive)

### Task 2.1: Create `InferenceNoLik`, `InferenceQuasiLik`, `InferencePartialLik`, `InferenceFullLik`

**Files:**
- Create: `EDI/R/inference_all_abstract_no_lik.R`
- Create: `EDI/R/inference_all_abstract_quasi_lik.R`
- Create: `EDI/R/inference_all_abstract_partial_lik.R`
- Create: `EDI/R/inference_all_abstract_full_lik.R`
- Test: `EDI/tests/testthat/test-inference-tier-classes.R` (new)

**Interfaces:**
- Consumes: `InferenceAsympLik` (parent for all four).
- Produces: `InferenceNoLik`, `InferenceQuasiLik`, `InferencePartialLik`, `InferenceFullLik` — each instantiable-but-abstract R6 classes (never constructed directly by user code, only subclassed), each fixing `private$supports_likelihood_tests()` and (for `InferenceFullLik` only) `private$supports_fisher_information()` to their tier default. No existing concrete class inherits from these yet — this task is purely additive and cannot break anything.

- [ ] **Step 1: Write the failing test**

Verified against the actual codebase first: `grep -rln "R6::R6Class(\"Test\|inherit = InferenceAsympLik\b" EDI/tests/testthat/` returns **no files** — this repo has no precedent for unit-testing an abstract `Inference*` class in isolation via a throwaway subclass; every test exercises abstract behavior indirectly through a real concrete kernel (Poisson, OLS, etc.). Do not invent that pattern here. Instead, since these four tier classes only override zero-argument private methods (`supports_likelihood_tests`), test them via the R6 generator's `$private_methods` list directly — this needs no instantiation at all, so there's no `des_obj`-argument problem to work around:

```r
# EDI/tests/testthat/test-inference-tier-classes.R
test_that("the four tier classes exist", {
  expect_true(exists("InferenceNoLik"))
  expect_true(exists("InferenceQuasiLik"))
  expect_true(exists("InferencePartialLik"))
  expect_true(exists("InferenceFullLik"))
})

test_that("InferenceNoLik and InferenceQuasiLik default supports_likelihood_tests() to FALSE", {
  expect_false(InferenceNoLik$private_methods$supports_likelihood_tests())
  expect_false(InferenceQuasiLik$private_methods$supports_likelihood_tests())
})

test_that("InferencePartialLik and InferenceFullLik default supports_likelihood_tests() to TRUE", {
  expect_true(InferencePartialLik$private_methods$supports_likelihood_tests())
  expect_true(InferenceFullLik$private_methods$supports_likelihood_tests())
})

test_that("all four tiers inherit from InferenceAsympLik", {
  expect_identical(InferenceNoLik$get_inherit(), InferenceAsympLik)
  expect_identical(InferenceQuasiLik$get_inherit(), InferenceAsympLik)
  expect_identical(InferencePartialLik$get_inherit(), InferenceAsympLik)
  expect_identical(InferenceFullLik$get_inherit(), InferenceAsympLik)
})
```

`$private_methods$supports_likelihood_tests()` calls the raw unbound function directly — valid here specifically because this particular private method takes no arguments and never references `self`/`private`/`super` internally (confirmed in Steps 3–6 below); this shortcut is not safe to use on private methods that do.

- [ ] **Step 2: Run the test file to verify it fails**

Run: `Rscript -e 'devtools::load_all("EDI"); testthat::test_file("EDI/tests/testthat/test-inference-tier-classes.R")'`
Expected: FAIL — `InferenceNoLik` etc. don't exist yet.

- [ ] **Step 3: Write `InferenceNoLik`**

```r
# EDI/R/inference_all_abstract_no_lik.R
#' No-Likelihood Tier
#'
#' @name InferenceNoLik
#' @description Tier 0 of the likelihood taxonomy: no likelihood object of any
#' kind backs this family (rank-based tests, sign test, KM-based survival
#' tests, exact tests, g-computation without a model likelihood, quantile
#' regression). Wald/permutation inference only — never score, gradient, LR,
#' or Bartlett-corrected tests, and never parametric-bootstrap simulation.
#'
#' @keywords internal
InferenceNoLik = R6::R6Class("InferenceNoLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(),
	private = list(
		supports_likelihood_tests = function() FALSE
	)
)
```

- [ ] **Step 4: Write `InferenceQuasiLik`**

```r
# EDI/R/inference_all_abstract_quasi_lik.R
#' Quasi-Likelihood Tier
#'
#' @name InferenceQuasiLik
#' @description Tier 1 of the likelihood taxonomy: an estimating-equation /
#' working-variance construct (GEE, quasi-Poisson, robust/sandwich-corrected
#' regression). Not a true density, so there is no valid likelihood-ratio
#' statistic without ad hoc correction. Mechanically identical to
#' \code{InferenceNoLik} (Wald-only via sandwich variance) but semantically
#' distinct, and the intended home for quasi-likelihood-specific shared logic
#' (e.g. quasi-deviance helpers) as it gets built out.
#'
#' @keywords internal
InferenceQuasiLik = R6::R6Class("InferenceQuasiLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(),
	private = list(
		supports_likelihood_tests = function() FALSE
	)
)
```

- [ ] **Step 5: Write `InferencePartialLik`**

```r
# EDI/R/inference_all_abstract_partial_lik.R
#' Partial/Conditional-Likelihood Tier
#'
#' @name InferencePartialLik
#' @description Tier 2 of the likelihood taxonomy: a genuine likelihood
#' function with rigorously valid score/Wald/gradient/LR asymptotics (Cox
#' 1975 proves this for the partial likelihood case), but one that does not
#' fully specify the joint density — a nuisance function is left unspecified
#' (the baseline hazard for Cox, the conditioning distribution for
#' conditional logistic regression). Houses Cox PH, stratified Cox,
#' conditional logistic regression, adjacent-category conditional logit,
#' LWA-Cox. May additionally splice \code{InferenceMixinParamBootstrapSimulate}
#' (see Phase 3) if a semiparametric \code{simulate_under_lik_null()} is
#' implemented — that mixin's gate is "at least this tier," not "exactly
#' \code{InferenceFullLik}."
#'
#' @keywords internal
InferencePartialLik = R6::R6Class("InferencePartialLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(),
	private = list(
		supports_likelihood_tests = function() TRUE
	)
)
```

- [ ] **Step 6: Write `InferenceFullLik`**

```r
# EDI/R/inference_all_abstract_full_lik.R
#' Full-Likelihood Tier
#'
#' @name InferenceFullLik
#' @description Tier 3 of the likelihood taxonomy: a fully specified
#' generative density. Superset of \code{InferencePartialLik}'s capability,
#' plus trivial parametric simulation from p(y|theta,X), plus eligibility for
#' exact analytic corrections (Fisher information, Bartlett-exact) since the
#' density is known in closed form. Houses Poisson, NegBin, Hurdle/ZI
#' mixtures, logistic, probit, cauchit, cloglog, ordered-probit,
#' proportional-odds, stereotype-logit, log-binomial, beta,
#' zero-one-inflated-beta, OLS, Weibull, Clayton-copula-with-parametric-
#' margins, dependent-censoring transform.
#'
#' @keywords internal
InferenceFullLik = R6::R6Class("InferenceFullLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(),
	private = list(
		supports_likelihood_tests = function() TRUE
	)
)
```

Note: do NOT set `supports_fisher_information` to `TRUE` unconditionally here — per the two capability lattices, not every Full-tier family has Fisher information wired up (e.g. Poisson does, NegBin/OLS/Weibull don't yet). Leave the default `FALSE` from `InferenceAsympLik` and let individual Full-tier daughters opt in, same as today.

- [ ] **Step 7: Run the test file to verify it passes**

Run: `Rscript -e 'devtools::load_all("EDI"); testthat::test_file("EDI/tests/testthat/test-inference-tier-classes.R")'`
Expected: PASS.

- [ ] **Step 8: Run the full existing test suite to confirm zero regression**

This task adds four new classes that nothing else references yet — the entire existing suite must pass byte-for-byte unchanged. Run the project's full test entry point (`comprehensive_tests.R` or `devtools::test("EDI")`, whichever this repo uses) and confirm no output changed anywhere.

- [ ] **Step 9: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_all_abstract_no_lik.R EDI/R/inference_all_abstract_quasi_lik.R EDI/R/inference_all_abstract_partial_lik.R EDI/R/inference_all_abstract_full_lik.R EDI/tests/testthat/test-inference-tier-classes.R EDI/NAMESPACE EDI/man/
git commit -m "Add InferenceNoLik/QuasiLik/PartialLik/FullLik tier classes (additive, unused by any concrete class yet)"
```

---

## Phase 3: Convert `InferenceParamBootstrap` into a mixin behind a deprecated pass-through

**Why this shape:** ~50+ existing concrete classes inherit `InferenceParamBootstrap` directly today. Converting it into `InferenceMixinParamBootstrapSimulate` while keeping `InferenceParamBootstrap` alive as a thin subclass that splices the new mixin means **zero existing files need to change in this phase** — they keep inheriting `InferenceParamBootstrap` exactly as before, and it keeps working identically. Only new/migrated Partial-tier and Full-tier classes need to reach for the mixin directly.

### Task 3.1: Extract `InferenceParamBootstrap`'s public/private surface into `InferenceMixinParamBootstrapSimulate`

**Files:**
- Create: `EDI/R/inference_mixin_param_bootstrap_simulate.R`
- Modify: `EDI/R/inference_all_abstract_param_boot.R`

**Interfaces:**
- Consumes: `InferenceExtBartlettApprox`, `InferenceExtParamBootstrapEstimate` (from Phase 0, Task 0.3) — the new mixin still splices these in, unchanged.
- Produces: `InferenceMixinParamBootstrapSimulate` — a Pattern-1 mixin (`list(public=, private=)`) containing everything `InferenceParamBootstrap` currently defines directly (its own `public` methods like `compute_lik_ratio_bootstrap_two_sided_pval`, `get_last_param_bootstrap_diagnostics`, etc., and its own `private` methods like `supports_lik_ratio_param_bootstrap = function() FALSE`, `simulate_param_boot_*`, `simulate_under_lik_null`, `run_param_bootstrap_replicates`, etc.), PLUS the already-spliced `InferenceExtBartlettApprox`/`InferenceExtParamBootstrapEstimate` content folded in at the top level (so a class splicing `InferenceMixinParamBootstrapSimulate` gets everything in one splice, matching how `InferenceMixinKKPassThroughCompound` already folds `InferenceMixinKKPassThrough` conceptually — though here just flatten the three source lists into one).

- [ ] **Step 1: Read the full current `InferenceParamBootstrap` class body**

```bash
cat EDI/R/inference_all_abstract_param_boot.R
```
Copy every `public` method and every `private` method/field verbatim — this is a pure cut-and-paste, no logic changes.

- [ ] **Step 2: Create the new mixin file with the extracted content**

```r
# EDI/R/inference_mixin_param_bootstrap_simulate.R
#' Mixin for Parametric-Bootstrap Likelihood-Ratio Simulation
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private}
#' slots) providing parametric-bootstrap LR calibration
#' (\code{compute_lik_ratio_bootstrap_two_sided_pval()},
#' \code{compute_lik_ratio_bootstrap_confidence_interval()}) and the
#' Monte-Carlo Bartlett correction, for any \code{InferencePartialLik} or
#' \code{InferenceFullLik} daughter that implements
#' \code{simulate_under_lik_null()}. Splice into a daughter class via
#' \code{public = c(InferenceMixinParamBootstrapSimulate$public, list(...))}
#' and \code{private = c(InferenceMixinParamBootstrapSimulate$private, list(...))}.
#'
#' Requires at least \code{InferencePartialLik} tier — a genuine likelihood is
#' needed to define the LR statistic, but the density need not be fully
#' specified (Cox's Breslow-plug-in semiparametric simulation qualifies, same
#' as a fully parametric \code{InferenceFullLik} simulation).
#'
#' Capability flag: \code{private$supports_lik_ratio_param_bootstrap()}.
#'
#' @keywords internal
#' @noRd
InferenceMixinParamBootstrapSimulate = list(
	public = list(
		# <paste every public method from the current InferenceParamBootstrap
		#  class body here, verbatim, unchanged>
	),
	private = c(InferenceExtBartlettApprox$private, InferenceExtParamBootstrapEstimate$private, list(
		# <paste every private method/field from the current
		#  InferenceParamBootstrap class body here, verbatim, unchanged,
		#  including supports_lik_ratio_param_bootstrap = function() FALSE>
	))
)
```

The two placeholder comments above must be filled with the **actual current file contents** when this step is executed — copy-paste exactly, do not summarize or abbreviate any method body.

- [ ] **Step 3: Replace `InferenceParamBootstrap`'s body with a thin pass-through that splices the new mixin**

```r
# EDI/R/inference_all_abstract_param_boot.R (replacing the entire prior body)
#' Parametric-Bootstrap-Capable Likelihood Inference (deprecated pass-through)
#'
#' @description \strong{Deprecated shape.} All behavior now lives in
#' \code{InferenceMixinParamBootstrapSimulate}. This class exists only so
#' that every pre-restructuring concrete daughter keeps inheriting
#' \code{InferenceParamBootstrap} unchanged. New Partial-tier or Full-tier
#' classes should splice \code{InferenceMixinParamBootstrapSimulate} directly
#' instead of inheriting this class — see the tier-migration batches in this
#' plan's later phases.
#'
#' @keywords internal
InferenceParamBootstrap = R6::R6Class("InferenceParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = InferenceMixinParamBootstrapSimulate$public,
	private = InferenceMixinParamBootstrapSimulate$private
)
```

- [ ] **Step 4: Verify byte-identical behavior**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
Run the FULL parametric-bootstrap test suite (every Poisson/NegBin/OLS/Cox-OneLik/Copula-OneLik/Frailty-OneLik/CondLogit-OneLik test that exercises `compute_lik_ratio_bootstrap_two_sided_pval` etc.) — every single one of these classes still inherits `InferenceParamBootstrap`, so output must be bit-identical to before this task. This is the highest-risk task in Phase 3 because it touches the most-inherited abstract class in the package; do not proceed to Phase 4 until this is green.

- [ ] **Step 5: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_mixin_param_bootstrap_simulate.R EDI/R/inference_all_abstract_param_boot.R EDI/NAMESPACE EDI/man/
git commit -m "Extract InferenceParamBootstrap into InferenceMixinParamBootstrapSimulate behind a deprecated pass-through"
```

---

## Phase 4: Convert `InferenceCountLikelihood` into a mixin behind deprecated pass-throughs

**Why this shape:** 4 existing concrete classes (`InferenceCountHurdleNegBin`, `InferenceCountNegBin`, `InferenceCountPoisson`, the zero-augmented-Poisson family) inherit `InferenceCountLikelihood` directly today. Same reasoning as Phase 3: extract the shared body into a real mixin, keep both existing classes as thin pass-throughs so nothing else needs to change in this phase.

### Task 4.1: Extract `inference_count_likelihood_public`/`private` into `InferenceMixinCountLikelihoodPlumbing`

**Files:**
- Create: `EDI/R/inference_mixin_count_likelihood_plumbing.R`
- Modify: `EDI/R/inference_all_abstract_count_likelihood.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `InferenceMixinCountLikelihoodPlumbing` (a real `list(public=, private=)` mixin) containing exactly what `inference_count_likelihood_public`/`inference_count_likelihood_private` currently contain — count-specific parameter packing, warm starts, and the `mark_count_likelihood_block_asymp_nonestimable()`-gated dispatch overrides for asymp/wald/score/gradient/lik_ratio/bootstrap paths.

- [ ] **Step 1: Create the new mixin file with the exact current content**

```r
# EDI/R/inference_mixin_count_likelihood_plumbing.R
#' Mixin for Count-Specific Likelihood Plumbing
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private}
#' slots) providing count-specific parameter packing, warm starts, and
#' likelihood-test dispatch shared by every count-based Full-likelihood
#' family (Poisson, Negative Binomial, Zero-Inflated, Hurdle). Splice into a
#' daughter class via
#' \code{public = c(InferenceMixinCountLikelihoodPlumbing$public, list(...))}
#' and
#' \code{private = c(InferenceMixinCountLikelihoodPlumbing$private, list(...))}.
#'
#' @keywords internal
#' @noRd
InferenceMixinCountLikelihoodPlumbing = list(
	public = list(
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			if (private$testing_type == "wald") {
				private$shared(estimate_only = FALSE)
				if (is.finite(private$cached_values$s_beta_hat_T %||% NA_real_)) {
					return(private$compute_z_or_t_ci_from_s_and_df(alpha))
				}
			}
			super$compute_asymp_confidence_interval(alpha)
		},
		compute_asymp_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			if (private$testing_type == "wald") {
				private$shared(estimate_only = FALSE)
				if (is.finite(private$cached_values$s_beta_hat_T %||% NA_real_)) {
					return(private$compute_z_or_t_two_sided_pval_from_s_and_df(delta))
				}
			}
			super$compute_asymp_two_sided_pval(delta)
		},
		compute_wald_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_wald_two_sided_pval(delta)
		},
		compute_wald_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_wald_confidence_interval(alpha)
		},
		compute_score_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_score_two_sided_pval(delta)
		},
		compute_score_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_score_confidence_interval(alpha)
		},
		compute_lik_ratio_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_lik_ratio_two_sided_pval(delta)
		},
		compute_lik_ratio_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_lik_ratio_confidence_interval(alpha)
		},
		compute_gradient_two_sided_pval = function(delta = 0){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_gradient_two_sided_pval(delta)
		},
		compute_gradient_confidence_interval = function(alpha = 0.05){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_gradient_confidence_interval(alpha)
		},
		compute_lik_ratio_bootstrap_two_sided_pval = function(delta = 0, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			super$compute_lik_ratio_bootstrap_two_sided_pval(
				delta = delta,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate
			)
		},
		compute_lik_ratio_bootstrap_confidence_interval = function(alpha = 0.05, B = 199, show_progress = FALSE, min_number_usable_samples = 5L, max_attempts_per_replicate = 2L, root_tolerance = NULL, max_root_iterations = 8L){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			super$compute_lik_ratio_bootstrap_confidence_interval(
				alpha = alpha,
				B = B,
				show_progress = show_progress,
				min_number_usable_samples = min_number_usable_samples,
				max_attempts_per_replicate = max_attempts_per_replicate,
				root_tolerance = root_tolerance,
				max_root_iterations = max_root_iterations
			)
		}
	),
	private = list(
		shared = function(estimate_only = FALSE){
			if (estimate_only &&
					!is.null(private$cached_values$beta_hat_T) &&
					is.finite(as.numeric(private$cached_values$beta_hat_T)[1L])) {
				return(invisible(NULL))
			}
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T) && is.finite(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			model_output = private$generate_mod(estimate_only = estimate_only)
			private$cached_mod = model_output

			if (is.null(model_output)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = model_output$beta_hat_T %||% model_output$b[2]

			if (!is.null(model_output$params) || !is.null(model_output$b)) {
				private$set_fit_warm_start(
					as.numeric(model_output$params %||% model_output$b),
					type = if (!is.null(model_output$params)) "params" else "beta",
					fisher = model_output$fisher_information %||% model_output$XtWX,
					weights = model_output$w %||% model_output$mu,
					force_pd = TRUE
				)
			}

			if (estimate_only) return(invisible(NULL))

			ssq = model_output$ssq_b_j %||% model_output$ssq_b_2
			if (!is.null(ssq) && is.finite(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$df = model_output$df %||% Inf
		},
		get_standard_error = function(){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$shared(estimate_only = FALSE)
			if (isTRUE(private$supports_information_preference())) {
				se = tryCatch(private$compute_standard_error_from_information_matrix(), error = function(e) NA_real_)
				if (is.finite(se)) return(se)
			}
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df %||% Inf
		},
		get_backend_warm_start_args = function(expected_length, expected_fisher_dim = expected_length) {
			private$get_optimal_warm_start_config(expected_length, expected_fisher_dim)
		},
		supports_lik_ratio_param_bootstrap = function() TRUE,
		supports_likelihood_tests = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			NULL
		},
		compute_score_two_sided_pval_impl = function(delta){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(NA_real_)
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		compute_score_confidence_interval_impl = function(alpha){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			private$invert_test_pval_confidence_interval(alpha)
		},
		compute_gradient_confidence_interval_impl = function(alpha){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			private$invert_test_pval_confidence_interval(alpha)
		},
		compute_lik_ratio_confidence_interval_impl = function(alpha){
			if (private$mark_count_likelihood_block_asymp_nonestimable()) return(private$count_likelihood_missing_ci(alpha))
			private$invert_test_pval_confidence_interval(alpha)
		},
		count_likelihood_block_asymp_unsupported = function(){
			private$jackknife_block_size_gt_one_unsupported(unit = "auto")
		},
		mark_count_likelihood_block_asymp_nonestimable = function(){
			if (!private$count_likelihood_block_asymp_unsupported()) return(FALSE)
			if (is.null(private$cached_values$beta_hat_T)) {
				try(private$shared(estimate_only = TRUE), silent = TRUE)
			}
			private$cache_nonestimable_se("count_likelihood_asymp_block_size_gt_one_not_supported")
			TRUE
		},
		count_likelihood_missing_ci = function(alpha = 0.05){
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	)
)
```

This is an exact, verbatim transcription of the current `inference_count_likelihood_public`/`inference_count_likelihood_private` bodies — no logic changes.

- [ ] **Step 2: Replace the two class definitions with deprecated pass-throughs**

```r
# EDI/R/inference_all_abstract_count_likelihood.R (replacing everything below the mixin extraction)
#' Count-Specific Likelihood Inference (deprecated pass-through)
#'
#' @name InferenceCountLikelihood
#' @description \strong{Deprecated shape.} All behavior now lives in
#' \code{InferenceMixinCountLikelihoodPlumbing}. Kept so pre-restructuring
#' daughters (Poisson, NegBin, Hurdle, Zero-Inflated/Zero-Augmented) keep
#' inheriting this class unchanged. New Full-tier count families should
#' splice \code{InferenceMixinCountLikelihoodPlumbing} directly onto
#' \code{InferenceFullLik} instead.
#'
#' @keywords internal
InferenceCountLikelihood = R6::R6Class("InferenceCountLikelihood",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = InferenceMixinCountLikelihoodPlumbing$public,
	private = InferenceMixinCountLikelihoodPlumbing$private
)

#' Count-Specific Likelihood Inference Without Parametric LR Bootstrap (deprecated pass-through, unused)
#'
#' @name InferenceCountLikelihoodNoParamBootstrap
#' @description \strong{Deprecated shape, currently zero concrete daughters}
#' (verified via \code{grep -rln "inherit = InferenceCountLikelihoodNoParamBootstrap"}
#' returning no files as of this plan's writing). Kept unchanged/unremoved
#' per this plan's Global Constraints (don't delete disconnected classes
#' without being asked) — do not delete in this task even though it is dead.
#'
#' @keywords internal
#' @noRd
InferenceCountLikelihoodNoParamBootstrap = R6::R6Class("InferenceCountLikelihoodNoParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = InferenceMixinCountLikelihoodPlumbing$public,
	private = InferenceMixinCountLikelihoodPlumbing$private
)
```

- [ ] **Step 3: Verify byte-identical behavior**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
Run every test covering `InferenceCountHurdleNegBin`, `InferenceCountNegBin`, `InferenceCountPoisson`, and the zero-augmented-Poisson family (Hurdle-Poisson, ZIP, ZINB) — confirm bit-identical output, including the Jackknife-block-size-gt-one nonestimability path (`mark_count_likelihood_block_asymp_nonestimable`) and the parametric-bootstrap LR path.

- [ ] **Step 4: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_mixin_count_likelihood_plumbing.R EDI/R/inference_all_abstract_count_likelihood.R EDI/NAMESPACE EDI/man/
git commit -m "Extract InferenceCountLikelihood into InferenceMixinCountLikelihoodPlumbing behind deprecated pass-throughs"
```

---

## Phase 5: Convert `InferenceAsympLikStdModCache` into a mixin behind deprecated pass-throughs

**Why this shape:** this is the single largest branch in the package (17 direct daughters of `InferenceAsympLikStdModCache`, 2 of `InferenceAsympLikStdModCacheNoParamBootstrap` — logistic, probit, cauchit, cloglog, ordered-probit, proportional-odds, stereotype-logit, adjacent-category-logit, log-binomial, beta, zero-one-inflated-beta, Cox, stratified Cox, Weibull, dependent-censoring transform, incidence risk-diff, proportion fractional-logit). Same reasoning as Phases 3–4.

### Task 5.1: Extract `inference_asymp_lik_std_mod_cache_public`/`private` into `InferenceMixinStdModCache`

**Files:**
- Create: `EDI/R/inference_mixin_std_mod_cache.R`
- Modify: `EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `InferenceMixinStdModCache` containing exactly what `inference_asymp_lik_std_mod_cache_public`/`inference_asymp_lik_std_mod_cache_private` currently contain — the generic GLM/KM `shared()`/caching logic, warm-start-aware bootstrap-worker hooks, and the score/gradient/LR dispatch built on `get_likelihood_test_spec()`.

- [ ] **Step 1: Create the new mixin file with the exact current content**

```r
# EDI/R/inference_mixin_std_mod_cache.R
#' Mixin for Standard-Model-Cache GLM/KM Plumbing
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private}
#' slots) providing the generic cached-model-fit plumbing shared by every
#' non-KK GLM-style Full-likelihood family (logistic, probit, cauchit,
#' cloglog, ordered-probit, proportional-odds, stereotype-logit,
#' adjacent-category-logit, log-binomial, beta, zero-one-inflated-beta, Cox,
#' stratified Cox, Weibull, dependent-censoring transform) — this is the
#' single largest branch in the package (see
#' \code{package_metadata/non_ivwc_capability_lattice.html}). Splice into a
#' daughter class via
#' \code{public = c(InferenceMixinStdModCache$public, list(...))} and
#' \code{private = c(InferenceMixinStdModCache$private, list(...))}.
#'
#' @keywords internal
#' @noRd
InferenceMixinStdModCache = list(
	public = list(
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (private$testing_type == "wald") {
				private$shared(estimate_only = FALSE)
				if (is.finite(private$cached_values$s_beta_hat_T %||% NA_real_)) {
					return(private$compute_z_or_t_ci_from_s_and_df(alpha))
				}
			}
			super$compute_asymp_confidence_interval(alpha)
		},
		compute_asymp_two_sided_pval = function(delta = 0){
			if (private$testing_type == "wald") {
				private$shared(estimate_only = FALSE)
				if (is.finite(private$cached_values$s_beta_hat_T %||% NA_real_)) {
					return(private$compute_z_or_t_two_sided_pval_from_s_and_df(delta))
				}
			}
			super$compute_asymp_two_sided_pval(delta)
		}
	),
	private = list(
		supports_likelihood_tests = function(){
			TRUE
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		generate_mod = function(estimate_only = FALSE) stop(class(self)[1], " must implement generate_mod()"),
		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},
		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},
		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			if (isTRUE(private$supports_information_preference())) {
				se = private$compute_standard_error_from_information_matrix()
				if (is.finite(se)) return(se)
			}
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		},
		compute_score_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "score")
		},
		compute_gradient_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		compute_lik_ratio_two_sided_pval_impl = function(delta){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "lik_ratio")
		},
		get_likelihood_test_spec = function(){
			NULL
		},
		make_warm_fit_null_wrapper = function(spec, cache_key){
			last_start = NULL
			last_delta = NULL
			ci_inversion = grepl("_ci$", as.character(cache_key)[1L])
			fit_null_formals = tryCatch(names(formals(spec$fit_null)), error = function(e) character())
			accepts_start = "start" %in% fit_null_formals
			function(delta){
				warm_enabled = isTRUE(private$null_fit_warm_start_enabled) && !ci_inversion
				cache_state = if (warm_enabled) private$get_likelihood_null_warm_state(cache_key) else NULL
				start = if (warm_enabled) last_start else NULL
				if (warm_enabled && is.null(start) && !is.null(cache_state)) start = cache_state$start
				fit = tryCatch(
					if (accepts_start) spec$fit_null(delta, start = start) else spec$fit_null(delta),
					error = function(e) NULL
				)
				extract_start = spec$extract_start %||% function(fit_obj) NULL
				last_start <<- if (warm_enabled && accepts_start) tryCatch(extract_start(fit), error = function(e) NULL) else NULL
				last_delta <<- delta
				if (warm_enabled && accepts_start) {
					private$set_likelihood_null_warm_state(cache_key, delta = delta, start = last_start)
				}
				fit
			}
		},
		compute_likelihood_test_two_sided_pval = function(delta, testing_type, bartlett_B = NULL){
			spec = private$get_likelihood_test_spec()
			if (is.null(spec)) {
				stop(class(self)[1], " does not expose a likelihood-test specification.", call. = FALSE)
			}
			p_value = private$get_memoized_likelihood_test_pval(
				delta = delta,
				testing_type = testing_type,
				spec = spec,
				warm_cache_key = paste0("likelihood_test:", testing_type),
				bartlett_B = bartlett_B
			)
			if (!is.finite(p_value) && !isTRUE(self$is_nonestimable("estimate"))) {
				private$cache_nonestimable_se(paste0(testing_type, "_test_unavailable"))
			}
			p_value
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			has_cached_se = !is.null(private$cached_values$s_beta_hat_T) &&
				length(private$cached_values$s_beta_hat_T) == 1L &&
				isTRUE(is.finite(private$cached_values$s_beta_hat_T))
			if (isTRUE(!is.null(private$cached_values$beta_hat_T) && (estimate_only || has_cached_se))) return(invisible(NULL))

			model_output = private$generate_mod(estimate_only = estimate_only)

			private$cached_mod = model_output
			if (is.null(model_output)) {
				private$cache_nonestimable_estimate("model_fit_unavailable")
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			beta_hat_T = as.numeric(model_output$beta_hat_T %||% model_output$b[2])[1L]
			if (!is.finite(beta_hat_T)) {
				private$cache_nonestimable_estimate("model_treatment_estimate_unavailable")
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T = beta_hat_T

			if (!is.null(model_output$b)) {
				private$set_fit_warm_start(
					as.numeric(model_output$params %||% model_output$b),
					type = if (!is.null(model_output$params)) "params" else "beta",
					fisher = model_output$fisher_information %||% model_output$XtWX,
					weights = model_output$w %||% model_output$mu,
					force_pd = TRUE
				)
			}
			if (estimate_only) return(invisible(NULL))
			ssq = model_output$ssq_b_2 %||% model_output$ssq_b_j
			ssq = if (length(ssq) >= 1L) as.numeric(ssq)[1L] else NA_real_
			private$cached_values$df = model_output$df %||% NA_real_
			if (is.finite(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
				private$clear_nonestimable_state()
			} else {
				private$cache_nonestimable_se("model_standard_error_unavailable")
			}
		},
		get_backend_warm_start_args = function(expected_length, expected_fisher_dim = expected_length) {
			private$get_optimal_warm_start_config(expected_length, expected_fisher_dim)
		}
	)
)
```

- [ ] **Step 2: Replace the two class definitions with deprecated pass-throughs**

```r
# EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R (replacing everything below the mixin extraction)
#' GLM and Kaplan-Meier Inference (deprecated pass-through)
#'
#' @description \strong{Deprecated shape.} All behavior now lives in
#' \code{InferenceMixinStdModCache}. Kept so the 17 pre-restructuring
#' daughters keep inheriting this class unchanged. New Full-tier GLM-style
#' families should splice \code{InferenceMixinStdModCache} directly onto
#' \code{InferenceFullLik} instead.
#'
#' @keywords internal
InferenceAsympLikStdModCache = R6::R6Class("InferenceAsympLikStdModCache",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = InferenceMixinStdModCache$public,
	private = InferenceMixinStdModCache$private
)

#' GLM and Kaplan-Meier Inference Without Parametric LR Bootstrap (deprecated pass-through)
#'
#' @description \strong{Deprecated shape.} All behavior now lives in
#' \code{InferenceMixinStdModCache}. Kept so the 2 pre-restructuring
#' daughters (\code{InferenceIncidRiskDiff}, \code{InferencePropFractionalLogit})
#' keep inheriting this class unchanged.
#'
#' @keywords internal
#' @noRd
InferenceAsympLikStdModCacheNoParamBootstrap = R6::R6Class("InferenceAsympLikStdModCacheNoParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = InferenceMixinStdModCache$public,
	private = InferenceMixinStdModCache$private
)
```

- [ ] **Step 3: Verify byte-identical behavior**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
This is the highest-traffic rename in the whole restructuring (19 total direct daughters). Run the full logistic/probit/cauchit/cloglog/ordered-probit/proportional-odds/stereotype-logit/adjacent-category-logit/log-binomial/beta/zero-one-inflated-beta/Cox/stratified-Cox/Weibull/dependent-censoring/incidence-risk-diff/fractional-logit test suite and confirm bit-identical output before moving on. Do not proceed to Phase 6 until every one of these is green.

- [ ] **Step 4: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_mixin_std_mod_cache.R EDI/R/inference_all_abstract_asymp_lik_std_mod_cache.R EDI/NAMESPACE EDI/man/
git commit -m "Extract InferenceAsympLikStdModCache into InferenceMixinStdModCache behind deprecated pass-throughs"
```

---

## Phase 6: Point `InferenceKKPassThroughCompound(NoParamBootstrap)` at the tier classes

**Why this shape is different from Phases 3–5:** `InferenceMixinKKPassThroughCompound` is already a real mixin (already correctly spliced, not embedded in the host class body) — nothing needs extracting. The only change is what the two host classes inherit from.

### Task 6.1: Re-point the two compound host classes

**Files:**
- Modify: `EDI/R/inference_all_abstract_KK_passthrough_compound.R`

**Interfaces:**
- Consumes: `InferenceFullLik` (Phase 2), `InferenceNoLik` (Phase 2).
- Produces: `InferenceKKPassThroughCompound` now inherits `InferenceFullLik` instead of `InferenceParamBootstrap`; `InferenceKKPassThroughCompoundNoParamBootstrap` now inherits `InferenceNoLik` instead of `InferenceAsympLik`. Every one of the six already-migrated IVWC classes (Task 2 of the six-file migration, done directly in this conversation) plus the 8 originally-correct IVWC classes automatically move with it — no per-file changes needed for any of those 14 classes.

- [ ] **Step 1: Change the two `inherit=` lines**

```r
# before
InferenceKKPassThroughCompound = R6::R6Class("InferenceKKPassThroughCompound",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = inference_kk_passthrough_compound_public,
	private = inference_kk_passthrough_compound_private
)
# after
InferenceKKPassThroughCompound = R6::R6Class("InferenceKKPassThroughCompound",
	lock_objects = FALSE,
	inherit = InferenceFullLik,
	public = inference_kk_passthrough_compound_public,
	private = inference_kk_passthrough_compound_private
)
```

```r
# before
InferenceKKPassThroughCompoundNoParamBootstrap = R6::R6Class("InferenceKKPassThroughCompoundNoParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = inference_kk_passthrough_compound_public,
	private = inference_kk_passthrough_compound_private
)
# after
InferenceKKPassThroughCompoundNoParamBootstrap = R6::R6Class("InferenceKKPassThroughCompoundNoParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceNoLik,
	public = inference_kk_passthrough_compound_public,
	private = inference_kk_passthrough_compound_private
)
```

Note: `InferenceMixinKKPassThroughCompound$private` already sets `supports_likelihood_tests = function() FALSE` (see `inference_kk_passthrough_compound_private`'s splice), which is consistent with — and now redundant with, but not contradicted by — `InferenceNoLik`'s own default of `FALSE`. Leave the explicit override in place; do not remove it in this task (removing redundant flags is exactly the kind of change that should get its own reviewed task, not ride along here).

- [ ] **Step 2: Verify byte-identical behavior for all 14 current daughters**

```bash
Rscript -e 'devtools::load_all("EDI")'
```
Run the full test suite for every class confirmed to inherit `InferenceKKPassThroughCompound` or `InferenceKKPassThroughCompoundNoParamBootstrap` (the 8 originally-correct IVWC classes plus the 6 migrated in this conversation's six-file IVWC migration) and confirm bit-identical output. This is a two-line change with a wide blast radius (14 daughters) — treat the verification step as seriously as the change is small.

- [ ] **Step 3: Regenerate docs and commit**

```bash
Rscript fast_roxygenize.R
git add EDI/R/inference_all_abstract_KK_passthrough_compound.R EDI/NAMESPACE EDI/man/
git commit -m "Point InferenceKKPassThroughCompound(NoParamBootstrap) at InferenceFullLik/InferenceNoLik tiers"
```

---

## Phase 7 and beyond: concrete kernel migration, not detailed task-by-task in this document

The remaining work — migrating every one of the ~130 concrete kernels' `inherit=` line off `InferenceParamBootstrap`/`InferenceAsympLik` and onto the tier classes directly (or onto the now tier-based `InferenceCountLikelihood`/`InferenceAsympLikStdModCache`/`InferenceKKPassThroughCompound(NoParamBootstrap)` pass-throughs, which is actually optional busywork once Phases 3–6 land, since those pass-throughs now already resolve to the correct tier transitively) — is too large to hand-author as individual bite-sized tasks in one document; that would mean thousands of near-identical lines and would violate this plan format's own "no placeholder / no `Similar to Task N`" rule the moment it got abbreviated.

Instead, once Phases 0–6 are merged and green:

1. **Use the two published capability lattices** (`package_metadata/ivwc_capability_lattice.html`, `package_metadata/non_ivwc_capability_lattice.html`) as the authoritative (class, target-tier) mapping — every row in both tables already has an unambiguous tier assignment from the design doc's "Where the original counterexamples land now" section.
2. **Generate per-family migration batches from that mapping** (e.g. "all Full-tier count-family kernels," "all Partial-tier survival kernels") sized at 5–10 files each, and run each batch through this same writing-plans skill separately — do not attempt one giant plan for all ~130 at once. Each batch plan should follow the six-file IVWC migration pattern already executed directly in this conversation (change `inherit=`, remove now-redundant tier-flag overrides only, verify byte-identical test output) as its worked template.
3. **Only after every concrete kernel has moved off the deprecated pass-through classes** should those pass-throughs (`InferenceParamBootstrap`, `InferenceCountLikelihood`, `InferenceAsympLikStdModCache`) be deleted — a final, separate cleanup task, not before. (`InferenceKKPassThroughCompound(NoParamBootstrap)` are not pass-throughs to delete — after Phase 6 they are the permanent tier-anchored compound base classes.)
