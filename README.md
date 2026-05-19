# EDI

Code for the R package `EDI` and the code to reproduce the simulations in the papers listed below. You can get started by cloning and running `R CMD INSTALL EDI` from the root of this repo. Then try running the code in the folder `package_tests`.

For a local performance build, you can opt into native compiler flags without
changing the default portable build:

```sh
EDI_NATIVE_SPEED=1 R CMD INSTALL EDI
```

This enables `-O3 -march=native -mtune=native` when your compiler supports
those flags. Use it for local benchmarking only; the default build remains the
portable one.

To compare plain native vs native plus link-time optimization:

```sh
EDI_NATIVE_SPEED=1 R CMD INSTALL EDI
EDI_NATIVE_SPEED=1 EDI_NATIVE_LTO=1 R CMD INSTALL EDI
```

`EDI_NATIVE_LTO=1` adds `-flto` on top of the native-speed flags.

To benchmark the three builds back-to-back, run:

```sh
bash scripts/benchmark_build_modes.sh
```

To compare the current working tree, the current tree with vectorization
disabled, and the last committed `HEAD` snapshot across several hot C++ kernels,
run:

```sh
bash scripts/benchmark_simd_matrix.sh
```

This uses `EDI_DISABLE_VECTORIZATION=1` to add `-DEIGEN_DONT_VECTORIZE` and
`-fno-tree-vectorize` for the no-vectorization build.

To benchmark a randomization CI workload at `num_cores = 3` across the same
three build modes, run:

```sh
bash scripts/benchmark_randomization_ci_build_modes.sh
```

### Setting a seed for reproducible output

All three layers of the package accept a seed for deterministic output.

**Design classes** — pass `seed` to `$new()`:

```r
# Fixed design: same seed → same draw_ws_according_to_design() every call
des = DesignFixedBernoulli$new(n = 100, response_type = "continuous", seed = 42)
des$add_all_subjects_to_experiment(X)
w1 = des$draw_ws_according_to_design(r = 500)
w2 = des$draw_ws_according_to_design(r = 500)
identical(w1, w2)  # TRUE

# Sequential design: same seed → same assignment sequence
des = DesignSeqOneByOneBernoulli$new(n = 100, response_type = "continuous", seed = 42)
for (i in seq_len(100)) des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
```

**Inference classes** — call `$set_seed()` after construction:

```r
inf = InferenceAllSimpleMeanDiff$new(des)
inf$set_seed(42)

# Same seed + same num_cores → identical p-value / CI / Bayesian-bootstrap distribution
p  = inf$compute_rand_two_sided_pval(r = 999, show_progress = FALSE)
ci = inf$compute_bootstrap_confidence_interval(B = 999, show_progress = FALSE)
```

Reproducibility is guaranteed only when `num_cores` is the same across runs (cross-core
determinism is out of scope). Set the number of parallel workers before running:

```r
set_num_cores(4)   # or unset_num_cores() for serial
```

**SimulationFramework** — pass `seed` to `$new()`:

```r
sim = SimulationFramework$new(
    response_type = "continuous",
    design_classes_and_params = list(DesignFixedBernoulli),
    inference_classes_and_params = list(InferenceAllSimpleMeanDiff),
    inference_types_and_params = list(asymp_pval = list()),
    n = 50L, Nrep = 200L, seed = 321, ...
)
sim$run()
```

**Note on `duplicate()`:** Objects produced by `$duplicate()` have their seed cleared intentionally. This prevents all parallel inference workers from resetting to the same RNG stream and producing duplicate treatment allocations.

---

### Experimental findings
The `scripts/benchmark_randomization_ci_ordinal_ppo.R`/`scripts/benchmark_randomization_ci_cases.R` experiments show the native-speed flags are beneficial for the heavier Eigen/OpenMP workloads but not universally faster:

- **Ordinal PPO (`InferenceOrdinalMultiPartialProportionalOddsRegr`)** with `r=201` and `reps=3`: portable `≈19.2s`, native `≈16.6s`, native+LTO `≈16.7s`.
- **KK mean-difference IVWC (`InferenceAllKKMeanDiffIVWC`)**: portable ≈13.4s, native ≈13.1s, native+LTO ≈13.8s.
- **Proportion fractional logit & simple Poisson (`InferencePropMultiFractionalLogit`, `InferenceCountUnivPoissonRegr`)**: portable was slightly faster than both native and native+LTO on those lightweight cases.

Bottom line: use `EDI_NATIVE_SPEED`/`EDI_NATIVE_LTO` to benchmark and tune the expensive ordinal/KK regressions locally, but keep the default portable flags for general development/distribution.

### Bootstrap diagnostics for modified Poisson incidence inference
Trimmed versions of the `cars`/`FixedCluster` workload (e.g., `scripts/diagnose_modified_poisson_bootstrap.R`) used to hit `Bootstrap confidence interval returned NA bounds` because the reduced design matrix had as many covariates as rows. The inference object now falls back to the univariate modified Poisson fit whenever the multivariate design is underdetermined, so the bootstrap diagnostics script now reports 25/25 finite replicates and `prop_illegal_values = 0.000` while still sharing the same treatment estimate. Run that script to reproduce the failure mode and confirm the fallback path locally.

To cite please use

Kapelner, A., & Krieger, A. (2023). A Matching Procedure for Sequential Experiments that Iteratively Learns which Covariates Improve Power. Biometrics. https://doi.org/10.1111/biom.13561

Kapelner, A. (2026). A Sequential Experiment Design for General Response Type and Unequal Allocation that Iteratively Learns Important Covariates. working paper
