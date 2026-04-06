# SeqExpMatch

Code for the R package `SeqExpMatch` and the code to reproduce the simulations in the papers listed below. You can get started by cloning and running `R CMD INSTALL SeqExpMatch` from the root of this repo. Then try running the code in the folder `package_tests`.

For a local performance build, you can opt into native compiler flags without
changing the default portable build:

```sh
EDI_NATIVE_SPEED=1 R CMD INSTALL SeqExpMatch
```

This enables `-O3 -march=native -mtune=native` when your compiler supports
those flags. Use it for local benchmarking only; the default build remains the
portable one.

To compare plain native vs native plus link-time optimization:

```sh
EDI_NATIVE_SPEED=1 R CMD INSTALL SeqExpMatch
EDI_NATIVE_SPEED=1 EDI_NATIVE_LTO=1 R CMD INSTALL SeqExpMatch
```

`EDI_NATIVE_LTO=1` adds `-flto` on top of the native-speed flags.

To benchmark the three builds back-to-back, run:

```sh
bash scripts/benchmark_build_modes.sh
```

To benchmark a randomization CI workload at `num_cores = 3` across the same
three build modes, run:

```sh
bash scripts/benchmark_randomization_ci_build_modes.sh
```

To cite please use

Kapelner, A., & Krieger, A. (2023). A Matching Procedure for Sequential Experiments that Iteratively Learns which Covariates Improve Power. Biometrics. https://doi.org/10.1111/biom.13561

Kapelner, A. (2026). A Sequential Experiment Design for General Response Type and Unequal Allocation that Iteratively Learns Important Covariates. working paper
