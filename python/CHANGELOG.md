# Changelog

All notable changes to `edi_kernels` are documented here. The version
number tracks `EDI/DESCRIPTION`'s `Version` field (see
`package_metadata/python_bindings_package_spec.md`'s "Versioning"
checklist item) — a `.postN` suffix is used for Python-packaging-only
changes that don't touch `EDI/src/*.cpp`.

## [1.0.0] - 2026-08-05

Initial release.

### Added

- pybind11 bindings for all 33 `EDI/src` model-fitting kernels (37 bound
  Python functions total, counting the `_with_var` secondary entry points
  for `fast_log_binomial_regression`/`fast_identity_binomial_regression`),
  covering the continuous, binary, count, proportion, ordinal, incidence/
  GEE, survival, and GLMM/CLMM/LMM response-type families.
- An `EDI_CORE_ONLY` build path: every bound kernel compiles directly out
  of `EDI/src/*.cpp` with zero R/Rcpp/Rmath dependency, against vanilla
  Eigen + LBFGSpp fetched from their own upstream repositories. No kernel
  logic is duplicated or reimplemented in Python or in the `pybind11`
  binding layer.
- Every `Rcpp::Nullable<T>` argument exposed as a Python keyword argument
  with a working `std::optional`-based default — full optionality, not
  just the arguments a typical user needs.
- `python/tests/`: an R-fixture parity test for every bound kernel
  (`atol=1e-9, rtol=1e-9`) plus a dedicated omitted-argument test per
  kernel proving every optional parameter has a working default when left
  unset. 176 tests, all passing.
- `benchmarks/baselines.py`: a canonical-baseline registry (statsmodels,
  lifelines, scikit-survival, etc.) for the response-type families that
  have a clean Python equivalent, plus explicit Baseline Gap documentation
  for the families that don't (GLMM/CLMM/LMM, adjacent-category/
  continuation-ratio/stereotype ordinal, zero-one-inflated beta, Weibull
  frailty, KK-combined estimators) rather than silently omitting them.
- `package_metadata/benchmark_model_fits_python.html`: a generated
  benchmark report in the same column shape as the R package's own
  `benchmark_model_fits_R.html`.
- Packaging: `python/README.md`, `python/LICENSE` (GPL-3.0-only, matching
  `EDI/DESCRIPTION`'s `License: GPL-3`), portable-wheel support via
  `[tool.cibuildwheel]` (manylinux_2_28, macOS x86_64+arm64, Windows with a
  `scipy-openblas32` BLAS fallback), and a CI pipeline
  (`.github/workflows/build-wheels.yml`) that builds every platform's wheel
  plus an sdist and publishes to PyPI via Trusted Publishing (OIDC) on a
  `py-v*` tag push.

### Known issues

- `fast_zero_one_inflated_beta`'s R-facing return list has no
  `converged`/`iterations`/`gradient_norm` fields, even though the
  underlying `LikelihoodFitResult` the Python binding surfaces does carry
  them — an R-side gap, not a Python binding bug.
- `estimate_only=True` omits the `vcov`/`params` keys entirely on some
  kernels (`fast_gaussian_lmm`, `fast_coxph_regression`,
  `fast_weibull_regression`) rather than setting them to `None` the way
  `fast_poisson_glmm` does — inconsistent across kernels, documented
  per-file in the affected parity tests.
