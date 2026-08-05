# edi_kernels

Python bindings for EDI's C++ model-fitting kernels, via pybind11.

This package has **no R or Rcpp dependency**. It compiles the same
`EDI/src/*.cpp` model-fitting kernels the R package
[`EDI`](https://github.com/kapelner/matching_on_the_fly_designs_R_package_and_paper_repr/tree/main/EDI)
(not yet on CRAN) uses, built under an `EDI_CORE_ONLY` preprocessor guard
that swaps out `RcppEigen`/`Rmath` for a vanilla `Eigen` + `LBFGSpp`, both
fetched directly from their own upstream repositories (see
`CMakeLists.txt`). Nothing under `python/` is a copy of a `EDI/src/*.cpp`
or `*.h` file — the compiled extension `#include`s them directly from the
R package's own source tree.

See `../package_metadata/python_bindings_package_spec.md` for the full
design spec, kernel-by-kernel scope, and baseline-benchmarking methodology.

See `../package_metadata/benchmark_model_fits_python.html` for the
generated benchmark report comparing each bound kernel's speed and
correctness against its Python canonical baseline (statsmodels, lifelines,
scikit-survival, etc., where one exists), in the same column shape as the
R package's own `benchmark_model_fits_R.html`. Kernels without a clean
canonical baseline are marked as explicit Baseline Gaps rather than
omitted.

## Install

```bash
pip install .
```

## Usage

```python
import numpy as np
from edi_kernels import fast_ols

X = np.column_stack([np.ones(100), np.random.default_rng(0).normal(size=100)])
y = X @ [1.0, 0.5] + np.random.default_rng(1).normal(size=100)
result = fast_ols(X, y)
print(result["b"])
```

Every bound kernel returns a `dict` (see `python/cpp/result_map_pybind.h`)
mirroring the R package's own `edi::ResultMap` -> list conversion field for
field.
