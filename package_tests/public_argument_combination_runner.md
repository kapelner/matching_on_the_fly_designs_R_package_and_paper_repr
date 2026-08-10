# Public Argument Combination Runner

The public argument-combination runner complements `package_tests/comprehensive_tests.R`.

Use `comprehensive_tests.R` for end-to-end statistical workflows across datasets, response families, designs, and inference classes. Its durable result keys are `dataset`, `response_type`, `design`, `inference_class`, and `function_run`.

Use the public argument-combination runner when the question is whether legal public API argument values also work together. It generates fixture-backed calls from the public API inventory, extracted `checkmate` contracts, the reviewed argument registry, and cross-argument constraints.

The two result families intentionally stay separate:

- Statistical comprehensive results keep the `comprehensive_tests_results_*` naming pattern.
- Argument-combination results keep the `public_argument_combination_*` naming pattern.
- Integrated argument-combination files add join-compatible labels without replacing statistical results.

Primary commands:

```sh
Rscript package_tests/generate_public_argument_combinations.R
Rscript package_tests/run_public_argument_combinations.R
Rscript package_tests/analyze_public_argument_combinations.R
Rscript package_tests/public_argument_combination_integration.R
Rscript package_tests/check_public_argument_combination_quality_gates.R report
```

Entry-point commands:

```sh
Rscript package_tests/run_public_argument_combination_ci.R
Rscript package_tests/run_public_argument_combination_nightly.R
```

When `EDI` is not installed in the active R library, the entry points fall back to `EDI.Rcheck/EDI` if that checked package directory exists.

The CI entry point runs the bounded smoke target set and writes the standard argument-combination artifacts. The nightly entry point defaults to the same bounded target set; pass `ALL` as the first argument to generate from every registry target, or pass a regular expression to select a target subset.

Integrated artifacts:

- `package_tests/public_argument_combination_cases_integrated.csv`
- `package_tests/public_argument_combination_results_integrated.csv`
- `package_tests/public_argument_combination_coverage_integrated.csv`

Integrated rows include `response_type`, `design_type`, `design`, `inference_class`, `function_run`, `dataset_name`, and `dataset` where fixture or target metadata can provide them. `dataset_name` and `dataset` use the deterministic fixture id because argument-combination cases are fixture-backed rather than dataset-backed.

To promote a failure into a focused `testthat` regression test:

1. Open `package_tests/public_argument_combination_failures.csv`.
2. Copy the failing `target`, `fixture_id`, `argument_signature`, `status`, `invariant_status`, and `error_message`.
3. Build the same fixture with `build_public_argument_fixtures(fixture_ids = "<fixture_id>")`.
4. Recreate the public call from `argument_signature`.
5. Assert the documented behavior in the relevant `EDI/tests/testthat/` file.
