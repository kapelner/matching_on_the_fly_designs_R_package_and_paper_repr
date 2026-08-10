# Public Argument Combination Quality Gates

The quality-gate checker enforces the argument-combination rollout without blocking broad coverage work too early.

Command:

```sh
Rscript package_tests/check_public_argument_combination_quality_gates.R report
Rscript package_tests/check_public_argument_combination_quality_gates.R ci
Rscript package_tests/check_public_argument_combination_quality_gates.R strict
```

Modes:

- `report`: writes gate CSVs and exits nonzero only for always-hard integrity failures.
- `ci`: same hard failures as `report`; this is intended for the bounded smoke/CI entry point.
- `strict`: also promotes future hard gates such as uncovered high-priority public APIs.

Outputs:

- `package_tests/public_argument_combination_quality_gates.csv`
- `package_tests/public_argument_combination_quality_gate_summary.csv`

Gate rows include `gate`, `severity`, `target`, `arg`, `value_expr`, `detail`, and `exempted`.

Implemented gates:

- `export_without_inventory`: exported `EDI` namespace symbols missing from the public API inventory.
- `multi_arg_without_combination_coverage`: public APIs with more than one configurable argument and no legal-combination case.
- `contract_without_registry_coverage`: public `checkmate` contracts with no registry target/argument row.
- `registry_value_never_selected`: registry values not selected by any valid generated case in the current artifact set.
- `unknown_constraint_status`: generated cases with an unrecognized constraint status.
- `ci_tier_unexpected_error`: smoke/CI result rows with unexpected errors or invariant failures.

Temporary exemptions can be added to `package_tests/public_argument_combination_gate_exemptions.csv` with exact `gate`, `target`, `arg`, and `value_expr` keys plus a reason. Empty exemption files are valid.
