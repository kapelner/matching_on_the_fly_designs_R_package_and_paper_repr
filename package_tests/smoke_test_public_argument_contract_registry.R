#!/usr/bin/env Rscript

rds_path = "package_tests/public_argument_contract_registry.rds"
csv_path = "package_tests/public_argument_contract_registry.csv"
if (!file.exists(rds_path) || !file.exists(csv_path)) {
	stop("Missing public argument contract registry artifacts. Run package_tests/public_argument_contract_registry.R first.", call. = FALSE)
}

registry = readRDS(rds_path)
flat = read.csv(csv_path, stringsAsFactors = FALSE)

if (!length(registry)) stop("Registry is empty.", call. = FALSE)
if (!nrow(flat)) stop("Flattened registry is empty.", call. = FALSE)

required_entry_fields = c("domains", "constraints", "invariants", "tier_metadata")
bad_entries = names(registry)[vapply(registry, function(entry) {
	any(!required_entry_fields %in% names(entry))
}, logical(1))]
if (length(bad_entries)) {
	stop("Registry entries missing required fields: ", paste(head(bad_entries), collapse = ", "), call. = FALSE)
}

required_sources = c("checkmate_choices", "checkmate_count", "checkmate_flag", "formal_default", "formula_rule", "runtime_tier")
missing_sources = setdiff(required_sources, unique(flat$source))
if (length(missing_sources)) {
	stop("Missing expected registry value sources: ", paste(missing_sources, collapse = ", "), call. = FALSE)
}

required_tiers = c("smoke", "ci", "nightly", "release")
missing_tiers = setdiff(required_tiers, unique(flat$tier))
if (length(missing_tiers)) {
	stop("Missing expected tiers: ", paste(missing_tiers, collapse = ", "), call. = FALSE)
}

high_priority_methods = c(
	"compute_asymp_two_sided_pval",
	"compute_asymp_confidence_interval",
	"compute_bootstrap_two_sided_pval",
	"compute_bootstrap_confidence_interval",
	"compute_bayesian_bootstrap_two_sided_pval",
	"compute_bayesian_bootstrap_confidence_interval",
	"compute_rand_two_sided_pval",
	"compute_rand_confidence_interval"
)
missing_hp = setdiff(high_priority_methods, unique(flat$method_name))
if (length(missing_hp)) {
	stop("Missing high-priority methods with argument contracts: ", paste(missing_hp, collapse = ", "), call. = FALSE)
}

if (!any(grepl("min_number_usable_samples <= B", flat$constraints, fixed = TRUE))) {
	stop("Expected min_number_usable_samples <= B combination constraint.", call. = FALSE)
}
if (!any(grepl("ci_length_2", flat$invariants, fixed = TRUE))) {
	stop("Expected confidence interval invariants.", call. = FALSE)
}
if (!any(grepl("pval_in_0_1", flat$invariants, fixed = TRUE))) {
	stop("Expected p-value invariants.", call. = FALSE)
}

message(
	"Phase 2 registry smoke test passed with ",
	length(registry), " entries and ", nrow(flat), " flattened value rows."
)
