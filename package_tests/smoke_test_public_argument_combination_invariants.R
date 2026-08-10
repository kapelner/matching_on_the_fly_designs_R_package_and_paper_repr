#!/usr/bin/env Rscript

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!length(script_path) || is.na(script_path)) script_path = "package_tests/smoke_test_public_argument_combination_invariants.R"
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "public_argument_combination_invariants.R"))

expect_true = function(x, msg) {
	if (!isTRUE(x)) stop(msg, call. = FALSE)
}

all_names = names(public_argument_invariant_functions)
expected_names = c(
	"no_unexpected_error",
	"numeric_finite_or_documented_missing",
	"pval_in_0_1",
	"ci_length_2",
	"ci_ordered",
	"estimate_shape_valid",
	"distribution_length_valid",
	"documented_nonestimable_ok",
	"seeded_call_deterministic",
	"output_names_stable_when_documented"
)
expect_true(identical(sort(all_names), sort(expected_names)), "Invariant registry is missing required invariants.")

ok_results = evaluate_public_argument_invariants(
	output = c(lower = 0.1, upper = 0.2),
	execution_status = "ok",
	invariants_expr = "no_unexpected_error;numeric_finite_or_documented_missing;ci_length_2;ci_ordered;output_names_stable_when_documented",
	context = list(expected_names = c("lower", "upper"))
)
expect_true(all(ok_results$status %in% valid_invariant_statuses), "Invariant status vocabulary is invalid.")
expect_true(identical(summarize_public_argument_invariants(ok_results), "ok"), "Expected all invariant checks to pass.")

failed_pval = evaluate_public_argument_invariants(1.5, "ok", "pval_in_0_1")
expect_true(identical(summarize_public_argument_invariants(failed_pval), "failed:pval_in_0_1"), "Bad p-value was not detected.")

failed_ci = evaluate_public_argument_invariants(c(2, 1), "ok", "ci_length_2;ci_ordered")
expect_true(identical(summarize_public_argument_invariants(failed_ci), "failed:ci_ordered"), "Bad CI ordering was not detected.")

missing_ok = evaluate_public_argument_invariants(c(NA_real_, NA_real_), "ok", "numeric_finite_or_documented_missing;documented_nonestimable_ok", list(documented_nonestimable_behavior = TRUE))
expect_true(identical(summarize_public_argument_invariants(missing_ok), "ok"), "Documented missing output was not accepted.")

distribution_ok = evaluate_public_argument_invariants(1:5, "ok", "distribution_length_valid", list(expected_distribution_length = 5L))
expect_true(identical(summarize_public_argument_invariants(distribution_ok), "ok"), "Valid distribution length was rejected.")

estimate_ok = evaluate_public_argument_invariants(0.25, "ok", "estimate_shape_valid")
expect_true(identical(summarize_public_argument_invariants(estimate_ok), "ok"), "Valid estimate shape was rejected.")

deterministic_ok = evaluate_public_argument_invariants(1:3, "ok", "seeded_call_deterministic", list(repeated_output = 1:3))
expect_true(identical(summarize_public_argument_invariants(deterministic_ok), "ok"), "Deterministic seeded output was rejected.")

name_fail = evaluate_public_argument_invariants(c(a = 1), "ok", "output_names_stable_when_documented", list(expected_names = "b"))
expect_true(identical(summarize_public_argument_invariants(name_fail), "failed:output_names_stable_when_documented"), "Name drift was not detected.")

skipped = evaluate_public_argument_invariants(NULL, "error", "no_unexpected_error;ci_length_2")
expect_true(grepl("^failed:", summarize_public_argument_invariants(skipped)), "Unexpected error was not treated as invariant failure.")

message("public_argument_combination_invariants smoke test passed")
