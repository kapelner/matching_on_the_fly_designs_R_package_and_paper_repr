#!/usr/bin/env Rscript

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!length(script_path) || is.na(script_path)) script_path = "package_tests/smoke_test_generate_public_argument_combinations.R"
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "generate_public_argument_combinations.R"))

expect_true = function(x, msg) {
	if (!isTRUE(x)) stop(msg, call. = FALSE)
}

cases = generate_public_argument_combinations(
	tier = "smoke",
	target_filter = "^InferenceAllSimpleMeanDiff::compute_bootstrap_confidence_interval$|^DesignFixedBernoulli::add_all_subject_responses$",
	fixture_ids = c("fixed_bernoulli_continuous_smoke", "fixed_bernoulli_survival_censored_smoke")
)

required_cols = c(
	"case_id", "api_kind", "class_name", "function_name", "method_name",
	"fixture_id", "tier", "argument_signature", "argument_coverage_kind",
	"source_contracts", "constraint_status", "constraint_reason"
)
expect_true(all(required_cols %in% names(cases)), "Generated cases are missing required schema columns.")
expect_true(nrow(cases) > 0L, "No cases generated.")
expect_true(!anyDuplicated(cases$case_id), "case_id values are not unique.")
expect_true(any(cases$argument_coverage_kind == "default"), "Default cases are missing.")
expect_true(any(cases$argument_coverage_kind == "one_non_default"), "One-non-default cases are missing.")
expect_true(any(cases$argument_coverage_kind == "pairwise"), "Pairwise cases are missing.")
expect_true(any(cases$constraint_status != "valid"), "Rejected candidates are missing from generated cases.")
expect_true(all(cases$constraint_status %in% valid_constraint_statuses), "Generated cases contain invalid constraint statuses.")

coverage = coverage_from_cases(cases)
expect_true(nrow(coverage) > 0L, "Coverage table is empty.")
expect_true(all(coverage$valid_cases <= coverage$total_cases), "Coverage valid counts exceed total counts.")

tmp_dir = file.path(tempdir(), paste0("public_arg_combo_smoke_", Sys.getpid()))
dir.create(tmp_dir, showWarnings = FALSE)
write_generation_outputs(
	cases,
	cases_file = file.path(tmp_dir, "cases.csv"),
	rejected_file = file.path(tmp_dir, "rejected.csv"),
	coverage_file = file.path(tmp_dir, "coverage.csv")
)
expect_true(file.exists(file.path(tmp_dir, "cases.csv")), "Cases CSV was not written.")
expect_true(file.exists(file.path(tmp_dir, "rejected.csv")), "Rejected CSV was not written.")
expect_true(file.exists(file.path(tmp_dir, "coverage.csv")), "Coverage CSV was not written.")

cases_again = generate_public_argument_combinations(
	tier = "smoke",
	target_filter = "^InferenceAllSimpleMeanDiff::compute_bootstrap_confidence_interval$|^DesignFixedBernoulli::add_all_subject_responses$",
	fixture_ids = c("fixed_bernoulli_continuous_smoke", "fixed_bernoulli_survival_censored_smoke")
)
expect_true(identical(cases, cases_again), "Generation is not deterministic.")

message("public argument combination generator smoke test passed with ", nrow(cases), " candidates")
