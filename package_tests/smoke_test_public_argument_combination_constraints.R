#!/usr/bin/env Rscript

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!length(script_path) || is.na(script_path)) script_path = "package_tests/smoke_test_public_argument_combination_constraints.R"
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "public_argument_combination_constraints.R"))

expect_true = function(x, msg) {
	if (!isTRUE(x)) stop(msg, call. = FALSE)
}

expect_status = function(case, expected_status) {
	res = constraint_context_result(
		args = case$args,
		fixture = case$fixture,
		tier = case$tier,
		registry_entry = case$registry_entry
	)
	if (!identical(res$status, expected_status)) {
		stop(sprintf("Expected status %s, got %s: %s", expected_status, res$status, res$reason), call. = FALSE)
	}
	invisible(res)
}

fixture = default_constraint_fixture()
base_case = list(args = list(), fixture = fixture, tier = "ci", registry_entry = list())
schema_res = expect_status(base_case, "valid")$results
expect_true(all(c("constraint_name", "is_valid", "status", "reason", "source") %in% names(schema_res)), "constraint result schema is incomplete")
expect_true(all(schema_res$status %in% valid_constraint_statuses), "returned status is outside the allowed vocabulary")
expect_true(all(c(
	"constraint_min_usable_samples_le_B",
	"constraint_bca_requires_jackknife",
	"constraint_studentized_requires_finite_se",
	"constraint_use_rcpp_false_optional_backend",
	"constraint_show_progress_true_tier",
	"constraint_transform_responses",
	"constraint_formula_columns",
	"constraint_strata_cols_valid",
	"constraint_cluster_col_valid",
	"constraint_strata_cluster_no_conflict",
	"constraint_fixed_blocking_n_compatible",
	"constraint_survival_censoring_supported",
	"constraint_optional_dependency_paths",
	"constraint_known_slow_paths"
) %in% schema_res$constraint_name), "not every required constraint was evaluated")

cases = example_constraint_cases()
expect_status(cases$valid, "valid")
expect_status(cases$invalid_registry, "invalid_registry")
expect_status(cases$unsupported, "unsupported")
expect_status(cases$skipped_dependency, "skipped_dependency")
expect_status(cases$skipped_slow, "skipped_slow")

expect_status(list(args = list(type = "bca"), fixture = within(fixture, n_exchangeable_units <- 2L), tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(type = "studentized"), fixture = within(fixture, has_finite_standard_errors <- FALSE), tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(type = "studentized"), fixture = within(within(fixture, has_finite_standard_errors <- FALSE), documented_nonestimable_behavior <- TRUE), tier = "ci", registry_entry = list()), "valid")
expect_status(list(args = list(transform_responses = "sqrt"), fixture = fixture, tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(model_formula = "~ x1 + missing_col"), fixture = fixture, tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(strata_cols = "missing_col"), fixture = fixture, tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(strata_cols = "x1"), fixture = fixture, tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(cluster_col = "x1"), fixture = fixture, tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(strata_cols = "stratum", cluster_col = "stratum"), fixture = within(fixture, cluster_columns <- c("stratum")), tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(n = 10L, block_size = 4L), fixture = fixture, tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(), fixture = within(within(fixture, response_type <- "survival"), has_censoring <- TRUE), tier = "ci", registry_entry = list()), "unsupported")
expect_status(list(args = list(), fixture = fixture, tier = "ci", registry_entry = list(optional_dependencies = "optional_non_rcpp_backend")), "skipped_dependency")
expect_status(list(args = list(B = 101L), fixture = fixture, tier = "smoke", registry_entry = list(slow_path_labels = "bootstrap")), "skipped_slow")

message("public_argument_combination_constraints smoke test passed")
