#!/usr/bin/env Rscript

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
if (!length(script_path) || is.na(script_path)) script_path = "package_tests/smoke_test_analyze_public_argument_combinations.R"
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "analyze_public_argument_combinations.R"))

expect_true = function(x, msg) {
	if (!isTRUE(x)) stop(msg, call. = FALSE)
}

out = analyze_public_argument_combinations()

expect_true(nrow(out$coverage) > 0L, "Coverage report is empty.")
expect_true(all(c("target", "total_cases", "valid_candidate_cases", "executed_ok_cases", "covered_argument_values", "covered_argument_pairs") %in% names(out$coverage)), "Coverage report schema is incomplete.")
expect_true(any(out$coverage$has_legal_combination_case), "Analyzer found no target with legal combination cases.")
expect_true(nrow(out$uncovered) >= 0L && all(c("target", "has_checkmate_contract", "has_legal_combination_case") %in% names(out$uncovered)), "Uncovered API report is malformed.")
expect_true(all(c("drift_type", "target", "arg", "detail") %in% names(out$drift)), "Drift report schema is incomplete.")
expect_true(all(c("total_unexpected_errors", "total_invariant_failures", "ci_should_fail") %in% names(out$ci_summary)), "CI failure summary schema is incomplete.")
expect_true(!isTRUE(out$ci_summary$ci_should_fail[1L]), "Smoke analyzer found CI-blocking failures.")
expect_true(file.exists(paths$coverage), "Coverage CSV was not written.")
expect_true(file.exists(paths$failures), "Failures CSV was not written.")
expect_true(file.exists(paths$drift), "Drift CSV was not written.")
expect_true(file.exists(paths$slow), "Slow-case CSV was not written.")
expect_true(file.exists(paths$report), "HTML report was not written.")

message("public argument combination analyzer smoke test passed with ", nrow(out$coverage), " coverage rows")
