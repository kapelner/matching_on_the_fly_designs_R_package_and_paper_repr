#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/smoke_test_public_argument_combination_quality_gates.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "check_public_argument_combination_quality_gates.R"))

out = check_public_argument_combination_quality_gates("report")

stopifnot(nrow(out$report) > 0L)
stopifnot(all(c("gate", "severity", "target", "arg", "value_expr", "detail", "exempted") %in% names(out$report)))
stopifnot(all(c("mode", "total_rows", "active_hard_rows", "ci_should_fail") %in% names(out$summary)))
stopifnot(all(file.exists(unlist(gate_paths[c("report", "summary")]))))
stopifnot(!any(out$report$gate == "unknown_constraint_status"))

message("public argument combination quality gates smoke test passed")
