#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/run_public_argument_combination_ci.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

local_check_lib = file.path(repo_root, "EDI.Rcheck")
if (!requireNamespace("EDI", quietly = TRUE) && dir.exists(file.path(local_check_lib, "EDI"))) {
	.libPaths(c(normalizePath(local_check_lib, mustWork = TRUE), .libPaths()))
}

source(file.path(repo_root, "package_tests", "generate_public_argument_combinations.R"))
source(file.path(repo_root, "package_tests", "run_public_argument_combinations.R"))
source(file.path(repo_root, "package_tests", "analyze_public_argument_combinations.R"))
source(file.path(repo_root, "package_tests", "public_argument_combination_integration.R"))
source(file.path(repo_root, "package_tests", "check_public_argument_combination_quality_gates.R"))

main = function() {
	args = commandArgs(trailingOnly = TRUE)
	tier = if (length(args) >= 1L && nzchar(args[1L])) args[1L] else "smoke"
	class_filter = if (length(args) >= 2L && nzchar(args[2L])) args[2L] else NULL
	method_filter = if (length(args) >= 3L && nzchar(args[3L])) args[3L] else NULL
	timeout_sec = if (length(args) >= 4L && nzchar(args[4L])) as.numeric(args[4L]) else 20
	cases = generate_public_argument_combinations(
		tier = tier,
		target_filter = default_smoke_target_filter,
		fixture_ids = default_smoke_fixture_ids
	)
	write_generation_outputs(cases)
	results = run_public_argument_combinations(
		tier = tier,
		class_filter = class_filter,
		method_filter = method_filter,
		timeout_sec = timeout_sec
	)
	analysis = analyze_public_argument_combinations()
	integrated = integrate_public_argument_combination_outputs()
	gates = check_public_argument_combination_quality_gates("ci")
	message("Public argument combination CI rows: cases=", nrow(cases), ", results=", nrow(results), ", coverage=", nrow(analysis$coverage), ", integrated_cases=", nrow(integrated$cases), ", active_hard_gates=", gates$summary$active_hard_rows[1L])
	if (isTRUE(analysis$ci_summary$ci_should_fail[1L]) || isTRUE(gates$summary$ci_should_fail[1L])) quit(status = 1L)
}

main()
