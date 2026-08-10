#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/smoke_test_public_argument_combination_integration.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "public_argument_combination_integration.R"))

required_join_keys = c("response_type", "design_type", "inference_class", "function_run", "dataset_name")

out = integrate_public_argument_combination_outputs()

stopifnot(nrow(out$cases) > 0L)
stopifnot(all(required_join_keys %in% names(out$cases)))
stopifnot(all(required_join_keys %in% names(out$results)))
stopifnot(all(required_join_keys %in% names(out$coverage)))
stopifnot(any(nzchar(out$cases$response_type)))
stopifnot(any(nzchar(out$cases$design_type)))
stopifnot(any(nzchar(out$cases$function_run)))
stopifnot(all(file.exists(unlist(integration_paths[c("cases_integrated", "results_integrated", "coverage_integrated")]))))
stopifnot(file.exists(file.path(repo_root, "package_tests", "comprehensive_tests.R")))

message("public argument combination integration smoke test passed")
