#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/run_public_argument_combinations.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "generate_public_argument_combinations.R"))
source(file.path(repo_root, "package_tests", "public_argument_combination_invariants.R"))

cases_path = file.path(repo_root, "package_tests", "public_argument_combination_cases.csv")
results_path = file.path(repo_root, "package_tests", "public_argument_combination_results.csv")
failures_path = file.path(repo_root, "package_tests", "public_argument_combination_failures.csv")

runner_statuses = c("ok", "error", "unsupported", "skipped_slow", "skipped_dependency", "invalid_registry")

parse_runner_args = function(args = commandArgs(trailingOnly = TRUE)) {
	list(
		tier = if (length(args) >= 1L && !is.na(args[1L]) && nzchar(args[1L])) args[1L] else "smoke",
		class_filter = if (length(args) >= 2L && !is.na(args[2L]) && nzchar(args[2L])) args[2L] else NULL,
		method_filter = if (length(args) >= 3L && !is.na(args[3L]) && nzchar(args[3L])) args[3L] else NULL,
		timeout_sec = if (length(args) >= 4L && !is.na(args[4L]) && nzchar(args[4L])) as.numeric(args[4L]) else 20
	)
}

load_cases = function(path = cases_path) {
	if (!file.exists(path)) {
		stop("Missing package_tests/public_argument_combination_cases.csv. Run Phase 5 first.", call. = FALSE)
	}
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

filter_cases = function(cases, tier = "smoke", class_filter = NULL, method_filter = NULL) {
	cases = cases[cases$tier == tier, , drop = FALSE]
	if (!is.null(class_filter)) {
		cases = cases[grepl(class_filter, cases$class_name) | grepl(class_filter, cases$function_name) | grepl(class_filter, cases$target), , drop = FALSE]
	}
	if (!is.null(method_filter)) {
		cases = cases[grepl(method_filter, cases$method_name) | grepl(method_filter, cases$function_name) | grepl(method_filter, cases$target), , drop = FALSE]
	}
	row.names(cases) = NULL
	cases
}

read_existing_results = function(path = results_path) {
	if (!file.exists(path)) return(NULL)
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

split_semicolon_signature = function(signature) {
	if (!nzchar(signature) || identical(signature, "<default>")) return(character())
	parts = strsplit(signature, ";", fixed = TRUE)[[1L]]
	parts[nzchar(parts)]
}

args_from_signature = function(signature, fixture) {
	parts = split_semicolon_signature(signature)
	if (!length(parts)) return(list())
	out = list()
	for (part in parts) {
		key = sub("=.*$", "", part)
		value_expr = sub("^[^=]*=", "", part)
		out[[key]] = eval_value_expr(value_expr, fixture)
	}
	out
}

output_type = function(x) {
	if (is.null(x)) return("NULL")
	paste(class(x), collapse = "/")
}

output_shape = function(x) {
	if (is.null(x)) return("NULL")
	dims = dim(x)
	if (!is.null(dims)) return(paste(dims, collapse = "x"))
	if (length(x) == 1L) return("scalar")
	paste0("length=", length(x))
}

capture_call = function(fun, args, timeout_sec) {
	warnings = character()
	messages = character()
	start = proc.time()[["elapsed"]]
	old_limit = getOption("expressions")
	result = NULL
	error_message = ""
	status = "ok"
	tryCatch({
		setTimeLimit(elapsed = timeout_sec, transient = TRUE)
		on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
		result <<- withCallingHandlers(
			do.call(fun, args),
			warning = function(w) {
				warnings <<- c(warnings, conditionMessage(w))
				invokeRestart("muffleWarning")
			},
			message = function(m) {
				messages <<- c(messages, conditionMessage(m))
				invokeRestart("muffleMessage")
			}
		)
	}, error = function(e) {
		status <<- "error"
		error_message <<- conditionMessage(e)
	})
	duration = proc.time()[["elapsed"]] - start
	options(expressions = old_limit)
	list(
		status = status,
		result = result,
		duration_time_sec = duration,
		warning_message = paste(unique(warnings), collapse = " | "),
		message_text = paste(unique(messages), collapse = " | "),
		error_message = error_message
	)
}

construct_r6_receiver = function(case_row, fixture) {
	if (case_row$class_name %in% class(fixture$design) && case_row$method_name %in% names(fixture$design)) {
		return(fixture$design)
	}
	constructor = tryCatch(getExportedValue("EDI", case_row$class_name), error = function(e) NULL)
	if (is.null(constructor) || is.null(constructor$new)) {
		stop("Cannot resolve exported R6 class: ", case_row$class_name, call. = FALSE)
	}
	obj = tryCatch(constructor$new(des_obj = fixture$design), error = function(e) NULL)
	if (is.null(obj)) obj = tryCatch(constructor$new(fixture$design), error = function(e) NULL)
	if (is.null(obj)) stop("Cannot construct receiver for public R6 class: ", case_row$class_name, call. = FALSE)
	if (!case_row$method_name %in% names(obj)) {
		stop("Public method not found on receiver: ", case_row$method_name, call. = FALSE)
	}
	obj
}

execute_valid_case = function(case_row, fixture, timeout_sec) {
	args = args_from_signature(case_row$argument_signature, fixture)
	if (identical(case_row$api_kind, "r6_public_method")) {
		receiver = construct_r6_receiver(case_row, fixture)
		return(capture_call(receiver[[case_row$method_name]], args, timeout_sec))
	}
	if (identical(case_row$api_kind, "function")) {
		fun = getExportedValue("EDI", case_row$function_name)
		return(capture_call(fun, args, timeout_sec))
	}
	list(
		status = "error",
		result = NULL,
		duration_time_sec = 0,
		warning_message = "",
		message_text = "",
		error_message = paste("Unsupported api_kind for execution:", case_row$api_kind)
	)
}

result_row = function(case_row, status, duration_time_sec = 0, warning_message = "", message_text = "", output = NULL, invariant_status = "skipped", error_message = "") {
	clean = function(x) {
		x = as.character(x %||% "")
		ifelse(is.na(x), "", x)
	}
	data.frame(
		case_id = clean(case_row$case_id),
		api_kind = clean(case_row$api_kind),
		class_name = clean(case_row$class_name),
		function_name = clean(case_row$function_name),
		method_name = clean(case_row$method_name),
		fixture_id = clean(case_row$fixture_id),
		argument_signature = clean(case_row$argument_signature),
		tier = clean(case_row$tier),
		status = status,
		duration_time_sec = as.numeric(duration_time_sec),
		warning_message = warning_message,
		message_text = message_text,
		output_type = output_type(output),
		output_shape = output_shape(output),
		invariant_status = invariant_status,
		error_message = error_message,
		stringsAsFactors = FALSE
	)
}

run_case = function(case_row, fixtures, timeout_sec = 20) {
	if (!identical(case_row$constraint_status, "valid")) {
		return(result_row(
			case_row,
			status = case_row$constraint_status,
			error_message = case_row$constraint_reason
		))
	}
	fixture = fixtures[[case_row$fixture_id]]
	if (is.null(fixture)) {
		return(result_row(case_row, status = "error", error_message = paste("Missing fixture:", case_row$fixture_id)))
	}
	executed = execute_valid_case(case_row, fixture, timeout_sec)
	invariant_status = check_invariants(executed$result, executed$status, case_row$invariants)
	result_row(
		case_row,
		status = executed$status,
		duration_time_sec = executed$duration_time_sec,
		warning_message = executed$warning_message,
		message_text = executed$message_text,
		output = executed$result,
		invariant_status = invariant_status,
		error_message = executed$error_message
	)
}

run_public_argument_combinations = function(tier = "smoke", class_filter = NULL, method_filter = NULL, timeout_sec = 20, cases_file = cases_path, results_file = results_path, failures_file = failures_path) {
	load_edi_for_fixtures()
	cases = filter_cases(load_cases(cases_file), tier, class_filter, method_filter)
	if (!nrow(cases)) stop("No generated cases selected.", call. = FALSE)
	existing = read_existing_results(results_file)
	completed = if (is.null(existing)) character() else existing$case_id
	to_run = cases[!cases$case_id %in% completed, , drop = FALSE]
	fixture_ids = unique(to_run$fixture_id)
	fixtures = if (length(fixture_ids)) build_public_argument_fixtures(tier = tier, fixture_ids = fixture_ids) else list()
	new_results = lapply(seq_len(nrow(to_run)), function(i) run_case(to_run[i, , drop = FALSE], fixtures, timeout_sec))
	new_results = if (length(new_results)) do.call(rbind, new_results) else NULL
	results = if (!is.null(existing) && !is.null(new_results)) rbind(existing, new_results) else existing %||% new_results
	results = results[match(unique(results$case_id), results$case_id), , drop = FALSE]
	row.names(results) = NULL
	write.csv(results, results_file, row.names = FALSE)
	failures = results[results$status == "error" | grepl("^failed:", results$invariant_status), , drop = FALSE]
	write.csv(failures, failures_file, row.names = FALSE)
	invisible(results)
}

main = function() {
	parsed = parse_runner_args()
	results = run_public_argument_combinations(
		tier = parsed$tier,
		class_filter = parsed$class_filter,
		method_filter = parsed$method_filter,
		timeout_sec = parsed$timeout_sec
	)
	message("Wrote ", nrow(results), " public argument combination result rows.")
	message("Unexpected errors: ", sum(results$status == "error"))
	if (parsed$tier %in% c("ci", "smoke") && any(results$status == "error")) {
		quit(status = 1L)
	}
}

called_as_runner_script = function() {
	file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% ""
	identical(basename(sub("^--file=", "", file_arg)), "run_public_argument_combinations.R")
}

if (called_as_runner_script()) {
	main()
}
