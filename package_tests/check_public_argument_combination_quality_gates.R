#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/check_public_argument_combination_quality_gates.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

gate_artifact = function(name) file.path(repo_root, "package_tests", name)

gate_paths = list(
	inventory = gate_artifact("public_api_inventory.csv"),
	contracts = gate_artifact("checkmate_argument_contracts.csv"),
	registry = gate_artifact("public_argument_contract_registry.csv"),
	cases = gate_artifact("public_argument_combination_cases.csv"),
	results = gate_artifact("public_argument_combination_results.csv"),
	coverage = gate_artifact("public_argument_combination_coverage.csv"),
	ci_summary = gate_artifact("public_argument_combination_ci_failure_summary.csv"),
	exemptions = gate_artifact("public_argument_combination_gate_exemptions.csv"),
	report = gate_artifact("public_argument_combination_quality_gates.csv"),
	summary = gate_artifact("public_argument_combination_quality_gate_summary.csv")
)

allowed_constraint_statuses = c("valid", "unsupported", "skipped_slow", "skipped_dependency", "invalid_registry")

required_schemas = list(
	inventory = c("export_name", "api_kind", "class_name", "method_name", "configurable_arg_count", "has_more_than_one_configurable_arg"),
	contracts = c("api_name", "api_kind", "class_name", "method_name", "arg", "assertion", "coverage_scope"),
	registry = c("target", "api_name", "api_kind", "class_name", "method_name", "arg", "value_expr", "tier"),
	cases = c("case_id", "target", "api_kind", "class_name", "function_name", "method_name", "fixture_id", "tier", "argument_signature", "constraint_status"),
	results = c("case_id", "tier", "status", "invariant_status", "error_message"),
	coverage = c("target", "api_kind", "class_name", "method_name", "has_more_than_one_configurable_arg", "has_legal_combination_case", "has_executed_legal_combination")
)

read_required_csv = function(path, label) {
	if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

read_optional_csv = function(path, cols) {
	if (!file.exists(path)) {
		out = as.data.frame(setNames(replicate(length(cols), character(), simplify = FALSE), cols), stringsAsFactors = FALSE)
		return(out)
	}
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

target_from_parts = function(api_name, method_name) {
	api_name = as.character(api_name %||% "")
	method_name = as.character(method_name %||% "")
	ifelse(nzchar(method_name), paste(api_name, method_name, sep = "::"), api_name)
}

ensure_schema = function(df, name) {
	missing = setdiff(required_schemas[[name]], names(df))
	if (!length(missing)) return(NULL)
	data.frame(
		gate = "schema_validity",
		severity = "hard",
		target = "",
		arg = "",
		value_expr = "",
		detail = paste("Missing columns in", name, ":", paste(missing, collapse = ", ")),
		stringsAsFactors = FALSE
	)
}

load_edi_exports = function() {
	local_check_lib = file.path(repo_root, "EDI.Rcheck")
	if (!requireNamespace("EDI", quietly = TRUE) && dir.exists(file.path(local_check_lib, "EDI"))) {
		.libPaths(c(normalizePath(local_check_lib, mustWork = TRUE), .libPaths()))
	}
	if (!requireNamespace("EDI", quietly = TRUE)) return(character())
	sort(getNamespaceExports("EDI"))
}

case_argument_values = function(cases) {
	parse_signature = function(signature) {
		if (is.na(signature) || !nzchar(signature) || identical(signature, "<default>")) {
			return(data.frame(arg = character(), value_expr = character(), stringsAsFactors = FALSE))
		}
		parts = strsplit(signature, ";", fixed = TRUE)[[1L]]
		parts = parts[nzchar(parts)]
		data.frame(
			arg = sub("=.*$", "", parts),
			value_expr = sub("^[^=]*=", "", parts),
			stringsAsFactors = FALSE
		)
	}
	rows = list()
	for (i in seq_len(nrow(cases))) {
		vals = parse_signature(cases$argument_signature[i])
		if (!nrow(vals)) next
		vals$target = cases$target[i]
		vals$case_id = cases$case_id[i]
		vals$constraint_status = cases$constraint_status[i]
		rows[[length(rows) + 1L]] = vals
	}
	if (!length(rows)) return(data.frame(target = character(), arg = character(), value_expr = character(), stringsAsFactors = FALSE))
	out = do.call(rbind, rows)
	row.names(out) = NULL
	out
}

empty_gate_rows = function() {
	data.frame(gate = character(), severity = character(), target = character(), arg = character(), value_expr = character(), detail = character(), stringsAsFactors = FALSE)
}

gate_rows = function(gate, severity, df, detail) {
	if (!nrow(df)) return(empty_gate_rows())
	if (!"target" %in% names(df)) df$target = ""
	if (!"arg" %in% names(df)) df$arg = ""
	if (!"value_expr" %in% names(df)) df$value_expr = ""
	data.frame(
		gate = gate,
		severity = severity,
		target = as.character(df$target),
		arg = as.character(df$arg),
		value_expr = as.character(df$value_expr),
		detail = detail,
		stringsAsFactors = FALSE
	)
}

apply_exemptions = function(rows, exemptions) {
	if (!nrow(rows) || !nrow(exemptions)) return(rows)
	for (col in c("gate", "target", "arg", "value_expr")) {
		if (!col %in% names(exemptions)) exemptions[[col]] = ""
	}
	exemptions$key = paste(exemptions$gate, exemptions$target, exemptions$arg, exemptions$value_expr, sep = "\r")
	rows$key = paste(rows$gate, rows$target, rows$arg, rows$value_expr, sep = "\r")
	rows$exempted = rows$key %in% exemptions$key
	rows$key = NULL
	rows
}

high_priority_target = function(target) {
	grepl("^Inference|^Design", target) &
		grepl("initialize|model_formula|bootstrap|bayesian_bootstrap|rand|asymp|wald|score|lik_ratio|gradient|jackknife|add_all_subject|assign_w|strata|cluster", target, ignore.case = TRUE)
}

build_quality_gate_report = function(inputs, exports = load_edi_exports()) {
	rows = list()
	for (name in names(required_schemas)) {
		rows[[length(rows) + 1L]] = ensure_schema(inputs[[name]], name)
	}
	if (length(exports)) {
		export_inventory = inputs$inventory[inputs$inventory$method_name == "" | is.na(inputs$inventory$method_name), , drop = FALSE]
		missing_exports = setdiff(exports, export_inventory$export_name)
		if (length(missing_exports)) {
			rows[[length(rows) + 1L]] = gate_rows(
				"export_without_inventory",
				"hard",
				data.frame(target = missing_exports, stringsAsFactors = FALSE),
				"Exported namespace symbol is absent from public API inventory."
			)
		}
	}
	uncovered = inputs$coverage[
		inputs$coverage$has_more_than_one_configurable_arg %in% c(TRUE, "TRUE") &
			!(inputs$coverage$has_legal_combination_case %in% c(TRUE, "TRUE")),
		,
		drop = FALSE
	]
	rows[[length(rows) + 1L]] = gate_rows(
		"multi_arg_without_combination_coverage",
		ifelse(high_priority_target(uncovered$target), "hard_later", "report"),
		uncovered,
		"Public API has multiple configurable arguments and no generated legal-combination case."
	)
	contracts = inputs$contracts
	contracts$target = target_from_parts(contracts$api_name, contracts$method_name)
	public_contracts = contracts[contracts$coverage_scope == "public_contract", , drop = FALSE]
	contract_keys = unique(public_contracts[, c("target", "arg"), drop = FALSE])
	registry_keys = unique(inputs$registry[, c("target", "arg"), drop = FALSE])
	missing_registry = contract_keys[!paste(contract_keys$target, contract_keys$arg, sep = "\r") %in% paste(registry_keys$target, registry_keys$arg, sep = "\r"), , drop = FALSE]
	rows[[length(rows) + 1L]] = gate_rows(
		"contract_without_registry_coverage",
		"report",
		missing_registry,
		"Checkmate-derived public argument contract has no registry target/argument coverage."
	)
	case_values = case_argument_values(inputs$cases)
	valid_case_values = case_values[case_values$constraint_status == "valid", , drop = FALSE]
	registry_value_keys = unique(inputs$registry[, c("target", "arg", "value_expr"), drop = FALSE])
	covered_value_key = paste(valid_case_values$target, valid_case_values$arg, valid_case_values$value_expr, sep = "\r")
	unselected_values = registry_value_keys[!paste(registry_value_keys$target, registry_value_keys$arg, registry_value_keys$value_expr, sep = "\r") %in% covered_value_key, , drop = FALSE]
	rows[[length(rows) + 1L]] = gate_rows(
		"registry_value_never_selected",
		"report",
		unselected_values,
		"Registry value is not selected by any valid generated case in the current tier/artifact set."
	)
	unknown_status = inputs$cases[!inputs$cases$constraint_status %in% allowed_constraint_statuses, , drop = FALSE]
	rows[[length(rows) + 1L]] = gate_rows(
		"unknown_constraint_status",
		"hard",
		data.frame(target = unknown_status$target, value_expr = unknown_status$case_id, stringsAsFactors = FALSE),
		paste("Generated case has constraint_status outside:", paste(allowed_constraint_statuses, collapse = ", "))
	)
	ci_errors = inputs$results[
		inputs$results$tier %in% c("smoke", "ci") &
			(inputs$results$status == "error" | grepl("^failed:", inputs$results$invariant_status)),
		,
		drop = FALSE
	]
	rows[[length(rows) + 1L]] = gate_rows(
		"ci_tier_unexpected_error",
		"hard",
		data.frame(target = ci_errors$case_id, value_expr = ci_errors$error_message, stringsAsFactors = FALSE),
		"Smoke/CI tier result has an unexpected error or invariant failure."
	)
	report = do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
	if (is.null(report) || !nrow(report)) report = empty_gate_rows()
	report
}

summarize_quality_gates = function(report, mode = "report") {
	active_hard = report$severity == "hard"
	if (identical(mode, "strict")) {
		active_hard = report$severity %in% c("hard", "hard_later")
	}
	if ("exempted" %in% names(report)) active_hard = active_hard & !report$exempted
	data.frame(
		mode = mode,
		total_rows = nrow(report),
		report_rows = sum(report$severity == "report"),
		future_hard_rows = sum(report$severity == "hard_later"),
		hard_rows = sum(report$severity == "hard"),
		active_hard_rows = sum(active_hard),
		ci_should_fail = sum(active_hard) > 0L,
		stringsAsFactors = FALSE
	)
}

check_public_argument_combination_quality_gates = function(mode = "report") {
	inputs = lapply(names(required_schemas), function(name) read_required_csv(gate_paths[[name]], name))
	names(inputs) = names(required_schemas)
	inputs$ci_summary = if (file.exists(gate_paths$ci_summary)) read_required_csv(gate_paths$ci_summary, "ci_summary") else data.frame()
	exemptions = read_optional_csv(gate_paths$exemptions, c("gate", "target", "arg", "value_expr", "reason"))
	report = apply_exemptions(build_quality_gate_report(inputs), exemptions)
	if (!"exempted" %in% names(report)) report$exempted = FALSE
	summary = summarize_quality_gates(report, mode)
	write.csv(report, gate_paths$report, row.names = FALSE)
	write.csv(summary, gate_paths$summary, row.names = FALSE)
	invisible(list(report = report, summary = summary))
}

main = function() {
	args = commandArgs(trailingOnly = TRUE)
	mode = if (length(args) >= 1L && nzchar(args[1L])) args[1L] else "report"
	out = check_public_argument_combination_quality_gates(mode)
	message("Wrote quality gate rows: ", nrow(out$report))
	message("Active hard gate rows: ", out$summary$active_hard_rows[1L])
	if (isTRUE(out$summary$ci_should_fail[1L])) quit(status = 1L)
}

called_as_gate_script = function() {
	file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% ""
	identical(basename(sub("^--file=", "", file_arg)), "check_public_argument_combination_quality_gates.R")
}

if (called_as_gate_script()) {
	main()
}
