#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/analyze_public_argument_combinations.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

artifact = function(name) file.path(repo_root, "package_tests", name)

paths = list(
	inventory = artifact("public_api_inventory.csv"),
	contracts = artifact("checkmate_argument_contracts.csv"),
	registry = artifact("public_argument_contract_registry.csv"),
	cases = artifact("public_argument_combination_cases.csv"),
	results = artifact("public_argument_combination_results.csv"),
	coverage = artifact("public_argument_combination_coverage.csv"),
	failures = artifact("public_argument_combination_failures.csv"),
	uncovered = artifact("public_argument_combination_uncovered_apis.csv"),
	drift = artifact("public_argument_combination_registry_drift.csv"),
	slow = artifact("public_argument_combination_slowest_cases.csv"),
	ci_summary = artifact("public_argument_combination_ci_failure_summary.csv"),
	report = artifact("public_argument_combination_report.html")
)

required_input_paths = c("inventory", "contracts", "registry", "cases", "results")

read_required_csv = function(path, label) {
	if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

target_from_parts = function(api_name, method_name) {
	api_name = as.character(api_name %||% "")
	method_name = as.character(method_name %||% "")
	ifelse(nzchar(method_name), paste(api_name, method_name, sep = "::"), api_name)
}

inventory_targets = function(inventory) {
	inventory$target = target_from_parts(inventory$export_name, inventory$method_name)
	inventory
}

contracts_targets = function(contracts) {
	contracts$target = target_from_parts(contracts$api_name, contracts$method_name)
	contracts
}

parse_signature = function(signature) {
	if (is.null(signature) || length(signature) == 0L || is.na(signature[1L]) || !nzchar(signature[1L]) || identical(signature[1L], "<default>")) {
		return(data.frame(arg = character(), value_expr = character(), stringsAsFactors = FALSE))
	}
	parts = strsplit(signature[1L], ";", fixed = TRUE)[[1L]]
	parts = parts[nzchar(parts)]
	data.frame(
		arg = sub("=.*$", "", parts),
		value_expr = sub("^[^=]*=", "", parts),
		stringsAsFactors = FALSE
	)
}

case_argument_values = function(cases) {
	rows = list()
	for (i in seq_len(nrow(cases))) {
		vals = parse_signature(cases$argument_signature[i])
		if (!nrow(vals)) next
		vals$case_id = cases$case_id[i]
		vals$target = cases$target[i]
		vals$constraint_status = cases$constraint_status[i]
		vals$argument_coverage_kind = cases$argument_coverage_kind[i]
		rows[[length(rows) + 1L]] = vals
	}
	if (!length(rows)) return(data.frame())
	out = do.call(rbind, rows)
	out = out[, c("target", "case_id", "arg", "value_expr", "argument_coverage_kind", "constraint_status"), drop = FALSE]
	row.names(out) = NULL
	out
}

case_pair_rows = function(cases) {
	rows = list()
	kinds = c("pairwise", "targeted_3way", "tiny_exhaustive")
	for (i in seq_len(nrow(cases))) {
		if (!cases$argument_coverage_kind[i] %in% kinds) next
		if (!identical(cases$constraint_status[i], "valid")) next
		vals = parse_signature(cases$argument_signature[i])
		if (nrow(vals) < 2L) next
		pairs = utils::combn(seq_len(nrow(vals)), 2L)
		for (j in seq_len(ncol(pairs))) {
			left = vals[pairs[1L, j], ]
			right = vals[pairs[2L, j], ]
			rows[[length(rows) + 1L]] = data.frame(
				target = cases$target[i],
				case_id = cases$case_id[i],
				arg_pair = paste(sort(c(left$arg, right$arg)), collapse = " + "),
				value_pair = paste(paste(c(left$arg, right$arg), c(left$value_expr, right$value_expr), sep = "="), collapse = " | "),
				coverage_kind = cases$argument_coverage_kind[i],
				stringsAsFactors = FALSE
			)
		}
	}
	if (!length(rows)) return(data.frame())
	out = unique(do.call(rbind, rows))
	row.names(out) = NULL
	out
}

count_by_target = function(df, target_col = "target") {
	if (!nrow(df)) return(data.frame(target = character(), n = integer(), stringsAsFactors = FALSE))
	tab = sort(table(df[[target_col]]), decreasing = FALSE)
	data.frame(target = names(tab), n = as.integer(tab), stringsAsFactors = FALSE)
}

value_coverage = function(cases, registry) {
	values = case_argument_values(cases)
	valid_values = values[values$constraint_status == "valid", , drop = FALSE]
	reg_keys = unique(paste(registry$target, registry$arg, registry$value_expr, sep = "\r"))
	covered_keys = unique(paste(valid_values$target, valid_values$arg, valid_values$value_expr, sep = "\r"))
	data.frame(
		total_registry_values = length(reg_keys),
		covered_registry_values = sum(reg_keys %in% covered_keys),
		coverage_fraction = if (length(reg_keys)) sum(reg_keys %in% covered_keys) / length(reg_keys) else NA_real_,
		stringsAsFactors = FALSE
	)
}

build_coverage_report = function(inventory, contracts, registry, cases, results) {
	targets = sort(unique(c(inventory$target, registry$target, cases$target)))
	arg_values = case_argument_values(cases)
	pairs = case_pair_rows(cases)
	first_rows = function(df, cols) {
		if (!nrow(df)) return(data.frame(target = character(), stringsAsFactors = FALSE))
		df = df[!duplicated(df$target), c("target", intersect(cols, names(df))), drop = FALSE]
		row.names(df) = NULL
		df
	}
	count_unique = function(df, key_cols, name) {
		if (!nrow(df)) return(data.frame(target = character(), stringsAsFactors = FALSE))
		keys = unique(df[, c("target", key_cols), drop = FALSE])
		out = aggregate(rep(1L, nrow(keys)), list(target = keys$target), sum)
		names(out)[2L] = name
		out
	}
	count_condition = function(df, condition, name) {
		if (!nrow(df)) return(data.frame(target = character(), stringsAsFactors = FALSE))
		tmp = df
		tmp$.hit = as.integer(condition)
		out = aggregate(tmp$.hit, list(target = tmp$target), sum)
		names(out)[2L] = name
		out
	}
	merge_left = function(base, extra) {
		if (!nrow(extra)) return(base)
		merge(base, extra, by = "target", all.x = TRUE, sort = FALSE)
	}
	fill_int = function(x) {
		if (is.null(x)) return(rep(0L, nrow(coverage)))
		x[is.na(x)] = 0L
		as.integer(x)
	}
	fill_chr = function(x) {
		x[is.na(x)] = ""
		as.character(x)
	}
	fill_lgl = function(x) {
		x[is.na(x)] = FALSE
		as.logical(x)
	}
	coverage = data.frame(target = targets, stringsAsFactors = FALSE)
	inv_meta = first_rows(inventory, c("api_kind", "class_name", "method_name"))
	reg_meta = first_rows(registry, c("api_kind", "class_name", "method_name"))
	case_meta = first_rows(cases, c("api_kind", "class_name", "method_name"))
	lookup_meta = function(df, col) {
		out = rep("", length(targets))
		if (!nrow(df) || !col %in% names(df)) return(out)
		idx = match(targets, df$target)
		hit = !is.na(idx)
		out[hit] = as.character(df[[col]][idx[hit]])
		out[is.na(out)] = ""
		out
	}
	pick_meta = function(primary, secondary, tertiary) {
		out = primary
		miss = is.na(out) | !nzchar(out)
		out[miss] = secondary[miss]
		miss = is.na(out) | !nzchar(out)
		out[miss] = tertiary[miss]
		fill_chr(out)
	}
	coverage$api_kind = pick_meta(lookup_meta(inv_meta, "api_kind"), lookup_meta(reg_meta, "api_kind"), lookup_meta(case_meta, "api_kind"))
	coverage$class_name = pick_meta(lookup_meta(inv_meta, "class_name"), lookup_meta(reg_meta, "class_name"), lookup_meta(case_meta, "class_name"))
	coverage$method_name = pick_meta(lookup_meta(inv_meta, "method_name"), lookup_meta(reg_meta, "method_name"), lookup_meta(case_meta, "method_name"))
	config = aggregate(as.integer(inventory$configurable_arg_count), list(target = inventory$target), max, na.rm = TRUE)
	names(config)[2L] = "configurable_arg_count"
	multi = aggregate(inventory$has_more_than_one_configurable_arg %in% c(TRUE, "TRUE"), list(target = inventory$target), any)
	names(multi)[2L] = "has_more_than_one_configurable_arg"
	coverage = merge_left(coverage, config)
	coverage = merge_left(coverage, multi)
	coverage$has_checkmate_contract = coverage$target %in% unique(contracts$target)
	coverage = merge_left(coverage, count_unique(registry, "arg", "registry_arg_count"))
	coverage = merge_left(coverage, count_unique(registry, c("arg", "value_expr"), "registry_value_count"))
	coverage = merge_left(coverage, count_condition(cases, rep(TRUE, nrow(cases)), "total_cases"))
	coverage = merge_left(coverage, count_condition(cases, cases$constraint_status == "valid", "valid_candidate_cases"))
	coverage = merge_left(coverage, count_condition(cases, cases$constraint_status == "valid" & cases$argument_coverage_kind %in% c("pairwise", "targeted_3way", "tiny_exhaustive"), "pairwise_or_higher_candidate_cases"))
	coverage = merge_left(coverage, count_condition(cases, cases$constraint_status == "valid" & cases$argument_coverage_kind == "targeted_3way", "targeted_3way_candidate_cases"))
	results_with_target = merge(results, cases[, c("case_id", "target"), drop = FALSE], by = "case_id", all.x = TRUE, sort = FALSE)
	coverage = merge_left(coverage, count_condition(results_with_target, results_with_target$status == "ok", "executed_ok_cases"))
	coverage = merge_left(coverage, count_condition(results_with_target, results_with_target$status == "error", "error_cases"))
	coverage = merge_left(coverage, count_condition(results_with_target, results_with_target$status == "unsupported", "unsupported_cases"))
	coverage = merge_left(coverage, count_condition(results_with_target, results_with_target$status == "skipped_slow", "skipped_slow_cases"))
	coverage = merge_left(coverage, count_condition(results_with_target, results_with_target$status == "skipped_dependency", "skipped_dependency_cases"))
	valid_values = arg_values[arg_values$constraint_status == "valid", , drop = FALSE]
	coverage = merge_left(coverage, count_unique(valid_values, c("arg", "value_expr"), "covered_argument_values"))
	coverage = merge_left(coverage, count_unique(pairs, c("arg_pair", "value_pair"), "covered_argument_pairs"))
	int_cols = c("configurable_arg_count", "registry_arg_count", "registry_value_count", "total_cases", "valid_candidate_cases", "pairwise_or_higher_candidate_cases", "targeted_3way_candidate_cases", "executed_ok_cases", "error_cases", "unsupported_cases", "skipped_slow_cases", "skipped_dependency_cases", "covered_argument_values", "covered_argument_pairs")
	for (col in int_cols) coverage[[col]] = fill_int(coverage[[col]])
	coverage$has_more_than_one_configurable_arg = fill_lgl(coverage$has_more_than_one_configurable_arg)
	coverage$has_legal_combination_case = coverage$valid_candidate_cases > 0L
	coverage$has_executed_legal_combination = coverage$executed_ok_cases > 0L
	row.names(coverage) = NULL
	coverage
}

build_uncovered_report = function(coverage) {
	coverage[
		coverage$has_more_than_one_configurable_arg &
			coverage$has_checkmate_contract &
			!coverage$has_legal_combination_case,
		,
		drop = FALSE
	]
}

build_drift_report = function(inventory, contracts, registry, cases) {
	contract_keys = unique(data.frame(target = contracts$target, arg = contracts$arg, stringsAsFactors = FALSE))
	registry_keys = unique(data.frame(target = registry$target, arg = registry$arg, stringsAsFactors = FALSE))
	inventory_targets_set = unique(inventory$target)
	case_targets_set = unique(cases$target)
	rows = list()
	missing_registry = contract_keys[!paste(contract_keys$target, contract_keys$arg, sep = "\r") %in% paste(registry_keys$target, registry_keys$arg, sep = "\r"), , drop = FALSE]
	if (nrow(missing_registry)) {
		rows[[length(rows) + 1L]] = data.frame(drift_type = "contract_missing_registry_arg", missing_registry, detail = "checkmate contract target/arg absent from registry", stringsAsFactors = FALSE)
	}
	extra_registry = registry_keys[!paste(registry_keys$target, registry_keys$arg, sep = "\r") %in% paste(contract_keys$target, contract_keys$arg, sep = "\r"), , drop = FALSE]
	if (nrow(extra_registry)) {
		rows[[length(rows) + 1L]] = data.frame(drift_type = "registry_arg_without_contract", extra_registry, detail = "registry target/arg absent from extracted contracts", stringsAsFactors = FALSE)
	}
	registry_missing_inventory = unique(registry["target"])[!unique(registry$target) %in% inventory_targets_set, , drop = FALSE]
	if (nrow(registry_missing_inventory)) {
		rows[[length(rows) + 1L]] = data.frame(drift_type = "registry_target_missing_inventory", target = registry_missing_inventory$target, arg = "", detail = "registry target absent from public API inventory", stringsAsFactors = FALSE)
	}
	registry_missing_cases = unique(registry["target"])[!unique(registry$target) %in% case_targets_set, , drop = FALSE]
	if (nrow(registry_missing_cases)) {
		rows[[length(rows) + 1L]] = data.frame(drift_type = "registry_target_missing_generated_case", target = registry_missing_cases$target, arg = "", detail = "registry target has no generated case in current case table", stringsAsFactors = FALSE)
	}
	if (!length(rows)) return(data.frame(drift_type = character(), target = character(), arg = character(), detail = character(), stringsAsFactors = FALSE))
	out = do.call(rbind, rows)
	out = out[, c("drift_type", "target", "arg", "detail"), drop = FALSE]
	row.names(out) = NULL
	out
}

build_failure_report = function(cases, results) {
	joined = merge(results, cases[, c("case_id", "target", "argument_coverage_kind", "constraint_status", "constraint_reason", "constraint_sources"), drop = FALSE], by = "case_id", all.x = TRUE)
	failures = joined[joined$status == "error" | grepl("^failed:", joined$invariant_status), , drop = FALSE]
	failures[order(failures$status, failures$target, failures$case_id), , drop = FALSE]
}

build_slow_report = function(cases, results, limit = 50L) {
	joined = merge(results, cases[, c("case_id", "target", "argument_coverage_kind"), drop = FALSE], by = "case_id", all.x = TRUE)
	joined = joined[order(-joined$duration_time_sec, joined$target, joined$case_id), , drop = FALSE]
	head(joined, limit)
}

build_ci_failure_summary = function(failures) {
	if (!nrow(failures)) {
		return(data.frame(total_unexpected_errors = 0L, total_invariant_failures = 0L, ci_should_fail = FALSE, stringsAsFactors = FALSE))
	}
	data.frame(
		total_unexpected_errors = sum(failures$status == "error"),
		total_invariant_failures = sum(grepl("^failed:", failures$invariant_status)),
		ci_should_fail = any(failures$status == "error" | grepl("^failed:", failures$invariant_status)),
		stringsAsFactors = FALSE
	)
}

write_html_report = function(coverage, failures, drift, slow, path = paths$report) {
	html = paste0(
		"<!doctype html><html><head><meta charset=\"utf-8\"><title>Public Argument Combination Report</title>",
		"<style>body{font-family:sans-serif;margin:2rem}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px 8px}</style></head><body>",
		"<h1>Public Argument Combination Report</h1>",
		"<p>Targets: ", nrow(coverage), "</p>",
		"<p>Targets with legal cases: ", sum(coverage$has_legal_combination_case), "</p>",
		"<p>Unexpected failures: ", nrow(failures), "</p>",
		"<p>Registry drift rows: ", nrow(drift), "</p>",
		"<p>Slowest case duration: ", if (nrow(slow)) max(slow$duration_time_sec, na.rm = TRUE) else 0, " sec</p>",
		"</body></html>"
	)
	writeLines(html, path)
	invisible(path)
}

analyze_public_argument_combinations = function() {
	inputs = lapply(required_input_paths, function(name) read_required_csv(paths[[name]], name))
	names(inputs) = required_input_paths
	inventory = inventory_targets(inputs$inventory)
	contracts = contracts_targets(inputs$contracts)
	registry = inputs$registry
	cases = inputs$cases
	results = inputs$results
	coverage = build_coverage_report(inventory, contracts, registry, cases, results)
	uncovered = build_uncovered_report(coverage)
	drift = build_drift_report(inventory, contracts, registry, cases)
	failures = build_failure_report(cases, results)
	slow = build_slow_report(cases, results)
	ci_summary = build_ci_failure_summary(failures)
	write.csv(coverage, paths$coverage, row.names = FALSE)
	write.csv(failures, paths$failures, row.names = FALSE)
	write.csv(uncovered, paths$uncovered, row.names = FALSE)
	write.csv(drift, paths$drift, row.names = FALSE)
	write.csv(slow, paths$slow, row.names = FALSE)
	write.csv(ci_summary, paths$ci_summary, row.names = FALSE)
	write_html_report(coverage, failures, drift, slow)
	invisible(list(
		coverage = coverage,
		failures = failures,
		uncovered = uncovered,
		drift = drift,
		slow = slow,
		ci_summary = ci_summary,
		value_coverage = value_coverage(cases, registry)
	))
}

main = function() {
	out = analyze_public_argument_combinations()
	message("Wrote coverage rows: ", nrow(out$coverage))
	message("Wrote failure rows: ", nrow(out$failures))
	message("Wrote drift rows: ", nrow(out$drift))
	message("Wrote uncovered API rows: ", nrow(out$uncovered))
	if (isTRUE(out$ci_summary$ci_should_fail[1L])) quit(status = 1L)
}

called_as_analyzer_script = function() {
	file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% ""
	identical(basename(sub("^--file=", "", file_arg)), "analyze_public_argument_combinations.R")
}

if (called_as_analyzer_script()) {
	main()
}
