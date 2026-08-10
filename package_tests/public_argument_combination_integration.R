#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/public_argument_combination_integration.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "public_argument_combination_fixtures.R"))

integration_artifact = function(name) file.path(repo_root, "package_tests", name)

integration_paths = list(
	cases = integration_artifact("public_argument_combination_cases.csv"),
	results = integration_artifact("public_argument_combination_results.csv"),
	coverage = integration_artifact("public_argument_combination_coverage.csv"),
	cases_integrated = integration_artifact("public_argument_combination_cases_integrated.csv"),
	results_integrated = integration_artifact("public_argument_combination_results_integrated.csv"),
	coverage_integrated = integration_artifact("public_argument_combination_coverage_integrated.csv")
)

read_optional_csv = function(path) {
	if (!file.exists(path)) return(data.frame(stringsAsFactors = FALSE))
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

clean_character_columns = function(df) {
	if (!nrow(df)) return(df)
	for (col in names(df)) {
		if (is.character(df[[col]])) {
			df[[col]][is.na(df[[col]]) | df[[col]] == "NA"] = ""
		}
	}
	df
}

public_argument_text_columns = c(
	"case_id", "api_kind", "class_name", "function_name", "method_name", "fixture_id",
	"tier", "argument_signature", "argument_labels", "argument_coverage_kind",
	"source_contracts", "constraint_status", "constraint_reason", "constraint_sources",
	"target", "invariants", "response_type", "design_type", "design", "dataset_name",
	"dataset", "dataset_source", "inference_class", "function_run", "status",
	"warning_message", "message_text", "output_type", "output_shape", "invariant_status",
	"error_message"
)

clean_text_columns = function(df, cols = public_argument_text_columns) {
	if (!nrow(df)) return(df)
	for (col in intersect(cols, names(df))) {
		df[[col]] = as.character(df[[col]])
		df[[col]][is.na(df[[col]]) | df[[col]] == "NA"] = ""
	}
	df
}

fixture_join_key_table = function() {
	specs = public_argument_fixture_specs()
	rows = lapply(specs, function(spec) {
		data.frame(
			fixture_id = spec$fixture_id,
			response_type = spec$response_type %||% "",
			design_type = spec$design_type %||% "",
			design = comprehensive_design_label(spec$class_name %||% "", spec$design_type %||% ""),
			dataset_name = spec$fixture_id,
			dataset = spec$fixture_id,
			dataset_source = "public_argument_fixture",
			stringsAsFactors = FALSE
		)
	})
	out = do.call(rbind, rows)
	row.names(out) = NULL
	out
}

comprehensive_design_label = function(class_name, design_type = "") {
	if (nzchar(class_name)) {
		label = sub("^DesignSeqOneByOne", "", class_name)
		label = sub("^Design", "", label)
		if (nzchar(label)) return(label)
	}
	switch(
		design_type,
		fixed = "FixedBernoulli",
		sequential = "Bernoulli",
		stratified = "FixedBlocking",
		clustered = "FixedCluster",
		blocked_cluster = "FixedBlockedCluster",
		matched = "FixedBinaryMatch",
		search = "FixedGreedy",
		design_type
	)
}

infer_function_run = function(df) {
	method = if ("method_name" %in% names(df)) df$method_name else rep("", nrow(df))
	function_name = if ("function_name" %in% names(df)) df$function_name else rep("", nrow(df))
	target = if ("target" %in% names(df)) df$target else rep("", nrow(df))
	out = ifelse(nzchar(method), method, function_name)
	missing = !nzchar(out) & nzchar(target)
	out[missing] = sub("^.*::", "", target[missing])
	out
}

infer_inference_class = function(df) {
	class_name = if ("class_name" %in% names(df)) df$class_name else rep("", nrow(df))
	ifelse(grepl("^Inference", class_name), class_name, "")
}

add_public_argument_join_keys = function(df, fixture_keys = fixture_join_key_table()) {
	if (!nrow(df)) return(df)
	df = clean_text_columns(clean_character_columns(df))
	if (!"fixture_id" %in% names(df)) {
		df$fixture_id = ""
	}
	out = merge(df, fixture_keys, by = "fixture_id", all.x = TRUE, sort = FALSE)
	for (col in c("response_type", "design_type", "design", "dataset_name", "dataset", "dataset_source")) {
		if (!col %in% names(out)) out[[col]] = ""
		out[[col]][is.na(out[[col]])] = ""
	}
	out$inference_class = infer_inference_class(out)
	out$function_run = infer_function_run(out)
	out$argument_combination_runner = TRUE
	out
}

collapse_join_values = function(x) {
	x = unique(as.character(x[!is.na(x) & nzchar(x)]))
	if (!length(x)) return("")
	paste(sort(x), collapse = ";")
}

coverage_join_keys = function(cases_with_keys) {
	if (!nrow(cases_with_keys)) {
		return(data.frame(target = character(), stringsAsFactors = FALSE))
	}
	cols = c("response_type", "design_type", "design", "inference_class", "function_run", "dataset_name", "dataset", "dataset_source")
	rows = lapply(split(cases_with_keys, cases_with_keys$target), function(chunk) {
		data.frame(
			target = chunk$target[1L],
			as.list(vapply(cols, function(col) collapse_join_values(chunk[[col]]), character(1))),
			stringsAsFactors = FALSE,
			check.names = FALSE
		)
	})
	out = do.call(rbind, rows)
	row.names(out) = NULL
	out
}

add_coverage_join_keys = function(coverage, cases_with_keys) {
	if (!nrow(coverage)) return(coverage)
	keys = coverage_join_keys(cases_with_keys)
	if (!nrow(keys)) return(coverage)
	out = merge(coverage, keys, by = "target", all.x = TRUE, sort = FALSE)
	for (col in setdiff(names(keys), "target")) {
		out[[col]][is.na(out[[col]])] = ""
	}
	out
}

integrate_public_argument_combination_outputs = function(
	cases_file = integration_paths$cases,
	results_file = integration_paths$results,
	coverage_file = integration_paths$coverage,
	cases_out = integration_paths$cases_integrated,
	results_out = integration_paths$results_integrated,
	coverage_out = integration_paths$coverage_integrated
) {
	fixture_keys = fixture_join_key_table()
	cases = add_public_argument_join_keys(read_optional_csv(cases_file), fixture_keys)
	results = read_optional_csv(results_file)
	coverage = read_optional_csv(coverage_file)
	if (nrow(results)) {
		results = clean_text_columns(clean_character_columns(results))
		result_meta = cases[, c("case_id", "target", "response_type", "design_type", "design", "inference_class", "function_run", "dataset_name", "dataset", "dataset_source"), drop = FALSE]
		results = merge(results, result_meta, by = "case_id", all.x = TRUE, sort = FALSE)
		for (col in setdiff(names(result_meta), "case_id")) {
			results[[col]][is.na(results[[col]])] = ""
		}
		results$argument_combination_runner = TRUE
	}
	coverage = clean_text_columns(clean_character_columns(add_coverage_join_keys(clean_text_columns(clean_character_columns(coverage)), cases)))
	cases = clean_text_columns(clean_character_columns(cases))
	results = clean_text_columns(clean_character_columns(results))
	write.csv(cases, cases_out, row.names = FALSE)
	write.csv(results, results_out, row.names = FALSE)
	write.csv(coverage, coverage_out, row.names = FALSE)
	invisible(list(cases = cases, results = results, coverage = coverage))
}

main = function() {
	out = integrate_public_argument_combination_outputs()
	message("Wrote integrated public argument case rows: ", nrow(out$cases))
	message("Wrote integrated public argument result rows: ", nrow(out$results))
	message("Wrote integrated public argument coverage rows: ", nrow(out$coverage))
}

called_as_integration_script = function() {
	file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% ""
	identical(basename(sub("^--file=", "", file_arg)), "public_argument_combination_integration.R")
}

if (called_as_integration_script()) {
	main()
}
