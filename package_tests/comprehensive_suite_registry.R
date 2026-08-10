#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/comprehensive_suite_registry.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

artifact = function(name) file.path(repo_root, "package_tests", name)

paths = list(
	inventory = artifact("public_api_inventory.csv"),
	argument_coverage = artifact("public_argument_combination_coverage.csv"),
	baseline_audit = artifact("comprehensive_suite_baseline_audit.csv"),
	exemptions = artifact("comprehensive_suite_exemptions.csv"),
	registry = artifact("comprehensive_suite_registry.csv")
)

read_required_csv = function(path, label) {
	if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

read_optional_csv = function(path) {
	if (!file.exists(path) || isTRUE(file.info(path)$size == 0)) return(data.frame())
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

target_from_inventory = function(inventory) {
	ifelse(nzchar(inventory$method_name), paste(inventory$export_name, inventory$method_name, sep = "::"), inventory$export_name)
}

clean_chr = function(x) {
	x = as.character(x %||% "")
	x[is.na(x)] = ""
	x
}

truthy = function(x) {
	x %in% c(TRUE, "TRUE", "true", "1", 1L)
}

target_text = function(target, class_name, method_name, export_name = "") {
	paste(clean_chr(target), clean_chr(class_name), clean_chr(method_name), clean_chr(export_name), sep = " ")
}

response_family_for = function(target, class_name, method_name, export_name = "") {
	txt = target_text(target, class_name, method_name, export_name)
	if (grepl("InferenceAll|SimpleMeanDiff|SimpleWilcox|MeanDiff", txt, ignore.case = TRUE)) return("all_response")
	if (grepl("Survival|surv|cox|weibull|gehan|logrank|km_|restrictedmean", txt, ignore.case = TRUE)) return("survival")
	if (grepl("Ordinal|ordinal|clmm|cloglog|cauchit|ridit|jonckheere|stereotype|adjacent|continuation", txt, ignore.case = TRUE)) return("ordinal")
	if (grepl("Count|count|Poisson|NegBin|Hurdle|ZeroInflated|ZeroAugmented|QuasiPoisson", txt, ignore.case = TRUE)) return("count")
	if (grepl("Prop|proportion|Beta|Fractional|ZeroOneInflated", txt, ignore.case = TRUE)) return("proportion")
	if (grepl("Incid|incidence|Logit|Probit|Binomial|Risk|Fisher|CMH|Robins|Newcombe|Miettinen", txt, ignore.case = TRUE)) return("incidence")
	if (grepl("Contin|continuous|OLS|Robust|Quantile|Lin|Wilcox|Bai", txt, ignore.case = TRUE)) return("continuous")
	if (grepl("^Design", clean_chr(class_name)) || grepl("^Design", clean_chr(export_name))) return("all_response")
	"unspecified"
}

design_family_for = function(target, class_name, method_name, export_name = "") {
	txt = target_text(target, class_name, method_name, export_name)
	if (grepl("BlockedCluster", txt, ignore.case = TRUE)) return("blocked_cluster")
	if (grepl("Cluster", txt, ignore.case = TRUE)) return("clustered")
	if (grepl("Blocking|Blocks|SPBR|strata|stratified|CMH|Robins", txt, ignore.case = TRUE)) return("stratified")
	if (grepl("BinaryMatch|Matching|matched|KK14|KK21|KKGEE|KKGLMM|CondLogit|OneLik|IVWC|Marginal", txt, ignore.case = TRUE)) return("matched")
	if (grepl("Greedy|Optimal|AOptimal|DOptimal|Rerandomization|PocockSimon|Atkinson|search", txt, ignore.case = TRUE)) return("optimal_search")
	if (grepl("SeqOneByOne|Bernoulli|Efron|Urn|iBCRD|RandomBlockSize|sequential", txt, ignore.case = TRUE)) return("sequential")
	if (grepl("Fixed", txt, ignore.case = TRUE)) return("fixed")
	"unspecified"
}

method_family_for = function(target, class_name, method_name, export_name = "") {
	txt = target_text(target, class_name, method_name, export_name)
	m = clean_chr(method_name)
	if (grepl("SimulationFramework", txt, fixed = TRUE)) return("simulation_framework")
	if (grepl("rand_bootstrap|randomization_bootstrap|brt", txt, ignore.case = TRUE)) return("randomization_bootstrap")
	if (grepl("bayesian_bootstrap|bayes_bootstrap", txt, ignore.case = TRUE)) return("bayesian_bootstrap")
	if (grepl("parametric_bootstrap|param_bootstrap", txt, ignore.case = TRUE)) return("parametric_bootstrap")
	if (grepl("m_out_of_n|subsampling|bootstrap", txt, ignore.case = TRUE)) return("bootstrap")
	if (grepl("bartlett", txt, ignore.case = TRUE)) return("bartlett")
	if (grepl("jackknife", txt, ignore.case = TRUE)) return("jackknife")
	if (grepl("rand|randomization|permutation", txt, ignore.case = TRUE)) return("randomization")
	if (grepl("exact|Fisher|Binomial|Jonckheere|SignTest", txt, ignore.case = TRUE)) return("exact")
	if (grepl("asymp|wald|score|lik_ratio|gradient|confidence_interval|two_sided_pval|pval", txt, ignore.case = TRUE)) return("asymptotic")
	if (grepl("estimate|beta_hat|mean|summary|predict|get_", txt, ignore.case = TRUE)) return("estimate")
	if (identical(m, "initialize") || grepl("^Design|^Inference", clean_chr(class_name))) return("workflow")
	"public_helper"
}

inference_family_for = function(api_kind, target, class_name, method_name, export_name = "") {
	txt = target_text(target, class_name, method_name, export_name)
	if (grepl("^Inference", clean_chr(class_name)) || grepl("^Inference", clean_chr(export_name))) return("inference")
	if (grepl("^Design", clean_chr(class_name)) || grepl("^Design", clean_chr(export_name))) return("design")
	if (grepl("SimulationFramework", txt, fixed = TRUE)) return("simulation")
	if (identical(clean_chr(api_kind), "function")) return("exported_function")
	if (identical(clean_chr(api_kind), "r6_class")) return("r6_class")
	"other"
}

source_runner_for = function(audit_row, exempted) {
	if (isTRUE(exempted)) return("exemption")
	runners = character()
	if (truthy(audit_row$has_argument_combination_coverage)) runners = c(runners, "argument_combinations")
	if (truthy(audit_row$has_comprehensive_workflow_coverage)) runners = c(runners, "comprehensive_tests")
	if (truthy(audit_row$has_focused_testthat_coverage)) runners = c(runners, "testthat")
	if (!length(runners)) return("unassigned")
	paste(unique(runners), collapse = ";")
}

runtime_tier_for = function(row) {
	if (identical(row$coverage_scope, "internal_safety_net")) return("ci")
	if (isTRUE(row$has_exemption)) return("exempted")
	if (truthy(row$has_argument_combination_coverage) || truthy(row$has_focused_testthat_coverage)) return("ci")
	if (truthy(row$has_comprehensive_workflow_coverage)) return("nightly")
	"unassigned"
}

expected_invariants_for = function(row) {
	if (identical(row$coverage_scope, "internal_safety_net")) return("focused internal testthat safety-net assertions pass")
	if (isTRUE(row$has_exemption)) return("documented exemption remains current and unexpired")
	if (truthy(row$has_argument_combination_coverage)) return("public call completes or returns documented supported/non-estimable status")
	if (truthy(row$has_comprehensive_workflow_coverage)) return("end-to-end design/inference workflow records analyzable result status")
	if (truthy(row$has_focused_testthat_coverage)) return("focused public testthat expectations pass")
	"coverage assignment required"
}

required_coverage_for = function(row) {
	if (identical(row$coverage_scope, "internal_safety_net")) return("internal_safety_net")
	if (isTRUE(row$has_exemption)) return("exemption")
	if (identical(row$api_kind, "r6_class") || identical(row$method_name, "initialize")) return("constructor_workflow")
	if (truthy(row$has_more_than_one_configurable_arg)) return("argument_combinations_or_workflow")
	if (grepl("^Inference|^Design|SimulationFramework", paste(row$class_name, row$export_name))) return("public_workflow")
	"focused_public_test_or_workflow"
}

merge_argument_coverage = function(public_rows, argument_coverage) {
	cols = c(
		"target", "has_checkmate_contract", "registry_arg_count", "registry_value_count",
		"total_cases", "valid_candidate_cases", "pairwise_or_higher_candidate_cases",
		"targeted_3way_candidate_cases", "executed_ok_cases", "error_cases",
		"unsupported_cases", "skipped_slow_cases", "skipped_dependency_cases",
		"covered_argument_values", "covered_argument_pairs", "has_legal_combination_case",
		"has_executed_legal_combination"
	)
	extra = argument_coverage[, intersect(cols, names(argument_coverage)), drop = FALSE]
	merge(public_rows, extra, by = "target", all.x = TRUE, sort = FALSE)
}

fill_registry_defaults = function(registry) {
	int_cols = c(
		"registry_arg_count", "registry_value_count", "total_cases", "valid_candidate_cases",
		"pairwise_or_higher_candidate_cases", "targeted_3way_candidate_cases", "executed_ok_cases",
		"error_cases", "unsupported_cases", "skipped_slow_cases", "skipped_dependency_cases",
		"covered_argument_values", "covered_argument_pairs"
	)
	for (col in int_cols) {
		if (!col %in% names(registry)) registry[[col]] = 0L
		registry[[col]][is.na(registry[[col]])] = 0L
		registry[[col]] = as.integer(registry[[col]])
	}
	lgl_cols = c(
		"has_checkmate_contract", "has_legal_combination_case", "has_executed_legal_combination",
		"has_argument_combination_coverage", "has_comprehensive_workflow_coverage",
		"has_focused_testthat_coverage", "has_more_than_one_configurable_arg", "has_exemption"
	)
	for (col in lgl_cols) {
		if (!col %in% names(registry)) registry[[col]] = FALSE
		registry[[col]] = truthy(registry[[col]])
	}
	registry
}

build_public_registry = function(inventory, argument_coverage, audit, exemptions) {
	public_audit = audit[audit$row_type == "public_api", , drop = FALSE]
	inventory$target = target_from_inventory(inventory)
	base = merge(public_audit, inventory[, c("target", "formal_args", "has_assertion_calls", "method_origin_class", "method_inheritance", "source_files"), drop = FALSE], by = "target", all.x = TRUE, sort = FALSE)
	base = merge_argument_coverage(base, argument_coverage)
	exemptions_key = exemptions[, c("target", "exemption_type", "reason", "expiry_date", "owner", "created_date"), drop = FALSE]
	names(exemptions_key) = c("target", "exemption_type", "exemption_reason", "exemption_expiry_date", "exemption_owner", "exemption_created_date")
	base = merge(base, exemptions_key, by = "target", all.x = TRUE, sort = FALSE)
	base$coverage_scope = "public_contract"
	base$has_exemption = nzchar(clean_chr(base$exemption_type))
	base = fill_registry_defaults(base)
	base$response_family = vapply(seq_len(nrow(base)), function(i) response_family_for(base$target[i], base$class_name[i], base$method_name[i], base$export_name[i]), character(1))
	base$design_family = vapply(seq_len(nrow(base)), function(i) design_family_for(base$target[i], base$class_name[i], base$method_name[i], base$export_name[i]), character(1))
	base$inference_family = vapply(seq_len(nrow(base)), function(i) inference_family_for(base$api_kind[i], base$target[i], base$class_name[i], base$method_name[i], base$export_name[i]), character(1))
	base$dataset_fixture_family = ifelse(base$response_family == "unspecified", "generic_public_api", paste(base$response_family, base$design_family, sep = "_"))
	base$method_family = vapply(seq_len(nrow(base)), function(i) method_family_for(base$target[i], base$class_name[i], base$method_name[i], base$export_name[i]), character(1))
	base$source_runner = vapply(seq_len(nrow(base)), function(i) source_runner_for(base[i, , drop = FALSE], base$has_exemption[i]), character(1))
	base$runtime_tier = vapply(seq_len(nrow(base)), function(i) runtime_tier_for(base[i, , drop = FALSE]), character(1))
	base$expected_invariants = vapply(seq_len(nrow(base)), function(i) expected_invariants_for(base[i, , drop = FALSE]), character(1))
	base$required_coverage = vapply(seq_len(nrow(base)), function(i) required_coverage_for(base[i, , drop = FALSE]), character(1))
	base
}

build_internal_registry = function(audit) {
	internal = audit[audit$row_type == "testthat_internal", , drop = FALSE]
	if (!nrow(internal)) return(data.frame())
	internal$export_name = ""
	internal$api_kind = "testthat_internal"
	internal$class_name = ""
	internal$method_name = ""
	internal$configurable_arg_count = 0L
	internal$has_more_than_one_configurable_arg = FALSE
	internal$formal_args = ""
	internal$has_assertion_calls = FALSE
	internal$method_origin_class = ""
	internal$method_inheritance = ""
	internal$source_files = internal$testthat_file
	internal$coverage_scope = "internal_safety_net"
	internal$has_argument_combination_coverage = FALSE
	internal$has_comprehensive_workflow_coverage = FALSE
	internal$has_focused_testthat_coverage = FALSE
	internal$has_checkmate_contract = FALSE
	internal$registry_arg_count = 0L
	internal$registry_value_count = 0L
	internal$total_cases = 0L
	internal$valid_candidate_cases = 0L
	internal$pairwise_or_higher_candidate_cases = 0L
	internal$targeted_3way_candidate_cases = 0L
	internal$executed_ok_cases = 0L
	internal$error_cases = 0L
	internal$unsupported_cases = 0L
	internal$skipped_slow_cases = 0L
	internal$skipped_dependency_cases = 0L
	internal$covered_argument_values = 0L
	internal$covered_argument_pairs = 0L
	internal$has_legal_combination_case = FALSE
	internal$has_executed_legal_combination = FALSE
	internal$exemption_type = ""
	internal$exemption_reason = ""
	internal$exemption_expiry_date = ""
	internal$exemption_owner = ""
	internal$exemption_created_date = ""
	internal$has_exemption = FALSE
	internal$response_family = "internal"
	internal$design_family = "internal"
	internal$inference_family = "internal"
	internal$dataset_fixture_family = "testthat_internal"
	internal$method_family = "internal_safety_net"
	internal$source_runner = "internal_safety_net"
	internal$runtime_tier = "ci"
	internal$expected_invariants = "focused internal testthat safety-net assertions pass"
	internal$required_coverage = "internal_safety_net"
	internal
}

registry_columns = c(
	"target", "coverage_scope", "api_kind", "export_name", "class_name", "method_name",
	"method_origin_class", "method_inheritance", "configurable_arg_count", "formal_args",
	"has_more_than_one_configurable_arg", "has_assertion_calls", "response_family",
	"design_family", "inference_family", "dataset_fixture_family", "method_family",
	"runtime_tier", "required_coverage", "expected_invariants", "source_runner",
	"coverage_status", "responsible_runner", "has_argument_combination_coverage",
	"has_comprehensive_workflow_coverage", "has_focused_testthat_coverage",
	"has_checkmate_contract", "registry_arg_count", "registry_value_count", "total_cases",
	"valid_candidate_cases", "pairwise_or_higher_candidate_cases", "targeted_3way_candidate_cases",
	"executed_ok_cases", "error_cases", "unsupported_cases", "skipped_slow_cases",
	"skipped_dependency_cases", "covered_argument_values", "covered_argument_pairs",
	"has_legal_combination_case", "has_executed_legal_combination", "has_exemption",
	"exemption_type", "exemption_reason", "exemption_expiry_date", "exemption_owner",
	"exemption_created_date", "source_files", "testthat_file", "internal_symbols",
	"internal_test_classification", "classification_reason"
)

build_comprehensive_suite_registry = function() {
	inventory = read_required_csv(paths$inventory, "public_api_inventory")
	argument_coverage = read_required_csv(paths$argument_coverage, "public_argument_combination_coverage")
	audit = read_required_csv(paths$baseline_audit, "comprehensive_suite_baseline_audit")
	exemptions = read_optional_csv(paths$exemptions)
	if (!nrow(exemptions)) {
		exemptions = data.frame(target = character(), exemption_type = character(), reason = character(), expiry_date = character(), owner = character(), created_date = character(), stringsAsFactors = FALSE)
	}
	public = build_public_registry(inventory, argument_coverage, audit, exemptions)
	internal = build_internal_registry(audit)
	registry = if (nrow(internal)) rbind(public[, names(internal), drop = FALSE], internal) else public
	for (col in setdiff(registry_columns, names(registry))) registry[[col]] = ""
	registry = registry[, registry_columns, drop = FALSE]
	registry = registry[order(registry$coverage_scope, registry$target), , drop = FALSE]
	row.names(registry) = NULL
	registry
}

validate_comprehensive_suite_registry = function(registry) {
	public = registry[registry$coverage_scope == "public_contract", , drop = FALSE]
	if (!nrow(public)) stop("Registry has no public_contract rows.", call. = FALSE)
	if (any(!nzchar(public$target))) stop("Public registry rows have blank targets.", call. = FALSE)
	if (any(!nzchar(public$required_coverage))) stop("Public registry rows have blank required_coverage.", call. = FALSE)
	if (any(!nzchar(public$source_runner))) stop("Public registry rows have blank source_runner.", call. = FALSE)
	uncovered = public[public$coverage_status == "uncovered", , drop = FALSE]
	if (nrow(uncovered) && any(!truthy(uncovered$has_exemption))) {
		stop("Uncovered public registry rows are missing exemption hooks.", call. = FALSE)
	}
	internal = registry[registry$coverage_scope == "internal_safety_net", , drop = FALSE]
	if (nrow(internal) && any(internal$source_runner != "internal_safety_net")) {
		stop("Internal registry rows must use source_runner=internal_safety_net.", call. = FALSE)
	}
	invisible(TRUE)
}

write_comprehensive_suite_registry = function(registry = build_comprehensive_suite_registry(), path = paths$registry) {
	validate_comprehensive_suite_registry(registry)
	write.csv(registry, path, row.names = FALSE)
	invisible(registry)
}

main = function() {
	registry = write_comprehensive_suite_registry()
	message("Wrote comprehensive suite registry rows: ", nrow(registry))
	message("Coverage scopes:")
	print(table(registry$coverage_scope, useNA = "ifany"))
	message("Source runners:")
	print(table(registry$source_runner, useNA = "ifany"))
	message("Required coverage:")
	print(table(registry$required_coverage, useNA = "ifany"))
}

main()
