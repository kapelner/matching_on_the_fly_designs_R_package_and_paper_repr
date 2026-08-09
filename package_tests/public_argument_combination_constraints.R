#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

valid_constraint_statuses = c(
	"valid",
	"unsupported",
	"skipped_slow",
	"skipped_dependency",
	"invalid_registry"
)

tier_rank = c(smoke = 1L, ci = 2L, nightly = 3L, release = 4L)

constraint_result = function(is_valid, status, reason, source) {
	if (!status %in% valid_constraint_statuses) {
		stop("Unknown constraint status: ", status, call. = FALSE)
	}
	data.frame(
		is_valid = isTRUE(is_valid),
		status = status,
		reason = as.character(reason),
		source = as.character(source),
		stringsAsFactors = FALSE
	)
}

valid_result = function(source, reason = "constraint not applicable or satisfied") {
	constraint_result(TRUE, "valid", reason, source)
}

arg_value = function(args, name, default = NULL) {
	if (is.null(args) || !name %in% names(args)) return(default)
	args[[name]]
}

fixture_value = function(fixture, name, default = NULL) {
	if (is.null(fixture) || !name %in% names(fixture)) return(default)
	fixture[[name]]
}

as_numberish = function(x) {
	if (is.null(x) || length(x) == 0L) return(NA_real_)
	if (is.numeric(x) || is.integer(x)) return(as.numeric(x[1L]))
	if (is.character(x)) return(suppressWarnings(as.numeric(sub("L$", "", x[1L]))))
	NA_real_
}

as_logicalish = function(x) {
	if (is.null(x) || length(x) == 0L) return(NA)
	if (is.logical(x)) return(isTRUE(x[1L]))
	if (is.character(x)) {
		val = toupper(trimws(x[1L]))
		if (val == "TRUE") return(TRUE)
		if (val == "FALSE") return(FALSE)
	}
	NA
}

as_character_values = function(x) {
	if (is.null(x) || length(x) == 0L) return(character())
	if (inherits(x, "formula")) return(deparse(x, width.cutoff = 500L))
	if (is.factor(x)) return(as.character(x))
	as.character(x)
}

normalize_formula_arg = function(x) {
	if (is.null(x) || length(x) == 0L) return(NULL)
	if (inherits(x, "formula")) return(x)
	if (!is.character(x)) return(NULL)
	val = trimws(x[1L])
	if (!nzchar(val) || val == "NULL" || grepl("^<fixture:", val)) return(NULL)
	if (!grepl("^~", val)) return(NULL)
	tryCatch(stats::as.formula(val), error = function(e) NULL)
}

formula_arg_names = function(args) {
	nms = names(args %||% list())
	nms[grepl("formula", nms, fixed = TRUE)]
}

constraint_min_usable_samples_le_B = function(args, fixture, tier, registry_entry) {
	min_usable = as_numberish(arg_value(args, "min_number_usable_samples"))
	B = as_numberish(arg_value(args, "B"))
	if (is.na(min_usable) || is.na(B)) {
		return(valid_result("constraint_min_usable_samples_le_B"))
	}
	if (min_usable <= B) {
		return(valid_result("constraint_min_usable_samples_le_B", "min_number_usable_samples is at most B"))
	}
	constraint_result(
		FALSE,
		"invalid_registry",
		sprintf("min_number_usable_samples (%s) exceeds B (%s)", min_usable, B),
		"constraint_min_usable_samples_le_B"
	)
}

constraint_bca_requires_jackknife = function(args, fixture, tier, registry_entry) {
	type = as_character_values(arg_value(args, "type"))
	if (!"bca" %in% type) return(valid_result("constraint_bca_requires_jackknife"))
	if (!isTRUE(fixture_value(fixture, "supports_jackknife", FALSE))) {
		return(constraint_result(FALSE, "unsupported", "type = 'bca' requires jackknife support", "constraint_bca_requires_jackknife"))
	}
	n_units = as_numberish(fixture_value(fixture, "n_exchangeable_units", fixture_value(fixture, "n")))
	if (!is.na(n_units) && n_units < 3L) {
		return(constraint_result(FALSE, "unsupported", "type = 'bca' requires at least three exchangeable units", "constraint_bca_requires_jackknife"))
	}
	valid_result("constraint_bca_requires_jackknife", "bca jackknife requirements are satisfied")
}

constraint_studentized_requires_finite_se = function(args, fixture, tier, registry_entry) {
	type = as_character_values(arg_value(args, "type"))
	if (!"studentized" %in% type) return(valid_result("constraint_studentized_requires_finite_se"))
	if (isTRUE(fixture_value(fixture, "has_finite_standard_errors", FALSE))) {
		return(valid_result("constraint_studentized_requires_finite_se", "studentized standard errors are finite"))
	}
	if (isTRUE(fixture_value(fixture, "documented_nonestimable_behavior", FALSE))) {
		return(valid_result("constraint_studentized_requires_finite_se", "non-estimable studentized behavior is documented"))
	}
	constraint_result(FALSE, "unsupported", "type = 'studentized' requires finite standard errors or documented non-estimable behavior", "constraint_studentized_requires_finite_se")
}

constraint_use_rcpp_false_optional_backend = function(args, fixture, tier, registry_entry) {
	use_rcpp = as_logicalish(arg_value(args, "use_rcpp"))
	if (is.na(use_rcpp) || isTRUE(use_rcpp)) return(valid_result("constraint_use_rcpp_false_optional_backend"))
	if (isTRUE(fixture_value(fixture, "optional_backend_available", FALSE))) {
		return(valid_result("constraint_use_rcpp_false_optional_backend", "optional non-Rcpp backend is available"))
	}
	constraint_result(FALSE, "skipped_dependency", "use_rcpp = FALSE requires the documented optional backend", "constraint_use_rcpp_false_optional_backend")
}

constraint_show_progress_true_tier = function(args, fixture, tier, registry_entry) {
	show_progress = as_logicalish(arg_value(args, "show_progress"))
	if (is.na(show_progress) || !isTRUE(show_progress)) return(valid_result("constraint_show_progress_true_tier"))
	if (isTRUE(fixture_value(fixture, "capture_progress_output", FALSE))) {
		return(valid_result("constraint_show_progress_true_tier", "progress output is intentionally captured"))
	}
	if (tier %in% c("nightly", "release")) {
		return(valid_result("constraint_show_progress_true_tier", "progress output is outside automated CI tier"))
	}
	constraint_result(FALSE, "skipped_slow", "show_progress = TRUE is not used in smoke/ci unless output is captured", "constraint_show_progress_true_tier")
}

constraint_transform_responses = function(args, fixture, tier, registry_entry) {
	transform_responses = as_character_values(arg_value(args, "transform_responses"))
	if (!length(transform_responses)) return(valid_result("constraint_transform_responses"))
	legal = as_character_values(fixture_value(fixture, "legal_response_transforms", character()))
	if (!length(legal)) {
		return(constraint_result(FALSE, "unsupported", "fixture does not declare legal response transforms", "constraint_transform_responses"))
	}
	missing = setdiff(transform_responses, legal)
	if (!length(missing)) {
		return(valid_result("constraint_transform_responses", "response transforms are compatible with fixture"))
	}
	constraint_result(FALSE, "unsupported", paste("transform_responses not supported by fixture:", paste(missing, collapse = ", ")), "constraint_transform_responses")
}

constraint_formula_columns = function(args, fixture, tier, registry_entry) {
	formula_args = formula_arg_names(args)
	if (!length(formula_args)) return(valid_result("constraint_formula_columns"))
	available_columns = as_character_values(fixture_value(fixture, "available_columns", character()))
	if (!length(available_columns)) {
		return(constraint_result(FALSE, "invalid_registry", "fixture must declare available_columns for formula arguments", "constraint_formula_columns"))
	}
	bad = character()
	for (nm in formula_args) {
		f = normalize_formula_arg(arg_value(args, nm))
		if (is.null(f)) next
		vars = all.vars(f)
		bad = c(bad, setdiff(vars, available_columns))
	}
	bad = unique(bad)
	if (!length(bad)) return(valid_result("constraint_formula_columns", "formula variables are present in fixture"))
	constraint_result(FALSE, "unsupported", paste("formula variables are missing from fixture:", paste(bad, collapse = ", ")), "constraint_formula_columns")
}

constraint_strata_cols_valid = function(args, fixture, tier, registry_entry) {
	strata_cols = as_character_values(arg_value(args, "strata_cols"))
	if (!length(strata_cols)) return(valid_result("constraint_strata_cols_valid"))
	available_columns = as_character_values(fixture_value(fixture, "available_columns", character()))
	missing = setdiff(strata_cols, available_columns)
	if (length(missing)) {
		return(constraint_result(FALSE, "unsupported", paste("strata_cols are missing from fixture:", paste(missing, collapse = ", ")), "constraint_strata_cols_valid"))
	}
	if (isTRUE(fixture_value(fixture, "strata_cols_must_be_factor", FALSE))) {
		factor_columns = as_character_values(fixture_value(fixture, "factor_columns", character()))
		not_factor = setdiff(strata_cols, factor_columns)
		if (length(not_factor)) {
			return(constraint_result(FALSE, "unsupported", paste("strata_cols must be factor/discretized:", paste(not_factor, collapse = ", ")), "constraint_strata_cols_valid"))
		}
	}
	valid_result("constraint_strata_cols_valid", "strata_cols are compatible with fixture")
}

constraint_cluster_col_valid = function(args, fixture, tier, registry_entry) {
	cluster_col = as_character_values(arg_value(args, "cluster_col"))
	if (!length(cluster_col) || !nzchar(cluster_col[1L])) return(valid_result("constraint_cluster_col_valid"))
	available_columns = as_character_values(fixture_value(fixture, "available_columns", character()))
	if (!cluster_col[1L] %in% available_columns) {
		return(constraint_result(FALSE, "unsupported", paste("cluster_col is missing from fixture:", cluster_col[1L]), "constraint_cluster_col_valid"))
	}
	cluster_columns = as_character_values(fixture_value(fixture, "cluster_columns", character()))
	if (length(cluster_columns) && !cluster_col[1L] %in% cluster_columns) {
		return(constraint_result(FALSE, "unsupported", paste("cluster_col is not declared as a valid cluster ID column:", cluster_col[1L]), "constraint_cluster_col_valid"))
	}
	valid_result("constraint_cluster_col_valid", "cluster_col is compatible with fixture")
}

constraint_strata_cluster_no_conflict = function(args, fixture, tier, registry_entry) {
	strata_cols = as_character_values(arg_value(args, "strata_cols"))
	cluster_col = as_character_values(arg_value(args, "cluster_col"))
	if (!length(strata_cols) || !length(cluster_col) || !nzchar(cluster_col[1L])) {
		return(valid_result("constraint_strata_cluster_no_conflict"))
	}
	if (!cluster_col[1L] %in% strata_cols) {
		return(valid_result("constraint_strata_cluster_no_conflict", "strata_cols and cluster_col do not overlap"))
	}
	if (isTRUE(fixture_value(fixture, "allows_strata_cluster_overlap", FALSE))) {
		return(valid_result("constraint_strata_cluster_no_conflict", "fixture permits strata/cluster overlap"))
	}
	constraint_result(FALSE, "unsupported", "strata_cols and cluster_col overlap but the fixture does not support that combination", "constraint_strata_cluster_no_conflict")
}

constraint_fixed_blocking_n_compatible = function(args, fixture, tier, registry_entry) {
	n = as_numberish(arg_value(args, "n", fixture_value(fixture, "n")))
	if (is.na(n)) return(valid_result("constraint_fixed_blocking_n_compatible"))
	block_size = as_numberish(arg_value(args, "block_size"))
	if (is.na(block_size)) block_size = as_numberish(arg_value(args, "block_sizes"))
	if (!is.na(block_size) && block_size > 0L && n %% block_size != 0L) {
		return(constraint_result(FALSE, "unsupported", sprintf("n (%s) is not divisible by block size (%s)", n, block_size), "constraint_fixed_blocking_n_compatible"))
	}
	strata_levels = as_numberish(fixture_value(fixture, "strata_level_count"))
	if (!is.na(strata_levels) && strata_levels > n) {
		return(constraint_result(FALSE, "unsupported", sprintf("strata levels (%s) exceed n (%s)", strata_levels, n), "constraint_fixed_blocking_n_compatible"))
	}
	valid_result("constraint_fixed_blocking_n_compatible", "fixed/blocking dimensions are compatible")
}

constraint_survival_censoring_supported = function(args, fixture, tier, registry_entry) {
	if (!identical(fixture_value(fixture, "response_type"), "survival")) {
		return(valid_result("constraint_survival_censoring_supported"))
	}
	if (!isTRUE(fixture_value(fixture, "has_censoring", FALSE))) {
		return(valid_result("constraint_survival_censoring_supported", "survival fixture is uncensored"))
	}
	if (isTRUE(fixture_value(fixture, "supports_censoring", FALSE))) {
		return(valid_result("constraint_survival_censoring_supported", "survival method supports censoring"))
	}
	constraint_result(FALSE, "unsupported", "survival fixture has censoring but target does not publicly support censoring", "constraint_survival_censoring_supported")
}

constraint_optional_dependency_paths = function(args, fixture, tier, registry_entry) {
	deps = registry_entry$optional_dependencies %||% character()
	deps = deps[!is.na(deps) & nzchar(as.character(deps)) & deps != "none"]
	if (!length(deps)) return(valid_result("constraint_optional_dependency_paths"))
	available = as_character_values(fixture_value(fixture, "available_optional_dependencies", character()))
	missing = setdiff(as.character(deps), available)
	if (!length(missing)) return(valid_result("constraint_optional_dependency_paths", "optional dependencies are available"))
	constraint_result(FALSE, "skipped_dependency", paste("optional dependencies are unavailable:", paste(missing, collapse = ", ")), "constraint_optional_dependency_paths")
}

constraint_known_slow_paths = function(args, fixture, tier, registry_entry) {
	labels = registry_entry$slow_path_labels %||% character()
	labels = labels[!is.na(labels) & nzchar(as.character(labels)) & labels != "none"]
	B = as_numberish(arg_value(args, "B"))
	r = as_numberish(arg_value(args, "r"))
	large_runtime_count = any(c(B, r) > c(smoke = 5L, ci = 25L)[tier] %||% Inf, na.rm = TRUE)
	if (!length(labels) && !large_runtime_count) return(valid_result("constraint_known_slow_paths"))
	min_tier = fixture_value(fixture, "min_slow_tier", "nightly")
	if (is.na(tier_rank[tier]) || is.na(tier_rank[min_tier]) || tier_rank[tier] >= tier_rank[min_tier]) {
		return(valid_result("constraint_known_slow_paths", "slow path is allowed at requested tier"))
	}
	reason = if (length(labels)) {
		paste("known slow path below allowed tier:", paste(labels, collapse = ", "))
	} else {
		"runtime count exceeds tier threshold"
	}
	constraint_result(FALSE, "skipped_slow", reason, "constraint_known_slow_paths")
}

public_argument_combination_constraints = list(
	constraint_min_usable_samples_le_B,
	constraint_bca_requires_jackknife,
	constraint_studentized_requires_finite_se,
	constraint_use_rcpp_false_optional_backend,
	constraint_show_progress_true_tier,
	constraint_transform_responses,
	constraint_formula_columns,
	constraint_strata_cols_valid,
	constraint_cluster_col_valid,
	constraint_strata_cluster_no_conflict,
	constraint_fixed_blocking_n_compatible,
	constraint_survival_censoring_supported,
	constraint_optional_dependency_paths,
	constraint_known_slow_paths
)

evaluate_public_argument_constraints = function(args = list(), fixture = list(), tier = "ci", registry_entry = list()) {
	rows = lapply(public_argument_combination_constraints, function(fn) {
		out = fn(args = args, fixture = fixture, tier = tier, registry_entry = registry_entry)
		out$constraint_name = out$source
		out
	})
	out = do.call(rbind, rows)
	out = out[, c("constraint_name", "is_valid", "status", "reason", "source"), drop = FALSE]
	row.names(out) = NULL
	out
}

overall_constraint_status = function(results) {
	if (!is.data.frame(results) || !nrow(results)) {
		return(list(is_valid = FALSE, status = "invalid_registry", reason = "constraint results are empty"))
	}
	for (status in c("invalid_registry", "skipped_dependency", "unsupported", "skipped_slow")) {
		hits = results[results$status == status, , drop = FALSE]
		if (nrow(hits)) {
			return(list(
				is_valid = FALSE,
				status = status,
				reason = paste(unique(hits$reason), collapse = " | ")
			))
		}
	}
	list(is_valid = TRUE, status = "valid", reason = "all constraints satisfied")
}

constraint_context_result = function(args = list(), fixture = list(), tier = "ci", registry_entry = list()) {
	results = evaluate_public_argument_constraints(args, fixture, tier, registry_entry)
	overall = overall_constraint_status(results)
	list(
		is_valid = overall$is_valid,
		status = overall$status,
		reason = overall$reason,
		results = results
	)
}

default_constraint_fixture = function() {
	list(
		fixture_id = "default_constraint_fixture",
		response_type = "continuous",
		n = 20L,
		n_exchangeable_units = 20L,
		supports_jackknife = TRUE,
		has_finite_standard_errors = TRUE,
		documented_nonestimable_behavior = FALSE,
		optional_backend_available = TRUE,
		capture_progress_output = FALSE,
		legal_response_transforms = c("identity", "rank", "log"),
		available_columns = c("y", "x1", "x2", "stratum", "cluster_id"),
		factor_columns = c("stratum"),
		cluster_columns = c("cluster_id"),
		strata_cols_must_be_factor = TRUE,
		allows_strata_cluster_overlap = FALSE,
		strata_level_count = 2L,
		has_censoring = FALSE,
		supports_censoring = FALSE,
		available_optional_dependencies = character(),
		min_slow_tier = "nightly"
	)
}

example_constraint_cases = function() {
	fixture = default_constraint_fixture()
	list(
		valid = list(args = list(B = 25L, min_number_usable_samples = 5L, type = "basic", show_progress = FALSE), fixture = fixture, tier = "ci", registry_entry = list()),
		invalid_registry = list(args = list(B = 5L, min_number_usable_samples = 10L), fixture = fixture, tier = "ci", registry_entry = list()),
		unsupported = list(args = list(type = "bca"), fixture = within(fixture, supports_jackknife <- FALSE), tier = "ci", registry_entry = list()),
		skipped_dependency = list(args = list(use_rcpp = FALSE), fixture = within(fixture, optional_backend_available <- FALSE), tier = "ci", registry_entry = list()),
		skipped_slow = list(args = list(show_progress = TRUE), fixture = fixture, tier = "ci", registry_entry = list())
	)
}
