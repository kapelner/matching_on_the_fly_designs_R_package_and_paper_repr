#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/public_argument_contract_registry.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

nz = function(x) {
	!is.null(x) && length(x) > 0L && !is.na(x[1L]) && nzchar(as.character(x[1L]))
}

split_semicolon = function(x) {
	if (!nz(x)) return(character())
	vals = strsplit(as.character(x), ";", fixed = TRUE)[[1L]]
	vals[nzchar(vals)]
}

unique_values = function(df) {
	df = unique(df)
	df[order(df$tier, df$label, df$value_expr), , drop = FALSE]
}

value_row = function(label, value_expr, tier = "ci", source = "registry", rationale = "") {
	data.frame(
		label = as.character(label),
		value_expr = as.character(value_expr),
		tier = as.character(tier),
		source = as.character(source),
		rationale = as.character(rationale),
		stringsAsFactors = FALSE
	)
}

value_rows = function(...) {
	rows = list(...)
	if (length(rows) == 1L && is.list(rows[[1L]]) && !is.data.frame(rows[[1L]])) {
		rows = rows[[1L]]
	}
	rows = rows[!vapply(rows, is.null, logical(1))]
	if (!length(rows)) {
		return(value_row("fixture", "<fixture>", "smoke", "fixture", "fixture-provided argument"))
	}
	unique_values(do.call(rbind, rows))
}

tier_for_choice = function(i, n) {
	if (n <= 4L) return("ci")
	if (i <= 2L) return("ci")
	if (i <= 6L) return("nightly")
	"release"
}

default_value_row = function(default_expr) {
	if (!nz(default_expr) || identical(default_expr, "<missing>")) return(NULL)
	value_row("default", default_expr, "smoke", "formal_default", "formal default expression")
}

numeric_domain = function(contract, arg) {
	rows = list(default_value_row(contract$default_expr))
	lower = suppressWarnings(as.numeric(contract$lower))
	upper = suppressWarnings(as.numeric(contract$upper))
	has_lower = length(lower) == 1L && is.finite(lower)
	has_upper = length(upper) == 1L && is.finite(upper)
	
	if (identical(arg, "alpha")) {
		rows = c(rows, list(
			value_row("narrow", "0.01", "ci", "hand_rule", "representative legal alpha"),
			value_row("standard", "0.05", "smoke", "hand_rule", "standard alpha"),
			value_row("wide", "0.10", "ci", "hand_rule", "representative legal alpha")
		))
	} else if (identical(arg, "delta")) {
		rows = c(rows, list(
			value_row("null", "0", "smoke", "hand_rule", "default null effect"),
			value_row("negative_small", "-0.5", "nightly", "hand_rule", "representative nonzero effect"),
			value_row("positive_small", "0.5", "ci", "hand_rule", "representative nonzero effect")
		))
	} else if (identical(arg, "pval_epsilon")) {
		rows = c(rows, list(
			value_row("standard", "0.01", "smoke", "hand_rule", "root-search p-value tolerance"),
			value_row("small", "0.001", "nightly", "hand_rule", "stricter root-search p-value tolerance")
		))
	} else if (has_lower || has_upper) {
		if (has_lower) {
			near_lower = if (lower <= 0) max(.Machine$double.eps, lower + .Machine$double.eps) else lower
			rows = c(rows, list(value_row("near_lower", format(near_lower, scientific = FALSE), "nightly", "checkmate_bound", "legal value near lower bound")))
		}
		if (has_upper) {
			near_upper = if (upper < 1) upper - .Machine$double.eps else upper
			rows = c(rows, list(value_row("near_upper", format(near_upper, scientific = FALSE), "nightly", "checkmate_bound", "legal value near upper bound")))
		}
		if (has_lower && has_upper && lower < upper) {
			rows = c(rows, list(value_row("interior", format((lower + upper) / 2, scientific = FALSE), "ci", "checkmate_bound", "interior legal numeric value")))
		}
	} else {
		rows = c(rows, list(value_row("zero", "0", "ci", "hand_rule", "generic finite numeric value")))
	}
	value_rows(rows)
}

count_domain = function(contract, arg) {
	rows = list(default_value_row(contract$default_expr))
	if (arg %in% c("B", "r", "nsim", "num_boot")) {
		rows = c(rows, list(
			value_row("tiny", "5L", "smoke", "runtime_tier", "tiny resampling count"),
			value_row("ci", "25L", "ci", "runtime_tier", "CI resampling count"),
			value_row("nightly", "101L", "nightly", "runtime_tier", "nightly resampling count")
		))
	} else if (grepl("iteration|maxit|attempt", arg, ignore.case = TRUE)) {
		rows = c(rows, list(
			value_row("one", "1L", "smoke", "checkmate_count", "minimum positive count"),
			value_row("small", "5L", "ci", "runtime_tier", "small iteration/attempt count")
		))
	} else {
		rows = c(rows, list(
			value_row("one", "1L", "smoke", "checkmate_count", "minimum positive count"),
			value_row("small", "5L", "ci", "checkmate_count", "small positive count")
		))
	}
	value_rows(rows)
}

choice_domain = function(contract) {
	choices = split_semicolon(contract$choices)
	rows = list(default_value_row(contract$default_expr))
	if (length(choices)) {
		rows = c(rows, lapply(seq_along(choices), function(i) {
			value_row(
				paste0("choice_", choices[i]),
				shQuote(choices[i]),
				tier_for_choice(i, length(choices)),
				"checkmate_choices",
				"literal checkmate choice"
			)
		}))
	}
	value_rows(rows)
}

flag_domain = function(contract, arg) {
	rows = list(default_value_row(contract$default_expr))
	true_tier = if (identical(arg, "show_progress")) "nightly" else "ci"
	rows = c(rows, list(
		value_row("false", "FALSE", "smoke", "checkmate_flag", "legal FALSE flag"),
		value_row("true", "TRUE", true_tier, "checkmate_flag", "legal TRUE flag")
	))
	value_rows(rows)
}

formula_domain = function(contract) {
	rows = list(default_value_row(contract$default_expr))
	null_tier = if (identical(tolower(contract$null_ok), "true")) "ci" else "release"
	rows = c(rows, list(
		value_row("null", "NULL", null_tier, "formula_rule", "legal NULL formula when allowed"),
		value_row("intercept_only", "~ 1", "smoke", "formula_rule", "intercept-only formula"),
		value_row("all_covariates", "~ .", "ci", "formula_rule", "all-covariate formula"),
		value_row("reduced_fixture", "~ x1", "nightly", "formula_rule", "fixture-specific reduced formula")
	))
	value_rows(rows)
}

fixture_domain = function(contract, kind) {
	label = paste0("fixture_", kind)
	value_row(label, paste0("<fixture:", kind, ">"), "smoke", "fixture", paste(kind, "provided by fixture layer"))
}

domain_for_contract = function(contract) {
	arg = as.character(contract$arg)
	assertion = as.character(contract$assertion)
	switch(
		assertion,
		assertChoice = choice_domain(contract),
		assertSubset = choice_domain(contract),
		assertResponseType = choice_domain(contract),
		assertClass = choice_domain(contract),
		assertNumeric = numeric_domain(contract, arg),
		assertIntegerish = numeric_domain(contract, arg),
		assertCount = count_domain(contract, arg),
		assertFlag = flag_domain(contract, arg),
		assertLogical = flag_domain(contract, arg),
		assertFormula = formula_domain(contract),
		assertList = fixture_domain(contract, "list"),
		assertDataFrame = fixture_domain(contract, "data.frame"),
		assertMatrix = fixture_domain(contract, "matrix"),
		assertString = value_rows(default_value_row(contract$default_expr), value_row("fixture_string", "<fixture:string>", "smoke", "fixture", "fixture-provided string")),
		assertNoCensoring = value_row("no_censoring", "<fixture:no_censoring>", "smoke", "unsupported_context", "public method requires uncensored data"),
		fixture_domain(contract, "value")
	)
}

combine_arg_domains = function(contracts_for_arg) {
	rows = lapply(seq_len(nrow(contracts_for_arg)), function(i) {
		dom = domain_for_contract(contracts_for_arg[i, ])
		dom$contract_assertion = contracts_for_arg$assertion[i]
		dom$contract_source_line = contracts_for_arg$source_line[i]
		dom
	})
	unique_values(do.call(rbind, rows))
}

target_key = function(api_name, method_name) {
	method_name = as.character(method_name %||% "")
	api_name = as.character(api_name %||% "")
	ifelse(nzchar(method_name), paste(api_name, method_name, sep = "::"), api_name)
}

constraints_for_entry = function(arg_names, value_table) {
	constraints = character()
	if (all(c("min_number_usable_samples", "B") %in% arg_names)) {
		constraints = c(constraints, "min_number_usable_samples <= B")
	}
	if ("type" %in% arg_names && any(grepl("bca", value_table$value_expr, fixed = TRUE))) {
		constraints = c(constraints, "type != 'bca' || fixture$supports_jackknife")
	}
	if ("type" %in% arg_names && any(grepl("studentized", value_table$value_expr, fixed = TRUE))) {
		constraints = c(constraints, "type != 'studentized' || fixture$has_finite_standard_errors")
	}
	if ("use_rcpp" %in% arg_names) {
		constraints = c(constraints, "use_rcpp != FALSE || fixture$optional_backend_available")
	}
	if ("show_progress" %in% arg_names) {
		constraints = c(constraints, "show_progress != TRUE || tier %in% c('nightly', 'release')")
	}
	if ("transform_responses" %in% arg_names) {
		constraints = c(constraints, "transform_responses %in% fixture$legal_response_transforms")
	}
	if (any(grepl("formula", arg_names, fixed = TRUE))) {
		constraints = c(constraints, "formula variables must be present in fixture$available_columns")
	}
	if ("assertNoCensoring" %in% value_table$contract_assertion) {
		constraints = c(constraints, "fixture$has_censoring == FALSE")
	}
	if (!length(constraints)) constraints = "no_extra_combination_constraints"
	unique(constraints)
}

invariants_for_entry = function(method_name) {
	nm = method_name %||% ""
	invariants = "no_unexpected_error"
	if (grepl("pval|p_value", nm, ignore.case = TRUE)) {
		invariants = c(invariants, "pval_in_0_1")
	}
	if (grepl("confidence_interval|\\bci\\b", nm, ignore.case = TRUE)) {
		invariants = c(invariants, "ci_length_2", "ci_ordered", "ci_numeric_or_documented_missing")
	}
	if (grepl("estimate", nm, ignore.case = TRUE)) {
		invariants = c(invariants, "estimate_shape_valid")
	}
	if (grepl("distribution|bootstrap|randomization|rand", nm, ignore.case = TRUE)) {
		invariants = c(invariants, "distribution_length_valid")
	}
	unique(invariants)
}

optional_dependencies_for_entry = function(arg_names, value_table) {
	deps = character()
	if ("use_rcpp" %in% arg_names && any(value_table$value_expr == "FALSE")) {
		deps = c(deps, "optional_non_rcpp_backend")
	}
	if (!length(deps)) deps = "none"
	unique(deps)
}

slow_path_labels_for_entry = function(arg_names, value_table) {
	labels = character()
	if (any(arg_names %in% c("B", "r", "nsim", "num_boot"))) labels = c(labels, "resampling_count")
	if (any(grepl("root|iteration|maxit", arg_names, ignore.case = TRUE))) labels = c(labels, "iteration_or_root_search")
	if (any(value_table$tier %in% c("nightly", "release"))) labels = c(labels, "higher_tier_values")
	if (!length(labels)) labels = "none"
	unique(labels)
}

unsupported_contexts_for_entry = function(value_table) {
	ctx = character()
	if ("assertNoCensoring" %in% value_table$contract_assertion) ctx = c(ctx, "censored_data")
	if (any(grepl("fixture:support", value_table$value_expr, fixed = TRUE))) ctx = c(ctx, "fixture_support_required")
	if (!length(ctx)) ctx = "none"
	unique(ctx)
}

build_entry = function(target_contracts) {
	api_name = target_contracts$api_name[1]
	method_name = target_contracts$method_name[1]
	arg_names = sort(unique(target_contracts$arg[nzchar(target_contracts$arg)]))
	domains = lapply(arg_names, function(arg) {
		combine_arg_domains(target_contracts[target_contracts$arg == arg, , drop = FALSE])
	})
	names(domains) = arg_names
	value_table = do.call(rbind, Map(function(arg, dom) {
		dom$arg = arg
		dom
	}, names(domains), domains))
	value_table = value_table[, c("arg", setdiff(names(value_table), "arg")), drop = FALSE]
	list(
		target = target_key(api_name, method_name),
		api_name = api_name,
		api_kind = target_contracts$api_kind[1],
		class_name = target_contracts$class_name[1],
		method_name = method_name,
		domains = domains,
		constraints = constraints_for_entry(arg_names, value_table),
		invariants = invariants_for_entry(method_name),
		tier_metadata = list(
			available_tiers = sort(unique(value_table$tier)),
			min_tier = if ("smoke" %in% value_table$tier) "smoke" else sort(unique(value_table$tier))[1]
		),
		optional_dependencies = optional_dependencies_for_entry(arg_names, value_table),
		slow_path_labels = slow_path_labels_for_entry(arg_names, value_table),
		unsupported_contexts = unsupported_contexts_for_entry(value_table),
		source_contract_count = nrow(target_contracts),
		source_contract_lines = sort(unique(target_contracts$source_line))
	)
}

flatten_registry = function(registry) {
	rows = list()
	for (entry in registry) {
		for (arg in names(entry$domains)) {
			dom = entry$domains[[arg]]
			for (i in seq_len(nrow(dom))) {
				rows[[length(rows) + 1L]] = data.frame(
					target = entry$target,
					api_name = entry$api_name,
					api_kind = entry$api_kind,
					class_name = entry$class_name,
					method_name = entry$method_name,
					arg = arg,
					label = dom$label[i],
					value_expr = dom$value_expr[i],
					tier = dom$tier[i],
					source = dom$source[i],
					rationale = dom$rationale[i],
					contract_assertion = dom$contract_assertion[i],
					constraints = paste(entry$constraints, collapse = ";"),
					invariants = paste(entry$invariants, collapse = ";"),
					optional_dependencies = paste(entry$optional_dependencies, collapse = ";"),
					slow_path_labels = paste(entry$slow_path_labels, collapse = ";"),
					unsupported_contexts = paste(entry$unsupported_contexts, collapse = ";"),
					source_contract_count = entry$source_contract_count,
					stringsAsFactors = FALSE
				)
			}
		}
	}
	if (!length(rows)) return(data.frame())
	do.call(rbind, rows)
}

high_priority_methods = c(
	"compute_estimate",
	"compute_asymp_two_sided_pval",
	"compute_asymp_confidence_interval",
	"compute_bootstrap_two_sided_pval",
	"compute_bootstrap_confidence_interval",
	"compute_bayesian_bootstrap_two_sided_pval",
	"compute_bayesian_bootstrap_confidence_interval",
	"compute_rand_two_sided_pval",
	"compute_rand_confidence_interval"
)

validate_registry = function(registry, flat) {
	if (!length(registry)) stop("Registry is empty.", call. = FALSE)
	required_entry_fields = c("domains", "constraints", "invariants", "tier_metadata")
	bad = names(registry)[vapply(registry, function(entry) any(!required_entry_fields %in% names(entry)), logical(1))]
	if (length(bad)) stop("Registry entries missing required fields: ", paste(head(bad), collapse = ", "), call. = FALSE)
	if (!nrow(flat)) stop("Flattened registry is empty.", call. = FALSE)
	hp_present = unique(flat$method_name[flat$method_name %in% high_priority_methods])
	if (length(hp_present) < 5L) {
		stop("Too few high-priority public inference methods represented in registry: ", paste(hp_present, collapse = ", "), call. = FALSE)
	}
	invisible(TRUE)
}

main = function() {
	contracts_path = file.path(repo_root, "package_tests", "checkmate_argument_contracts.csv")
	if (!file.exists(contracts_path)) {
		stop("Missing package_tests/checkmate_argument_contracts.csv. Run Phase 1 first.", call. = FALSE)
	}
	contracts = read.csv(contracts_path, stringsAsFactors = FALSE, na.strings = character())
	contracts = contracts[nzchar(contracts$arg), , drop = FALSE]
	contracts$key = target_key(contracts$api_name, contracts$method_name)
	groups = split(contracts, contracts$key)
	message("Building registry entries for ", length(groups), " public API targets.")
	registry = lapply(groups, build_entry)
	names(registry) = names(groups)
	message("Flattening registry domains.")
	flat = flatten_registry(registry)
	flat = flat[order(flat$target, flat$arg, flat$tier, flat$label), , drop = FALSE]
	row.names(flat) = NULL
	validate_registry(registry, flat)
	saveRDS(registry, file.path(repo_root, "package_tests", "public_argument_contract_registry.rds"))
	write.csv(flat, file.path(repo_root, "package_tests", "public_argument_contract_registry.csv"), row.names = FALSE)
	message("Wrote ", length(registry), " registry entries.")
	message("Wrote ", nrow(flat), " flattened registry value rows.")
	message("High-priority methods represented: ", paste(sort(unique(flat$method_name[flat$method_name %in% high_priority_methods])), collapse = ", "))
}

if (identical(environment(), globalenv())) {
	main()
}
