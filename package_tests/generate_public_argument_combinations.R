#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/generate_public_argument_combinations.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(repo_root, "package_tests", "public_argument_combination_constraints.R"))
source(file.path(repo_root, "package_tests", "public_argument_combination_fixtures.R"))

generator_tier_rank = c(smoke = 1L, ci = 2L, nightly = 3L, release = 4L)

registry_path = file.path(repo_root, "package_tests", "public_argument_contract_registry.rds")
cases_path = file.path(repo_root, "package_tests", "public_argument_combination_cases.csv")
rejected_path = file.path(repo_root, "package_tests", "public_argument_combination_rejected_candidates.csv")
coverage_path = file.path(repo_root, "package_tests", "public_argument_combination_coverage.csv")

read_registry = function(path = registry_path) {
	if (!file.exists(path)) {
		stop("Missing package_tests/public_argument_contract_registry.rds. Run Phase 2 first.", call. = FALSE)
	}
	readRDS(path)
}

tier_allows = function(value_tier, requested_tier) {
	value_rank = generator_tier_rank[[value_tier]]
	request_rank = generator_tier_rank[[requested_tier]]
	!is.null(value_rank) && !is.null(request_rank) && !is.na(value_rank) && !is.na(request_rank) && value_rank <= request_rank
}

target_response_family = function(entry) {
	txt = paste(entry$target %||% "", entry$class_name %||% "", entry$api_name %||% "", sep = " ")
	if (grepl("Survival|surv", txt, ignore.case = TRUE)) return("survival")
	if (grepl("Ordinal|ordinal", txt, ignore.case = TRUE)) return("ordinal")
	if (grepl("Count|count|Poisson|NegBin|Hurdle", txt, ignore.case = TRUE)) return("count")
	if (grepl("Prop|proportion|Beta", txt, ignore.case = TRUE)) return("proportion")
	if (grepl("Incid|incidence|Logit|Probit|Binomial|Risk", txt, ignore.case = TRUE)) return("incidence")
	if (grepl("Contin|continuous|Robust|Quantile|Mean", txt, ignore.case = TRUE)) return("continuous")
	NA_character_
}

target_design_type = function(entry) {
	txt = paste(entry$target %||% "", entry$class_name %||% "", entry$api_name %||% "", sep = " ")
	if (grepl("BlockedCluster|blocked_cluster", txt, ignore.case = TRUE)) return("blocked_cluster")
	if (grepl("Cluster", txt, ignore.case = TRUE)) return("clustered")
	if (grepl("Blocking|Blocks|CMH|Robins", txt, ignore.case = TRUE)) return("stratified")
	if (grepl("BinaryMatch|Matching", txt, ignore.case = TRUE)) return("matched")
	if (grepl("Greedy|Optimal|AOptimal", txt, ignore.case = TRUE)) return("search")
	if (grepl("SeqOneByOne|KK14|KK21|Efron|Atkinson", txt, ignore.case = TRUE)) return("sequential")
	if (grepl("Fixed", txt, ignore.case = TRUE)) return("fixed")
	NA_character_
}

is_high_priority_target = function(entry) {
	method = entry$method_name %||% ""
	method %in% high_priority_fixture_methods ||
		grepl("bootstrap|bayesian_bootstrap|rand|asymp", method, ignore.case = TRUE)
}

has_extra_constraints = function(entry) {
	constraints = entry$constraints %||% character()
	any(nzchar(constraints) & constraints != "no_extra_combination_constraints")
}

compatible_fixtures = function(entry, fixtures) {
	family = target_response_family(entry)
	design_type = target_design_type(entry)
	hits = fixtures
	if (!is.na(family)) {
		hits = hits[vapply(hits, function(fixture) identical(fixture$metadata$response_type, family), logical(1))]
	}
	if (!is.na(design_type)) {
		design_hits = hits[vapply(hits, function(fixture) identical(fixture$metadata$design_type, design_type), logical(1))]
		if (length(design_hits)) hits = design_hits
	}
	if (!length(hits)) {
		fallback = fixtures[vapply(fixtures, function(fixture) identical(fixture$metadata$fixture_id, "fixed_bernoulli_continuous_smoke"), logical(1))]
		if (length(fallback)) return(fallback)
		return(fixtures[seq_len(min(1L, length(fixtures)))])
	}
	hits[seq_len(min(2L, length(hits)))]
}

domain_values = function(domain, tier, include_any_if_empty = TRUE) {
	if (!nrow(domain)) return(domain)
	allowed = domain[vapply(domain$tier, tier_allows, logical(1), requested_tier = tier), , drop = FALSE]
	if (!nrow(allowed) && isTRUE(include_any_if_empty)) {
		allowed = domain[1L, , drop = FALSE]
	}
	row.names(allowed) = NULL
	allowed
}

default_domain_row = function(domain, tier) {
	if (!nrow(domain)) return(domain)
	default = domain[domain$label == "default", , drop = FALSE]
	if (nrow(default)) return(default[1L, , drop = FALSE])
	allowed = domain_values(domain, tier)
	allowed[1L, , drop = FALSE]
}

non_default_domain_rows = function(domain, tier, max_values = 3L) {
	rows = domain_values(domain, tier, include_any_if_empty = FALSE)
	rows = rows[rows$label != "default", , drop = FALSE]
	if (!nrow(rows)) return(rows)
	rows[seq_len(min(max_values, nrow(rows))), , drop = FALSE]
}

eval_value_expr = function(value_expr, fixture) {
	value_expr = trimws(as.character(value_expr))
	if (!nzchar(value_expr)) return(NULL)
	if (value_expr == "NULL") return(NULL)
	if (grepl("^<fixture:", value_expr)) {
		kind = sub("^<fixture:([^>]+)>$", "\\1", value_expr)
		if (kind == "data.frame") return(fixture$data)
		if (kind == "matrix") return(as.matrix(fixture$data))
		if (kind == "list") return(as.list(fixture$data))
		if (kind == "string") return(fixture$fixture_id)
		if (kind == "no_censoring") return(!fixture$metadata$has_censoring)
		return(fixture)
	}
	if (grepl("^~", value_expr)) return(stats::as.formula(value_expr))
	tryCatch(eval(parse(text = value_expr), envir = baseenv()), error = function(e) value_expr)
}

signature_from_exprs = function(exprs) {
	if (!length(exprs)) return("<default>")
	paste(paste(names(exprs), exprs, sep = "="), collapse = ";")
}

stable_checksum = function(x) {
	ints = utf8ToInt(enc2utf8(x))
	if (!length(ints)) return("00000000")
	acc = 0
	for (i in seq_along(ints)) {
		acc = (acc * 33 + ints[i]) %% 4294967291
	}
	sprintf("%08x", as.integer(acc %% 2147483647))
}

case_id_for = function(entry, fixture_id, coverage_kind, signature) {
	prefix = gsub("[^A-Za-z0-9]+", "_", paste(entry$target, fixture_id, coverage_kind, sep = "__"))
	prefix = gsub("_+", "_", prefix)
	prefix = substr(prefix, 1L, 120L)
	paste(prefix, stable_checksum(signature), sep = "_")
}

case_candidate = function(entry, fixture, tier, coverage_kind, value_rows_by_arg) {
	exprs = vapply(value_rows_by_arg, function(row) row$value_expr[1L], character(1))
	labels = vapply(value_rows_by_arg, function(row) row$label[1L], character(1))
	args = lapply(exprs, eval_value_expr, fixture = fixture)
	names(args) = names(exprs)
	signature = signature_from_exprs(exprs)
	constraint = constraint_context_result(args = args, fixture = fixture$metadata, tier = tier, registry_entry = entry)
	data.frame(
		case_id = case_id_for(entry, fixture$fixture_id, coverage_kind, signature),
		api_kind = entry$api_kind %||% "",
		class_name = entry$class_name %||% "",
		function_name = if (nzchar(entry$method_name %||% "")) "" else entry$api_name %||% "",
		method_name = entry$method_name %||% "",
		fixture_id = fixture$fixture_id,
		tier = tier,
		argument_signature = signature,
		argument_labels = signature_from_exprs(labels),
		argument_coverage_kind = coverage_kind,
		source_contracts = paste(entry$source_contract_lines %||% character(), collapse = ";"),
		constraint_status = constraint$status,
		constraint_reason = constraint$reason,
		constraint_sources = paste(unique(constraint$results$source[constraint$results$status != "valid"]), collapse = ";"),
		target = entry$target %||% "",
		invariants = paste(entry$invariants %||% character(), collapse = ";"),
		stringsAsFactors = FALSE
	)
}

default_case_values = function(entry, tier) {
	rows = lapply(entry$domains, default_domain_row, tier = tier)
	names(rows) = names(entry$domains)
	rows
}

one_nondefault_candidates = function(entry, tier) {
	defaults = default_case_values(entry, tier)
	out = list()
	for (arg in names(entry$domains)) {
		rows = non_default_domain_rows(entry$domains[[arg]], tier, max_values = 1L)
		if (!nrow(rows)) next
		vals = defaults
		vals[[arg]] = rows[1L, , drop = FALSE]
		out[[length(out) + 1L]] = vals
	}
	out
}

pairwise_candidates = function(entry, tier, max_values_per_arg = 3L, max_candidates = 24L) {
	if (length(entry$domains) < 2L) return(list())
	defaults = default_case_values(entry, tier)
	values = lapply(entry$domains, non_default_domain_rows, tier = tier, max_values = max_values_per_arg)
	values = values[vapply(values, nrow, integer(1)) > 0L]
	if (length(values) < 2L) return(list())
	out = list()
	arg_pairs = utils::combn(names(values), 2L, simplify = FALSE)
	for (pair in arg_pairs) {
		left = values[[pair[1L]]]
		right = values[[pair[2L]]]
		grid = expand.grid(i = seq_len(nrow(left)), j = seq_len(nrow(right)))
		for (k in seq_len(nrow(grid))) {
			vals = defaults
			vals[[pair[1L]]] = left[grid$i[k], , drop = FALSE]
			vals[[pair[2L]]] = right[grid$j[k], , drop = FALSE]
			out[[length(out) + 1L]] = vals
			if (length(out) >= max_candidates) return(out)
		}
	}
	out
}

targeted_three_way_candidates = function(entry, tier) {
	high_risk_sets = list(
		c("type", "B", "min_number_usable_samples"),
		c("type", "B", "show_progress"),
		c("model_formula", "transform_responses", "type"),
		c("strata_cols", "cluster_col", "B")
	)
	defaults = default_case_values(entry, tier)
	out = list()
	for (arg_set in high_risk_sets) {
		if (!all(arg_set %in% names(entry$domains))) next
		vals = defaults
		for (arg in arg_set) {
			rows = non_default_domain_rows(entry$domains[[arg]], tier, max_values = 1L)
			if (nrow(rows)) vals[[arg]] = rows[1L, , drop = FALSE]
		}
		out[[length(out) + 1L]] = vals
	}
	out
}

tiny_exhaustive_candidates = function(entry, tier, max_cases = 16L) {
	if (length(entry$domains) == 0L || length(entry$domains) > 3L) return(list())
	values = lapply(entry$domains, domain_values, tier = tier, include_any_if_empty = FALSE)
	if (any(vapply(values, nrow, integer(1)) == 0L)) return(list())
	total = prod(vapply(values, nrow, integer(1)))
	if (total > max_cases) return(list())
	grid = expand.grid(lapply(values, function(x) seq_len(nrow(x))))
	out = vector("list", nrow(grid))
	for (i in seq_len(nrow(grid))) {
		vals = Map(function(domain, idx) domain[idx, , drop = FALSE], values, grid[i, ])
		names(vals) = names(values)
		out[[i]] = vals
	}
	out
}

dedupe_cases = function(cases) {
	if (!nrow(cases)) return(cases)
	keys = paste(cases$case_id, cases$argument_coverage_kind, sep = "\r")
	cases[!duplicated(keys), , drop = FALSE]
}

generate_cases_for_entry = function(entry, fixtures, tier) {
	fixtures = compatible_fixtures(entry, fixtures)
	rows = list()
	expand_combinations = has_extra_constraints(entry) || is_high_priority_target(entry)
	for (fixture in fixtures) {
		rows[[length(rows) + 1L]] = case_candidate(entry, fixture, tier, "default", default_case_values(entry, tier))
		for (vals in one_nondefault_candidates(entry, tier)) {
			rows[[length(rows) + 1L]] = case_candidate(entry, fixture, tier, "one_non_default", vals)
		}
		if (expand_combinations) {
			for (vals in pairwise_candidates(entry, tier)) {
				rows[[length(rows) + 1L]] = case_candidate(entry, fixture, tier, "pairwise", vals)
			}
			for (vals in targeted_three_way_candidates(entry, tier)) {
				rows[[length(rows) + 1L]] = case_candidate(entry, fixture, tier, "targeted_3way", vals)
			}
			for (vals in tiny_exhaustive_candidates(entry, tier)) {
				rows[[length(rows) + 1L]] = case_candidate(entry, fixture, tier, "tiny_exhaustive", vals)
			}
		}
	}
	dedupe_cases(do.call(rbind, rows))
}

generate_public_argument_combinations = function(tier = "smoke", target_filter = NULL, fixture_ids = NULL) {
	registry = read_registry()
	if (!is.null(target_filter)) {
		keep = grepl(target_filter, names(registry)) |
			vapply(registry, function(entry) grepl(target_filter, entry$target %||% ""), logical(1))
		registry = registry[keep]
	}
	if (!length(registry)) stop("No registry entries selected.", call. = FALSE)
	fixtures = build_public_argument_fixtures(tier = tier, fixture_ids = fixture_ids)
	rows = lapply(registry, generate_cases_for_entry, fixtures = fixtures, tier = tier)
	cases = dedupe_cases(do.call(rbind, rows))
	cases = cases[order(cases$target, cases$fixture_id, cases$argument_coverage_kind, cases$argument_signature), , drop = FALSE]
	row.names(cases) = NULL
	cases
}

coverage_from_cases = function(cases) {
	valid = cases[cases$constraint_status == "valid", , drop = FALSE]
	pairwise = valid[valid$argument_coverage_kind %in% c("pairwise", "targeted_3way", "tiny_exhaustive"), , drop = FALSE]
	data.frame(
		target = sort(unique(cases$target)),
		total_cases = as.integer(tabulate(match(cases$target, sort(unique(cases$target))))),
		valid_cases = as.integer(tabulate(match(valid$target, sort(unique(cases$target))), nbins = length(sort(unique(cases$target))))),
		valid_pairwise_or_higher_cases = as.integer(tabulate(match(pairwise$target, sort(unique(cases$target))), nbins = length(sort(unique(cases$target))))),
		stringsAsFactors = FALSE
	)
}

write_generation_outputs = function(cases, cases_file = cases_path, rejected_file = rejected_path, coverage_file = coverage_path) {
	write.csv(cases, cases_file, row.names = FALSE)
	rejected = cases[cases$constraint_status != "valid", , drop = FALSE]
	write.csv(rejected, rejected_file, row.names = FALSE)
	write.csv(coverage_from_cases(cases), coverage_file, row.names = FALSE)
	invisible(list(cases = cases_file, rejected = rejected_file, coverage = coverage_file))
}

main = function() {
	args = commandArgs(trailingOnly = TRUE)
	tier = args[1L] %||% "smoke"
	target_filter = args[2L] %||% NULL
	cases = generate_public_argument_combinations(tier = tier, target_filter = target_filter)
	write_generation_outputs(cases)
	message("Wrote ", nrow(cases), " public argument combination candidates.")
	message("Valid candidates: ", sum(cases$constraint_status == "valid"))
	message("Rejected candidates: ", sum(cases$constraint_status != "valid"))
}

if (sys.nframe() == 0L) {
	main()
}
