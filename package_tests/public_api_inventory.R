#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/public_api_inventory.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

load_edi_namespace = function() {
	if (requireNamespace("pkgload", quietly = TRUE)) {
		pkgload::load_all(file.path(repo_root, "EDI"), quiet = TRUE, export_all = FALSE)
	} else {
		stop("Package 'pkgload' is required for public API inventory generation.", call. = FALSE)
	}
	asNamespace("EDI")
}

safe_deparse = function(x) {
	paste(deparse(x, width.cutoff = 500L), collapse = "\n")
}

collapse_chr = function(x) {
	x = unique(as.character(x))
	x = x[nzchar(x)]
	if (length(x) == 0L) "" else paste(sort(x), collapse = ";")
}

formal_names = function(fn) {
	if (!is.function(fn)) return(character())
	setdiff(names(formals(fn)), "...")
}

configurable_arg_count = function(fn) {
	length(formal_names(fn))
}

is_r6_generator = function(x) {
	inherits(x, "R6ClassGenerator")
}

r6_public_method_rows = function(export_name, generator) {
	rows = list()
	seen = character()
	current = generator
	level = 0L
	while (!is.null(current)) {
		methods = current$public_methods %||% list()
		for (method_name in names(methods)) {
			if (method_name %in% seen) next
			seen = c(seen, method_name)
			fn = methods[[method_name]]
			cache_key = paste(current$classname %||% "", method_name, sep = "::")
			if (exists(cache_key, envir = .method_meta_cache, inherits = FALSE)) {
				meta = get(cache_key, envir = .method_meta_cache, inherits = FALSE)
			} else {
				meta = list(
					configurable_arg_count = configurable_arg_count(fn),
					formal_args = collapse_chr(formal_names(fn)),
					has_assertion_calls = grepl("\\b(assert[A-Za-z0-9_.]*|check[A-Za-z0-9_.]*)\\s*\\(", safe_deparse(fn))
				)
				assign(cache_key, meta, envir = .method_meta_cache)
			}
			rows[[length(rows) + 1L]] = data.frame(
				export_name = export_name,
				api_kind = "r6_public_method",
				class_name = generator$classname %||% export_name,
				method_name = method_name,
				method_origin_class = current$classname %||% "",
				method_inheritance = if (level == 0L) "class_local" else "inherited",
				configurable_arg_count = meta$configurable_arg_count,
				formal_args = meta$formal_args,
				has_assertion_calls = meta$has_assertion_calls,
				stringsAsFactors = FALSE
			)
		}
		current = current$get_inherit()
		level = level + 1L
	}
	do.call(rbind, rows)
}

.method_meta_cache = new.env(parent = emptyenv())

object_kind = function(obj) {
	if (is_r6_generator(obj)) return("r6_class")
	if (is.function(obj)) return("function")
	if (is.data.frame(obj) || is.matrix(obj) || is.atomic(obj) || is.list(obj)) return("data")
	"other"
}

literal_present = function(text, literal) {
	if (!nzchar(literal)) return(FALSE)
	grepl(literal, text, fixed = TRUE)
}

cached_literal_present = function(text, literal, cache) {
	if (!nzchar(literal)) return(FALSE)
	if (exists(literal, envir = cache, inherits = FALSE)) {
		return(get(literal, envir = cache, inherits = FALSE))
	}
	value = literal_present(text, literal)
	assign(literal, value, envir = cache)
	value
}

find_source_files = function(symbol, file_lines) {
	if (length(file_lines) == 0L) return("")
	hits = vapply(file_lines, function(lines) {
		any(grepl(symbol, lines, fixed = TRUE))
	}, logical(1))
	collapse_chr(sub(paste0("^", repo_root, "/?"), "", names(file_lines)[hits]))
}

scan_reference_text = function(text, names_to_find) {
	if (!nzchar(text)) {
		return(setNames(rep(FALSE, length(names_to_find)), names_to_find))
	}
	setNames(vapply(names_to_find, function(nm) literal_present(text, nm), logical(1)), names_to_find)
}

main = function() {
	ns = load_edi_namespace()
	exports = sort(getNamespaceExports("EDI"))
	testthat_files = list.files(file.path(repo_root, "EDI", "tests", "testthat"), pattern = "[.]R$", full.names = TRUE)
	testthat_files = testthat_files[!grepl("(^|/)helper-|(^|/)setup-|(^|/)teardown-", testthat_files)]
	testthat_text = if (length(testthat_files)) paste(unlist(lapply(testthat_files, readLines, warn = FALSE), use.names = FALSE), collapse = "\n") else ""
	comprehensive_file = file.path(repo_root, "package_tests", "comprehensive_tests.R")
	comprehensive_text = if (file.exists(comprehensive_file)) paste(readLines(comprehensive_file, warn = FALSE), collapse = "\n") else ""
	testthat_refs = scan_reference_text(testthat_text, exports)
	comprehensive_ref_cache = new.env(parent = emptyenv())
	testthat_ref_cache = new.env(parent = emptyenv())
	source_file_map = setNames(rep("", length(exports)), exports)
	message("Loaded ", length(exports), " exported symbols.")
	
	export_rows = lapply(exports, function(nm) {
		obj = get(nm, envir = ns)
		fn = if (is.function(obj) && !is_r6_generator(obj)) obj else NULL
		body_text = if (is.null(fn)) "" else safe_deparse(fn)
		data.frame(
			export_name = nm,
			api_kind = object_kind(obj),
			class_name = if (is_r6_generator(obj)) obj$classname %||% nm else "",
			method_name = "",
			method_origin_class = "",
			method_inheritance = "",
			configurable_arg_count = configurable_arg_count(fn),
			formal_args = collapse_chr(formal_names(fn)),
			has_assertion_calls = if (is.null(fn)) FALSE else grepl("\\b(assert[A-Za-z0-9_.]*|check[A-Za-z0-9_.]*)\\s*\\(", body_text),
			has_more_than_one_configurable_arg = configurable_arg_count(fn) > 1L,
			reached_by_comprehensive_tests = cached_literal_present(comprehensive_text, nm, comprehensive_ref_cache),
			reached_by_one_off_testthat = isTRUE(testthat_refs[[nm]]) && !cached_literal_present(comprehensive_text, nm, comprehensive_ref_cache),
			source_files = source_file_map[[nm]] %||% "",
			stringsAsFactors = FALSE
		)
	})
	inventory = do.call(rbind, export_rows)
	
	r6_exports = exports[vapply(exports, function(nm) is_r6_generator(get(nm, envir = ns)), logical(1))]
	message("Enumerating public methods for ", length(r6_exports), " exported R6 classes.")
	method_rows = lapply(r6_exports, function(nm) {
		rows = r6_public_method_rows(nm, get(nm, envir = ns))
		rows$has_more_than_one_configurable_arg = rows$configurable_arg_count > 1L
		rows$reached_by_comprehensive_tests = vapply(
			seq_len(nrow(rows)),
			function(i) cached_literal_present(comprehensive_text, rows$export_name[i], comprehensive_ref_cache) ||
				cached_literal_present(comprehensive_text, rows$method_name[i], comprehensive_ref_cache),
			logical(1)
		)
		rows$reached_by_one_off_testthat = vapply(
			seq_len(nrow(rows)),
			function(i) {
				method_ref = cached_literal_present(testthat_text, rows$export_name[i], testthat_ref_cache) ||
					cached_literal_present(testthat_text, rows$method_name[i], testthat_ref_cache)
				method_ref && !rows$reached_by_comprehensive_tests[i]
			},
			logical(1)
		)
		rows$source_files = vapply(rows$export_name, function(nm) source_file_map[[nm]] %||% "", character(1))
		rows
	})
	if (length(method_rows)) {
		inventory = rbind(inventory, do.call(rbind, method_rows))
	}
	
	inventory = inventory[order(inventory$export_name, inventory$api_kind, inventory$method_name), ]
	row.names(inventory) = NULL
	
	# Phase 0 predates the legal-combination runner, so existing combination coverage
	# is imported if present; otherwise every multi-argument API is treated as a gap.
	combination_coverage_path = file.path(repo_root, "package_tests", "public_argument_combination_coverage.csv")
	if (file.exists(combination_coverage_path)) {
		comb_cov = read.csv(combination_coverage_path, stringsAsFactors = FALSE)
		covered_keys = unique(paste(comb_cov$export_name %||% comb_cov$class_name, comb_cov$method_name %||% "", sep = "::"))
	} else {
		covered_keys = character()
	}
	keys = paste(inventory$export_name, inventory$method_name, sep = "::")
	inventory$has_legal_combination_coverage = keys %in% covered_keys
	
	gap_report = inventory[
		inventory$has_more_than_one_configurable_arg &
			!inventory$has_legal_combination_coverage,
		,
		drop = FALSE
	]
	gap_report$gap_reason = ifelse(
		gap_report$has_assertion_calls,
		"multiple_configurable_args_with_unidimensional_assertions_only",
		"multiple_configurable_args_without_detected_assertion_contract"
	)
	gap_report = gap_report[
		order(gap_report$api_kind, gap_report$export_name, gap_report$method_name),
		c(
			"export_name", "api_kind", "class_name", "method_name",
			"method_origin_class", "method_inheritance", "configurable_arg_count",
			"formal_args", "has_assertion_calls", "reached_by_comprehensive_tests",
			"reached_by_one_off_testthat", "has_legal_combination_coverage",
			"gap_reason", "source_files"
		)
	]
	row.names(gap_report) = NULL
	
	write.csv(inventory, file.path(repo_root, "package_tests", "public_api_inventory.csv"), row.names = FALSE)
	write.csv(gap_report, file.path(repo_root, "package_tests", "public_argument_baseline_gap_report.csv"), row.names = FALSE)
	
	message("Wrote ", nrow(inventory), " inventory rows.")
	message("Wrote ", nrow(gap_report), " baseline gap rows.")
}

main()
