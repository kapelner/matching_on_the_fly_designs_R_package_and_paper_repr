#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

script_path = sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "package_tests/audit_comprehensive_suite_baseline.R")
repo_root = normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "EDI"))) {
	repo_root = normalizePath(getwd(), mustWork = TRUE)
}

artifact = function(name) file.path(repo_root, "package_tests", name)

paths = list(
	inventory = artifact("public_api_inventory.csv"),
	arg_coverage = artifact("public_argument_combination_coverage.csv"),
	baseline_audit = artifact("comprehensive_suite_baseline_audit.csv"),
	exemptions = artifact("comprehensive_suite_exemptions.csv"),
	comprehensive = artifact("comprehensive_tests.R"),
	testthat_dir = file.path(repo_root, "EDI", "tests", "testthat")
)

read_required_csv = function(path, label) {
	if (!file.exists(path)) stop("Missing ", label, ": ", path, call. = FALSE)
	read.csv(path, stringsAsFactors = FALSE, na.strings = character())
}

read_optional_csv = function(path) {
	if (!file.exists(path) || isTRUE(file.info(path)$size == 0)) return(data.frame())
	tryCatch(read.csv(path, stringsAsFactors = FALSE, na.strings = character()), error = function(e) data.frame())
}

read_text = function(path) {
	if (!file.exists(path) || isTRUE(file.info(path)$size == 0)) return("")
	paste(readLines(path, warn = FALSE), collapse = "\n")
}

target_from_inventory = function(inventory) {
	ifelse(nzchar(inventory$method_name), paste(inventory$export_name, inventory$method_name, sep = "::"), inventory$export_name)
}

literal_present = function(text, value) {
	if (!nzchar(value)) return(FALSE)
	grepl(paste0("\\b", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", value), "\\b"), text, perl = TRUE)
}

identifier_tokens = function(text) {
	unique(unlist(regmatches(text, gregexpr("[A-Za-z.][A-Za-z0-9._]*", text, perl = TRUE)), use.names = FALSE))
}

token_file_index = function(files_text) {
	tokens_by_file = lapply(files_text, identifier_tokens)
	all_tokens = unique(unlist(tokens_by_file, use.names = FALSE))
	setNames(lapply(all_tokens, function(tok) {
		names(tokens_by_file)[vapply(tokens_by_file, function(tokens) tok %in% tokens, logical(1))]
	}), all_tokens)
}

token_present = function(token_index, value) {
	nzchar(value) && value %in% names(token_index)
}

token_files = function(token_index, value) {
	if (!token_present(token_index, value)) return(character())
	token_index[[value]]
}

base_inference_class = function(x) {
	sub(" [\\(\\[].*$", "", as.character(x))
}

design_label_to_class = function(x) {
	map = c(
		Bernoulli = "DesignSeqOneByOneBernoulli",
		iBCRD = "DesignSeqOneByOneiBCRD",
		Efron = "DesignSeqOneByOneEfron",
		KK14 = "DesignSeqOneByOneKK14",
		KK21 = "DesignSeqOneByOneKK21",
		KK21stepwise = "DesignSeqOneByOneKK21stepwise",
		SPBR = "DesignSeqOneByOneSPBR",
		PocockSimon = "DesignSeqOneByOnePocockSimon",
		Urn = "DesignSeqOneByOneUrn",
		RandomBlockSize = "DesignSeqOneByOneRandomBlockSize",
		FixedBernoulli = "DesignFixedBernoulli",
		FixediBCRD = "DesignFixediBCRD",
		FixedBlocking = "DesignFixedBlocking",
		FixedCluster = "DesignFixedCluster",
		FixedBlockedCluster = "DesignFixedBlockedCluster",
		FixedBinaryMatch = "DesignFixedBinaryMatch",
		FixedGreedy = "DesignFixedGreedy",
		FixedRerandomization = "DesignFixedRerandomization",
		FixedMatchingGreedy = "DesignFixedMatchingGreedyPairSwitching",
		FixedDOptimal = "DesignFixedDOptimal",
		FixedAOptimal = "DesignFixedAOptimal"
	)
	unname(map[as.character(x)] %||% NA_character_)
}

read_comprehensive_results = function() {
	files = list.files(file.path(repo_root, "package_tests"), pattern = "^comprehensive_tests_results.*[.]csv$", full.names = TRUE)
	rows = lapply(files, function(path) {
		if (!file.exists(path) || isTRUE(file.info(path)$size == 0)) return(NULL)
		df = if (requireNamespace("data.table", quietly = TRUE)) {
			tryCatch(
				as.data.frame(data.table::fread(path, select = c("design", "inference_class", "function_run", "status"), showProgress = FALSE)),
				error = function(e) data.frame()
			)
		} else {
			read_optional_csv(path)
		}
		if (!nrow(df)) return(NULL)
		df$source_file = basename(path)
		df
	})
	rows = rows[!vapply(rows, is.null, logical(1))]
	if (!length(rows)) return(data.frame())
	do.call(rbind, rows)
}

derive_comprehensive_targets = function(results, comprehensive_text) {
	targets = character()
	if (nrow(results)) {
		if (all(c("inference_class", "function_run") %in% names(results))) {
			inf = base_inference_class(results$inference_class)
			method = as.character(results$function_run)
			ok = nzchar(inf) & nzchar(method)
			targets = c(targets, paste(inf[ok], method[ok], sep = "::"), inf[ok], paste(inf[ok], "initialize", sep = "::"))
		}
		if ("design" %in% names(results)) {
			design_classes = design_label_to_class(unique(results$design))
			design_classes = design_classes[!is.na(design_classes) & nzchar(design_classes)]
			targets = c(targets, design_classes, paste(design_classes, "initialize", sep = "::"))
		}
	}
	for (design_class in unique(design_label_to_class(c(
		"Bernoulli", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise", "SPBR", "PocockSimon", "Urn",
		"RandomBlockSize", "FixedBernoulli", "FixediBCRD", "FixedBlocking", "FixedCluster",
		"FixedBlockedCluster", "FixedBinaryMatch", "FixedGreedy", "FixedRerandomization",
		"FixedMatchingGreedy", "FixedDOptimal", "FixedAOptimal"
	)))) {
		if (!is.na(design_class) && literal_present(comprehensive_text, design_class)) {
			targets = c(targets, design_class, paste(design_class, "initialize", sep = "::"))
		}
	}
	unique(targets[nzchar(targets)])
}

public_testthat_coverage = function(inventory, testthat_text_by_file) {
	token_index = token_file_index(testthat_text_by_file)
	vapply(seq_len(nrow(inventory)), function(i) {
		row = inventory[i, , drop = FALSE]
		if (identical(row$api_kind, "function")) {
			return(token_present(token_index, row$export_name))
		}
		if (identical(row$api_kind, "r6_class")) {
			return(token_present(token_index, row$export_name) || token_present(token_index, row$class_name))
		}
		if (identical(row$api_kind, "r6_public_method")) {
			method_files = token_files(token_index, row$method_name)
			if (!length(method_files)) return(FALSE)
			class_files = unique(c(token_files(token_index, row$export_name), token_files(token_index, row$class_name)))
			return(length(intersect(method_files, class_files)) > 0L)
		}
		FALSE
	}, logical(1))
}

classify_internal_testthat_files = function(testthat_text_by_file) {
	if (!length(testthat_text_by_file)) {
		return(data.frame(
			row_type = character(),
			target = character(),
			testthat_file = character(),
			internal_symbols = character(),
			internal_test_classification = character(),
			classification_reason = character(),
			stringsAsFactors = FALSE
		))
	}
	rows = lapply(names(testthat_text_by_file), function(file) {
		text = testthat_text_by_file[[file]]
		internal_calls = unique(unlist(regmatches(text, gregexpr("EDI:::[A-Za-z.][A-Za-z0-9._]*|getFromNamespace\\([^)]+\\)", text, perl = TRUE)), use.names = FALSE))
		if (!length(internal_calls)) return(NULL)
		data.frame(
			row_type = "testthat_internal",
			target = paste0("testthat::", file),
			testthat_file = file,
			internal_symbols = paste(internal_calls, collapse = ";"),
			internal_test_classification = "internal_safety_net",
			classification_reason = "Direct namespace-internal access; retain as internal safety-net coverage, not public contract coverage.",
			stringsAsFactors = FALSE
		)
	})
	rows = rows[!vapply(rows, is.null, logical(1))]
	if (!length(rows)) {
		return(data.frame(
			row_type = character(),
			target = character(),
			testthat_file = character(),
			internal_symbols = character(),
			internal_test_classification = character(),
			classification_reason = character(),
			stringsAsFactors = FALSE
		))
	}
	do.call(rbind, rows)
}

coverage_status = function(has_arg, has_comp, has_testthat) {
	n = as.integer(has_arg) + as.integer(has_comp) + as.integer(has_testthat)
	ifelse(
		n == 0L, "uncovered",
		ifelse(
			n > 1L, "multiple_sources",
			ifelse(has_arg, "argument_combinations", ifelse(has_comp, "comprehensive_workflow", "focused_testthat"))
		)
	)
}

build_baseline_audit = function() {
	inventory = read_required_csv(paths$inventory, "public_api_inventory")
	arg_coverage = read_required_csv(paths$arg_coverage, "public_argument_combination_coverage")
	comprehensive_results = read_comprehensive_results()
	comprehensive_text = read_text(paths$comprehensive)
	testthat_files = list.files(paths$testthat_dir, pattern = "[.]R$", full.names = TRUE)
	testthat_files = testthat_files[!grepl("(^|/)helper-|(^|/)setup-|(^|/)teardown-", testthat_files)]
	testthat_text_by_file = setNames(lapply(testthat_files, read_text), basename(testthat_files))
	
	inventory$target = target_from_inventory(inventory)
	arg_targets = unique(arg_coverage$target[arg_coverage$has_executed_legal_combination %in% c(TRUE, "TRUE") | arg_coverage$has_legal_combination_case %in% c(TRUE, "TRUE")])
	comp_targets = derive_comprehensive_targets(comprehensive_results, comprehensive_text)
	has_comp_text = inventory$reached_by_comprehensive_tests %in% c(TRUE, "TRUE")
	has_comp_results = inventory$target %in% comp_targets
	has_testthat = public_testthat_coverage(inventory, testthat_text_by_file) | inventory$reached_by_one_off_testthat %in% c(TRUE, "TRUE")
	
	audit = data.frame(
		row_type = "public_api",
		target = inventory$target,
		api_kind = inventory$api_kind,
		export_name = inventory$export_name,
		class_name = inventory$class_name,
		method_name = inventory$method_name,
		configurable_arg_count = inventory$configurable_arg_count,
		has_more_than_one_configurable_arg = inventory$has_more_than_one_configurable_arg,
		has_argument_combination_coverage = inventory$target %in% arg_targets,
		has_comprehensive_workflow_coverage = has_comp_text | has_comp_results,
		has_focused_testthat_coverage = has_testthat,
		coverage_status = coverage_status(inventory$target %in% arg_targets, has_comp_text | has_comp_results, has_testthat),
		responsible_runner = "",
		testthat_file = "",
		internal_symbols = "",
		internal_test_classification = "",
		classification_reason = "",
		stringsAsFactors = FALSE
	)
	audit$responsible_runner = ifelse(
		audit$coverage_status == "uncovered", "exemption_required",
		ifelse(
			audit$coverage_status == "argument_combinations", "argument_combinations",
			ifelse(audit$coverage_status == "comprehensive_workflow", "comprehensive_tests", ifelse(audit$coverage_status == "focused_testthat", "testthat", "multiple"))
		)
	)
	internal_rows = classify_internal_testthat_files(testthat_text_by_file)
	if (nrow(internal_rows)) {
		for (col in setdiff(names(audit), names(internal_rows))) internal_rows[[col]] = ""
		internal_rows = internal_rows[, names(audit), drop = FALSE]
		audit = rbind(audit, internal_rows)
	}
	row.names(audit) = NULL
	audit
}

build_exemptions = function(audit) {
	uncovered = audit[audit$row_type == "public_api" & audit$coverage_status == "uncovered", , drop = FALSE]
	exemptions = data.frame(
		target = uncovered$target,
		api_kind = uncovered$api_kind,
		class_name = uncovered$class_name,
		method_name = uncovered$method_name,
		exemption_type = "phase0_uncovered_public_api",
		reason = "Initial Phase 0 baseline exemption: no argument-combination, comprehensive-workflow, or focused testthat coverage detected in current artifacts.",
		expiry_date = "",
		owner = "",
		created_date = as.character(Sys.Date()),
		stringsAsFactors = FALSE
	)
	row.names(exemptions) = NULL
	exemptions
}

main = function() {
	audit = build_baseline_audit()
	exemptions = build_exemptions(audit)
	write.csv(audit, paths$baseline_audit, row.names = FALSE)
	write.csv(exemptions, paths$exemptions, row.names = FALSE)
	message("Wrote baseline audit rows: ", nrow(audit))
	message("Wrote exemption rows: ", nrow(exemptions))
	message("Public API coverage status:")
	print(table(audit$coverage_status[audit$row_type == "public_api"], useNA = "ifany"))
	message("Internal testthat files classified: ", sum(audit$row_type == "testthat_internal"))
	invisible(list(audit = audit, exemptions = exemptions))
}

main()
