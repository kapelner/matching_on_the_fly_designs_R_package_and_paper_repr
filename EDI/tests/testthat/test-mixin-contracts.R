library(testthat)
library(EDI)

mixin_slot_names = function(mixin_names, slot){
	as.character(unlist(lapply(mixin_names, function(mixin_name) {
		mixin = get(mixin_name, envir = asNamespace("EDI"))
		names(mixin[[slot]])
	}), use.names = FALSE))
}

edi_r_source_files = function() {
	files = Sys.glob(file.path("EDI", "R", "*.R"))
	if (length(files) == 0L) {
		files = Sys.glob(file.path("R", "*.R"))
	}
	files
}

test_that("every mixin has a documented host contract and is collated after the registry", {
	contracts = EDI:::EDI_MIXIN_CONTRACTS
	mixin_names = ls(asNamespace("EDI"), pattern = "^InferenceMixin")
	expect_setequal(names(contracts), mixin_names)

	for (mixin_name in names(contracts)) {
		contract = contracts[[mixin_name]]
		expect_named(contract, c("file", "private_methods", "private_state"))
		expect_true(is.character(contract$file) && length(contract$file) == 1L)
		expect_true(is.character(contract$private_methods))
		expect_true(is.character(contract$private_state))
	}

	collate = strsplit(utils::packageDescription("EDI")$Collate, "[[:space:]]+")[[1L]]
	collate = gsub("'", "", collate, fixed = TRUE)
	registry_position = match("mixin_contracts.R", collate)
	expect_true(is.finite(registry_position))
	for (file in vapply(contracts, `[[`, character(1), "file")) {
		expect_gt(match(file, collate), registry_position)
	}
})

test_that("mixin composition has no undocumented method-name collisions", {
	for (target in names(EDI:::EDI_MIXIN_COMPOSITIONS)) {
		mixins = EDI:::EDI_MIXIN_COMPOSITIONS[[target]]
		allowed = EDI:::EDI_MIXIN_ALLOWED_COLLISIONS[[target]]
		if (is.null(allowed)) allowed = list(public = character(), private = character())
		for (slot in c("public", "private")) {
			methods = mixin_slot_names(mixins, slot)
			collisions = sort(unique(methods[duplicated(methods)]))
			expect_setequal(collisions, allowed[[slot]])
		}
	}
})

test_that("mixin composition declares dependencies and intended overrides", {
	expect_silent(EDI:::assert_valid_mixin_composition(
		target_name = "InferenceKKPassThroughCompound",
		mixin_names = EDI:::EDI_MIXIN_COMPOSITIONS$InferenceKKPassThroughCompound,
		public_overrides = c(
			"approximate_bootstrap_distribution_beta_hat_T",
			"compute_estimate_with_bootstrap_weights"
		)
	))
	expect_silent(EDI:::assert_valid_mixin_composition(
		target_name = "InferenceKKPassThroughCompoundNoParamBootstrap",
		mixin_names = EDI:::EDI_MIXIN_COMPOSITIONS$InferenceKKPassThroughCompoundNoParamBootstrap,
		public_overrides = c(
			"approximate_bootstrap_distribution_beta_hat_T",
			"compute_estimate_with_bootstrap_weights"
		)
	))
	expect_error(
		EDI:::assert_valid_mixin_composition(
			target_name = "BadCompound",
			mixin_names = c("InferenceMixinKKPassThrough", "InferenceMixinKKPassThrough")
		),
		"duplicate mixin"
	)
	expect_error(
		EDI:::assert_valid_mixin_composition(
			target_name = "BadCompound",
			mixin_names = "InferenceMixinKKPassThroughCompound"
		),
		"without required component"
	)
	expect_error(
		EDI:::assert_valid_mixin_composition(
			target_name = "BadCompound",
			mixin_names = c("InferenceMixinKKPassThrough", "InferenceMixinKKPassThroughCompound")
		),
		"undeclared private mixin collision"
	)
	expect_error(
		EDI:::assert_valid_mixin_composition(
			target_name = "BadCompound",
			mixin_names = "InferenceMixinKKPassThrough",
			public_overrides = "approximate_bootstrap_distribution_beta_hat_T"
		),
		"without declaration"
	)
})

test_that("compound descendants do not re-splice KK pass-through mixins", {
	files = edi_r_source_files()
	records = list()
	for (file in files) {
		lines = readLines(file, warn = FALSE)
		starts = grep("R6::R6Class", lines, fixed = TRUE)
		for (i in seq_along(starts)) {
			start = starts[[i]]
			end = if (i < length(starts)) starts[[i + 1L]] - 1L else length(lines)
			block = paste(lines[start:end], collapse = "\n")
			class_name = sub(".*R6Class\\(\"([^\"]+)\".*", "\\1", lines[[start]])
			inherit_line = grep("inherit =", lines[start:end], value = TRUE)[1L]
			inherit_name = if (length(inherit_line) == 0L || is.na(inherit_line)) {
				NA_character_
			} else {
				sub(".*inherit = ([A-Za-z0-9_]+).*", "\\1", inherit_line)
			}
			records[[class_name]] = list(
				file = basename(file),
				inherit = inherit_name,
				splices_pass_through =
					grepl("public = .*InferenceMixinKKPassThrough\\$public", block) ||
					grepl("private = .*InferenceMixinKKPassThrough\\$private", block)
			)
		}
	}
	inherits_compound_base = function(class_name) {
		seen = character()
		while (!is.na(class_name) && !(class_name %in% seen)) {
			if (class_name %in% c("InferenceKKPassThroughCompound", "InferenceKKPassThroughCompoundNoParamBootstrap")) {
				return(TRUE)
			}
			seen = c(seen, class_name)
			class_name = records[[class_name]]$inherit
			if (is.null(class_name)) return(FALSE)
		}
		FALSE
	}
	offenders = names(records)[vapply(names(records), function(class_name) {
		isTRUE(records[[class_name]]$splices_pass_through) && inherits_compound_base(class_name)
	}, logical(1L))]
	if (is.null(offenders)) offenders = character()
	expect_equal(offenders, character())
})

test_that("R6 generator private lists are not accessed from generator objects", {
	files = edi_r_source_files()
	generator_names = character()
	source_lines = list()
	for (file in files) {
		lines = readLines(file, warn = FALSE)
		source_lines[[file]] = lines
		generator_lines = grep("R6::R6Class", lines, fixed = TRUE, value = TRUE)
		generator_names = c(
			generator_names,
			sub("^\\s*([A-Za-z][A-Za-z0-9_.]*)\\s*(=|<-).*", "\\1", generator_lines)
		)
	}
	generator_names = sort(unique(generator_names))
	offenders = character()
	for (file in names(source_lines)) {
		lines = source_lines[[file]]
		code_lines = ifelse(grepl("^\\s*#", lines), "", lines)
		for (generator_name in generator_names) {
			pattern = paste0("(^|[^A-Za-z0-9_.])", generator_name, "\\$private([^A-Za-z0-9_]|$)")
			hits = grep(pattern, code_lines)
			if (length(hits) > 0L) {
				offenders = c(
					offenders,
					paste0(basename(file), ":", hits, ": ", trimws(lines[hits]))
				)
			}
		}
	}
	expect_equal(offenders, character())
})

test_that("KK OLS IVWC uses the compound bootstrap-weight estimator", {
	expect_null(EDI:::InferenceContinKKOLSIVWC$public_methods$compute_estimate_with_bootstrap_weights)
	expect_true(is.function(EDI:::InferenceKKPassThroughCompoundNoParamBootstrap$public_methods$compute_estimate_with_bootstrap_weights))
	expect_false(identical(
		body(EDI:::InferenceKKPassThroughCompoundNoParamBootstrap$public_methods$compute_estimate_with_bootstrap_weights),
		body(EDI:::InferenceMixinKKPassThrough$public$compute_estimate_with_bootstrap_weights)
	))
})

test_that("single-host protected bases are not decomposed into new mixins", {
	expect_false("InferenceCountZeroAugmentedPoissonAbstract" %in% names(EDI:::EDI_MIXIN_COMPOSITIONS))
	expect_false("InferenceMixinHurdlePoissonClosedForm" %in% names(EDI:::EDI_MIXIN_CONTRACTS))
	expect_false(exists("InferenceMixinHurdlePoissonClosedForm", envir = asNamespace("EDI"), inherits = FALSE))
})
