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

test_that("every mixin is registered as an inference component", {
	EDI:::populate_inference_component_registry()
	components = EDI:::inference_component_registry_as_list()
	mixin_names = ls(asNamespace("EDI"), pattern = "^InferenceMixin")

	expect_setequal(names(components), mixin_names)
	for (component_name in names(components)) {
		component = components[[component_name]]
		expect_silent(EDI:::validate_inference_component(component))
		expect_named(component, c(
			"name", "status", "source_name", "file", "public", "private",
			"dependencies", "provides_public_methods", "provides_private_methods",
			"owns_state", "requires_state", "requires_public_methods",
			"requires_private_methods", "optional_public_methods",
			"optional_private_methods", "requires_super_methods",
			"requires_capabilities", "provides_capabilities",
			"allowed_likelihood_tiers", "conflicts", "allowed_host_overrides",
			"forbidden_refs"
		))
		expect_true(component$status %in% c("active", "scaffold"))
		expect_true(is.list(component$public))
		expect_true(is.list(component$private))
		expect_true(is.character(component$owns_state))
		expect_true(is.character(component$requires_state))
		expect_true(is.character(component$requires_public_methods))
		expect_true(is.character(component$requires_private_methods))
		expect_true(is.character(component$optional_public_methods))
		expect_true(is.character(component$optional_private_methods))
		expect_true(is.character(component$dependencies))
		expect_true(is.character(component$requires_capabilities))
		expect_true(is.character(component$provides_capabilities))
		expect_true(is.character(component$conflicts))
		expect_true(all(component$allowed_likelihood_tiers %in% EDI:::EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS))
	}
})

test_that("component provided method metadata matches actual list names", {
	EDI:::populate_inference_component_registry()
	components = EDI:::inference_component_registry_as_list()
	for (component in components) {
		expect_identical(
			sort(component$provides_public_methods),
			sort(EDI:::component_public_names(component))
		)
		expect_identical(
			sort(component$provides_private_methods),
			sort(EDI:::component_private_names(component))
		)
	}
})

test_that("scaffold components are registered but cannot be effective class components", {
	EDI:::populate_inference_component_registry()
	components = EDI:::inference_component_registry_as_list()
	scaffolds = names(Filter(function(component) {
		identical(component$status, "scaffold")
	}, components))

	expect_setequal(scaffolds, c(
		"InferenceMixinCordeiroFerrariApprox",
		"InferenceMixinLemonteGradientApprox"
	))
	expect_silent(EDI:::validate_no_scaffold_effective_components())

	metadata = list(
		abstract = FALSE,
		exported = FALSE,
		response_types = "continuous",
		design_families = "all",
		compatibility = EDI:::always_compatible_inference_metadata,
		likelihood_tier = "none",
		required_packages = character()
	)
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)
	EDI:::register_inference_class(
		name = "InferenceTemporaryScaffoldComponentHost",
		parent = "Inference",
		metadata = metadata,
		direct_components = "InferenceMixinCordeiroFerrariApprox"
	)
	expect_error(
		EDI:::get_effective_components("InferenceTemporaryScaffoldComponentHost"),
		"Scaffold component"
	)
})

test_that("component body references are declared by component contracts", {
	EDI:::populate_inference_component_registry()
	components = EDI:::inference_component_registry_as_list()
	for (component in components) {
		expect_silent(EDI:::validate_component_body_references(component))
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

test_that("define_inference_class assembles component public and private members", {
	on.exit(EDI:::populate_inference_component_registry(), add = TRUE)
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryFactoryComponent",
		status = "active",
		file = "test",
		public = list(
			component_public = function() private$component_private()
		),
		private = list(
			component_state = "assembled",
			component_private = function() private$component_state
		),
		owns_state = "component_state",
		provides_capabilities = "temporary_factory_capability"
	))
	gen = EDI:::define_inference_class(
		classname = "InferenceTemporaryFactoryHost",
		components = "InferenceTemporaryFactoryComponent",
		public = list(
			host_public = function() "host"
		),
		private = list(
			host_private = function() "private"
		),
		metadata = list(
			likelihood_tier = "none"
		)
	)
	obj = gen$new()

	expect_false(gen$lock_objects)
	expect_identical(obj$component_public(), "assembled")
	expect_identical(obj$host_public(), "host")
	expect_true("component_public" %in% names(gen$public_methods))
	expect_true("component_private" %in% names(gen$private_methods))
})

test_that("factory validation rejects unsatisfied component contracts", {
	on.exit(EDI:::populate_inference_component_registry(), add = TRUE)
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryNeedsHostContract",
		status = "active",
		file = "test",
		requires_public_methods = "needed_public",
		requires_private_methods = "needed_private",
		requires_state = "needed_state"
	))

	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryMissingPublic",
			components = "InferenceTemporaryNeedsHostContract",
			private = list(
				needed_private = function() TRUE,
				needed_state = TRUE
			)
		),
		"missing public method"
	)
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryMissingPrivate",
			components = "InferenceTemporaryNeedsHostContract",
			public = list(needed_public = function() TRUE),
			private = list(needed_state = TRUE)
		),
		"missing private method"
	)
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryMissingState",
			components = "InferenceTemporaryNeedsHostContract",
			public = list(needed_public = function() TRUE),
			private = list(needed_private = function() TRUE)
		),
		"missing private state"
	)
	expect_silent(EDI:::define_inference_class(
		classname = "InferenceTemporarySatisfiedContract",
		components = "InferenceTemporaryNeedsHostContract",
		public = list(needed_public = function() TRUE),
		private = list(
			needed_private = function() TRUE,
			needed_state = TRUE
		)
	))
})

test_that("factory validation enforces likelihood tiers and capabilities", {
	on.exit(EDI:::populate_inference_component_registry(), add = TRUE)
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryProviderCapability",
		status = "active",
		file = "test",
		provides_capabilities = "temporary_needed_capability"
	))
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryNeedsCapability",
		status = "active",
		file = "test",
		requires_capabilities = "temporary_needed_capability"
	))

	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryWrongTier",
			components = "InferenceMixinOffOptimumLikelihoodEval",
			metadata = list(likelihood_tier = "none")
		),
		"disallowed likelihood tier"
	)
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryMissingCapability",
			components = "InferenceTemporaryNeedsCapability"
		),
		"missing capability"
	)
	expect_silent(EDI:::define_inference_class(
		classname = "InferenceTemporarySatisfiedCapability",
		components = c("InferenceTemporaryProviderCapability", "InferenceTemporaryNeedsCapability")
	))
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryCapabilityMissingMethod",
			components = "InferenceTemporaryProviderCapability",
			public_methods_for_capability = list(
				temporary_needed_capability = "capability_public_method"
			)
		),
		"without public method"
	)
	expect_silent(EDI:::define_inference_class(
		classname = "InferenceTemporaryCapabilityHasMethod",
		components = "InferenceTemporaryProviderCapability",
		public = list(capability_public_method = function() TRUE),
		public_methods_for_capability = list(
			temporary_needed_capability = "capability_public_method"
		)
	))
})

test_that("factory validation rejects collisions unless overrides declare them", {
	on.exit(EDI:::populate_inference_component_registry(), add = TRUE)
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryCollisionA",
		status = "active",
		file = "test",
		public = list(dup_public = function() "a"),
		private = list(dup_private = function() "a")
	))
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryCollisionB",
		status = "active",
		file = "test",
		public = list(dup_public = function() "b"),
		private = list(dup_private = function() "b")
	))
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryStateCollision",
		status = "active",
		file = "test",
		private = list(dup_private = "state")
	))

	expect_error(
		EDI:::assemble_public(
			"BadCollisionHost",
			c("InferenceTemporaryCollisionA", "InferenceTemporaryCollisionB"),
			resolve = FALSE
		),
		"undeclared public component collision"
	)
	expect_silent(EDI:::assemble_public(
		"DeclaredCollisionHost",
		c("InferenceTemporaryCollisionA", "InferenceTemporaryCollisionB"),
		overrides = list(public = "dup_public"),
		resolve = FALSE
	))
	expect_error(
		EDI:::assemble_private(
			"BadStateCollisionHost",
			c("InferenceTemporaryCollisionA", "InferenceTemporaryStateCollision"),
			resolve = FALSE
		),
		"method/state collision"
	)
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryPublicPrivateDup",
			public = list(same_name = function() TRUE),
			private = list(same_name = function() TRUE)
		),
		"public/private name duplication"
	)
	expect_silent(EDI:::validate_inference_class_definition(
		classname = "InferenceTemporaryDeclaredPublicPrivateDup",
		public = list(same_name = function() TRUE),
		private = list(same_name = function() TRUE),
		overrides = list(public_private = "same_name")
	))
	expect_error(
		EDI:::define_inference_class(
			classname = "InferenceTemporaryLocked",
			lock_objects = TRUE
		),
		"lock_objects = FALSE"
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
