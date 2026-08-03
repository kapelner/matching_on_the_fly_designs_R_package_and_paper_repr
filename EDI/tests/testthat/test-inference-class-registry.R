inference_registry_source_files = function() {
	files = Sys.glob(file.path("EDI", "R", "*"))
	if (length(files) == 0L) {
		files = Sys.glob(file.path("R", "*"))
	}
	files[file.info(files)$isdir %in% FALSE]
}

canonical_inference_generators = function() {
	ns = asNamespace("EDI")
	all_names = sort(ls(ns))
	Filter(function(name) {
		obj = get(name, envir = ns)
		EDI:::is_inference_r6_generator(obj) && identical(obj$classname, name)
	}, all_names)
}

temporary_inference_metadata = function(likelihood_tier = "none") {
	list(
		abstract = FALSE,
		exported = FALSE,
		response_types = "continuous",
		design_families = "all",
		compatibility = EDI:::always_compatible_inference_metadata,
		likelihood_tier = likelihood_tier,
		required_packages = character()
	)
}

test_that("inference class registry has one canonical metadata entry per generator", {
	EDI:::populate_inference_class_registry()
	registry = EDI:::inference_class_registry_as_list()
	generators = canonical_inference_generators()

	expect_setequal(names(registry), generators)
	expect_equal(length(registry), length(generators))
	expect_silent(EDI:::validate_inference_class_registry(registry))

	for (name in generators) {
		metadata = EDI:::get_inference_class_metadata(name)
		gen = get(name, envir = asNamespace("EDI"))
		parent_gen = gen$get_inherit()
		parent_name = if (is.null(parent_gen)) NULL else parent_gen$classname

		expect_identical(metadata$name, name)
		expect_identical(metadata$parent, parent_name)
		expect_true(is.logical(metadata$abstract))
		expect_true(is.logical(metadata$exported))
		expect_true(metadata$likelihood_tier %in% EDI:::EDI_INFERENCE_ALLOWED_LIKELIHOOD_TIERS)
		expect_identical(EDI:::get_direct_components(name), metadata$direct_components)
		expect_type(EDI:::get_effective_components(name), "character")
		expect_type(EDI:::get_effective_capabilities(name), "character")
	}
})

test_that("inference class metadata rejects missing fields and invalid tiers", {
	valid = EDI:::get_inference_class_metadata("Inference")
	expect_silent(EDI:::validate_inference_class_metadata(valid))

	missing_tier = valid
	missing_tier$likelihood_tier = NULL
	expect_error(
		EDI:::validate_inference_class_metadata(missing_tier),
		"missing required field"
	)

	bad_tier = valid
	bad_tier$likelihood_tier = "pseudo"
	expect_error(
		EDI:::validate_inference_class_metadata(bad_tier),
		"invalid `likelihood_tier`"
	)

	with_alias = valid
	with_alias$alias_of = "OtherInference"
	expect_error(
		EDI:::validate_inference_class_metadata(with_alias),
		"forbidden alias metadata"
	)
})

test_that("top-level Inference aliases do not exist in source", {
	files = inference_registry_source_files()
	alias_pattern = "^\\s*Inference[[:alnum:]_.]+\\s*(<-|=)\\s*Inference[[:alnum:]_.]+\\s*(#.*)?$"
	offenders = unlist(lapply(files, function(file) {
		lines = readLines(file, warn = FALSE)
		hits = grep(alias_pattern, lines)
		if (length(hits) == 0L) return(character())
		sprintf("%s:%s:%s", file, hits, trimws(lines[hits]))
	}), use.names = FALSE)

	if (length(offenders) > 0L) {
		stop(paste(offenders, collapse = "\n"), call. = FALSE)
	}
	expect_equal(length(offenders), 0L)
})

test_that("resolve_inference_components expands dependencies deterministically", {
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)
	EDI:::register_inference_class(
		name = "InferenceTemporaryComponentDependencyHost",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "InferenceMixinKKPassThroughCompound"
	)

	expect_identical(
		EDI:::resolve_inference_components("InferenceTemporaryComponentDependencyHost"),
		c("InferenceMixinKKPassThrough", "InferenceMixinKKPassThroughCompound")
	)
	expect_identical(
		EDI:::get_effective_capabilities("InferenceTemporaryComponentDependencyHost"),
		c("kk_passthrough", "nonparametric_bootstrap", "kk_compound")
	)
})

test_that("daughter classes inherit parent components and capabilities", {
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)
	EDI:::register_inference_class(
		name = "InferenceTemporaryParentComponentHost",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "InferenceMixinKKPassThrough"
	)
	EDI:::register_inference_class(
		name = "InferenceTemporaryChildComponentHost",
		parent = "InferenceTemporaryParentComponentHost",
		metadata = temporary_inference_metadata(likelihood_tier = "full"),
		direct_components = "InferenceMixinOffOptimumLikelihoodEval"
	)

	expect_identical(
		EDI:::get_direct_components("InferenceTemporaryChildComponentHost"),
		"InferenceMixinOffOptimumLikelihoodEval"
	)
	expect_identical(
		EDI:::resolve_inference_components("InferenceTemporaryChildComponentHost"),
		c("InferenceMixinKKPassThrough", "InferenceMixinOffOptimumLikelihoodEval")
	)
	expect_true("kk_passthrough" %in% EDI:::get_effective_capabilities("InferenceTemporaryChildComponentHost"))
	expect_true("nonparametric_bootstrap" %in% EDI:::get_effective_capabilities("InferenceTemporaryChildComponentHost"))
	expect_true("off_optimum_likelihood_eval" %in% EDI:::get_effective_capabilities("InferenceTemporaryChildComponentHost"))
})

test_that("parent components satisfy child component dependencies", {
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)
	EDI:::register_inference_class(
		name = "InferenceTemporaryDependencyParent",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "InferenceMixinKKPassThrough"
	)
	EDI:::register_inference_class(
		name = "InferenceTemporaryDependencyChild",
		parent = "InferenceTemporaryDependencyParent",
		metadata = temporary_inference_metadata(),
		direct_components = "InferenceMixinKKPassThroughCompound"
	)

	expect_identical(
		EDI:::resolve_inference_components("InferenceTemporaryDependencyChild"),
		c("InferenceMixinKKPassThrough", "InferenceMixinKKPassThroughCompound")
	)
})

test_that("component resolution rejects inherited re-listing and duplicate transitive components", {
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)
	EDI:::register_inference_class(
		name = "InferenceTemporaryRelistParent",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "InferenceMixinKKPassThrough"
	)
	EDI:::register_inference_class(
		name = "InferenceTemporaryRelistChild",
		parent = "InferenceTemporaryRelistParent",
		metadata = temporary_inference_metadata(),
		direct_components = "InferenceMixinKKPassThrough"
	)
	expect_error(
		EDI:::resolve_inference_components("InferenceTemporaryRelistChild"),
		"re-lists inherited component"
	)

	EDI:::register_inference_class(
		name = "InferenceTemporaryDuplicateTransitiveHost",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = c("InferenceMixinKKPassThrough", "InferenceMixinKKPassThroughCompound")
	)
	expect_error(
		EDI:::resolve_inference_components("InferenceTemporaryDuplicateTransitiveHost"),
		"duplicates transitive dependency"
	)
})

test_that("component resolution rejects dependency cycles", {
	on.exit({
		EDI:::populate_inference_component_registry()
		EDI:::populate_inference_class_registry()
	}, add = TRUE)
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryCycleA",
		status = "active",
		file = "test",
		dependencies = "InferenceTemporaryCycleB",
		allowed_likelihood_tiers = EDI:::EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS
	))
	EDI:::register_inference_component(EDI:::InferenceComponent(
		name = "InferenceTemporaryCycleB",
		status = "active",
		file = "test",
		dependencies = "InferenceTemporaryCycleA",
		allowed_likelihood_tiers = EDI:::EDI_COMPONENT_ALLOWED_LIKELIHOOD_TIERS
	))
	EDI:::register_inference_class(
		name = "InferenceTemporaryCycleHost",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "InferenceTemporaryCycleA"
	)

	expect_error(
		EDI:::resolve_inference_components("InferenceTemporaryCycleHost"),
		"dependency cycle"
	)
})
