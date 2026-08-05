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

temporary_inference_metadata = function(likelihood_tier = "none", capabilities = character()) {
	list(
		abstract = FALSE,
		exported = FALSE,
		response_types = "continuous",
		design_families = "all",
		compatibility = EDI:::always_compatible_inference_metadata,
		likelihood_tier = likelihood_tier,
		required_packages = character(),
		capabilities = capabilities
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

test_that("hierarchy migration manifest covers every current generator", {
	EDI:::populate_inference_class_registry()
	registry = EDI:::inference_class_registry_as_list()
	manifest = EDI:::inference_hierarchy_migration_manifest_as_list()

	expect_setequal(names(manifest), names(registry))
	for (name in names(registry)) {
		record = EDI:::get_inference_hierarchy_migration_record(name)
		metadata = registry[[name]]
		expect_silent(EDI:::validate_inference_hierarchy_migration_record(record))
		expect_identical(record$name, name)
		expect_identical(record$current_parent, metadata$parent)
		expect_identical(record$current_likelihood_tier, metadata$likelihood_tier)
		expect_identical(record$target_likelihood_tier, metadata$likelihood_tier)
		expect_false(any(grepl("^InferenceMixin", record$target_components)))
		expect_true(all(record$target_components %in% ls(EDI:::EDI_INFERENCE_COMPONENTS)))
	}
})

test_that("concrete migration targets are shallow and component-backed", {
	EDI:::populate_inference_class_registry()
	manifest = EDI:::inference_hierarchy_migration_manifest_as_list()
	concrete = Filter(function(record) {
		!isTRUE(record$current_abstract)
	}, manifest)

	expect_gt(length(concrete), 0L)
	for (record in concrete) {
		expect_identical(record$target_parent, "Inference")
		expect_identical(record$target_components, EDI:::get_effective_components(record$name))
		expect_identical(record$target_capabilities, EDI:::get_effective_capabilities(record$name))
	}
})

test_that("remaining algorithmic compatibility descendants are explicitly tracked", {
	EDI:::populate_inference_class_registry()
	manifest = EDI:::inference_hierarchy_migration_manifest_as_list()
	pending = Filter(function(record) {
		identical(record$migration_status, "pending")
	}, manifest)

	expect_gt(length(pending), 0L)
	for (record in pending) {
		expect_gt(length(record$algorithmic_compatibility_ancestors), 0L)
		expect_true(all(
			record$algorithmic_compatibility_ancestors %in%
				EDI:::EDI_INFERENCE_ALGORITHM_COMPATIBILITY_BASES
		))
	}
})

test_that("migrated concrete classes do not descend from algorithmic compatibility bases", {
	EDI:::populate_inference_class_registry()
	manifest = EDI:::inference_hierarchy_migration_manifest_as_list()
	migrated = Filter(function(record) {
		identical(record$migration_status, "migrated")
	}, manifest)

	expect_type(migrated, "list")
	for (record in migrated) {
		expect_identical(record$current_parent, "Inference")
		expect_length(record$algorithmic_compatibility_ancestors, 0L)
	}
})

test_that("mark_inference_class_migrated validates current structure before updating status", {
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)

	EDI:::register_inference_class(
		name = "InferenceTemporaryMigratable",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = character()
	)
	record = EDI:::build_inference_hierarchy_migration_record("InferenceTemporaryMigratable")
	record$migration_status = "pending"
	assign(record$name, record, envir = EDI:::EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)

	expect_silent(EDI:::mark_inference_class_migrated("InferenceTemporaryMigratable"))
	expect_identical(
		EDI:::get_inference_hierarchy_migration_record("InferenceTemporaryMigratable")$migration_status,
		"migrated"
	)

	EDI:::register_inference_class(
		name = "InferenceTemporaryPendingParent",
		parent = "InferenceAsymp",
		metadata = temporary_inference_metadata(),
		direct_components = character()
	)
	pending_parent_record = EDI:::build_inference_hierarchy_migration_record("InferenceTemporaryPendingParent")
	assign(pending_parent_record$name, pending_parent_record, envir = EDI:::EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
	expect_error(
		EDI:::mark_inference_class_migrated("InferenceTemporaryPendingParent"),
		"current parent"
	)
	expect_identical(
		EDI:::get_inference_hierarchy_migration_record("InferenceTemporaryPendingParent")$migration_status,
		"pending"
	)

	EDI:::register_inference_class(
		name = "InferenceTemporaryComponentMismatch",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "RandomizationTest"
	)
	component_mismatch_record = EDI:::build_inference_hierarchy_migration_record("InferenceTemporaryComponentMismatch")
	component_mismatch_record$migration_status = "pending"
	component_mismatch_record$target_components = character()
	assign(component_mismatch_record$name, component_mismatch_record, envir = EDI:::EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
	expect_error(
		EDI:::mark_inference_class_migrated("InferenceTemporaryComponentMismatch"),
		"effective components"
	)
	expect_identical(
		EDI:::get_inference_hierarchy_migration_record("InferenceTemporaryComponentMismatch")$migration_status,
		"pending"
	)

	EDI:::register_inference_class(
		name = "InferenceTemporaryCapabilityMismatch",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "RandomizationTest"
	)
	capability_mismatch_record = EDI:::build_inference_hierarchy_migration_record("InferenceTemporaryCapabilityMismatch")
	capability_mismatch_record$migration_status = "pending"
	capability_mismatch_record$target_capabilities = character()
	assign(capability_mismatch_record$name, capability_mismatch_record, envir = EDI:::EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
	expect_error(
		EDI:::mark_inference_class_migrated("InferenceTemporaryCapabilityMismatch"),
		"effective capabilities"
	)
	expect_identical(
		EDI:::get_inference_hierarchy_migration_record("InferenceTemporaryCapabilityMismatch")$migration_status,
		"pending"
	)

	EDI:::register_inference_class(
		name = "InferenceTemporaryMissingCapabilityMethods",
		parent = "Inference",
		metadata = temporary_inference_metadata(capabilities = "nonparametric_bootstrap"),
		direct_components = character()
	)
	missing_methods_record = EDI:::build_inference_hierarchy_migration_record("InferenceTemporaryMissingCapabilityMethods")
	missing_methods_record$migration_status = "pending"
	assign(missing_methods_record$name, missing_methods_record, envir = EDI:::EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
	expect_error(
		EDI:::mark_inference_class_migrated(
			"InferenceTemporaryMissingCapabilityMethods",
			public_method_names = character()
		),
		"missing public method"
	)
	expect_identical(
		EDI:::get_inference_hierarchy_migration_record("InferenceTemporaryMissingCapabilityMethods")$migration_status,
		"pending"
	)
})

test_that("migration manifest summary reports pending counts by migration dimensions", {
	EDI:::populate_inference_class_registry()
	manifest = EDI:::inference_hierarchy_migration_manifest_as_list()
	summary = EDI:::inference_hierarchy_migration_summary(manifest)
	pending = Filter(function(record) {
		identical(record$migration_status, "pending") && !isTRUE(record$current_abstract)
	}, manifest)

	expect_named(summary, c("by_current_parent", "by_likelihood_tier", "by_response_type", "by_kk_status"))
	for (table in summary) {
		expect_s3_class(table, "data.frame")
		expect_named(table, c("metric", "value", "n"))
	}
	expect_equal(sum(summary$by_current_parent$n), length(pending))
	expect_equal(sum(summary$by_likelihood_tier$n), length(pending))
	expect_equal(sum(summary$by_kk_status$n), length(pending))
	expect_equal(
		sum(summary$by_response_type$n),
		sum(vapply(pending, function(record) length(record$current_response_types), integer(1L)))
	)
})

test_that("no-likelihood concrete classes are grouped for piecewise migration", {
	EDI:::populate_inference_class_registry()
	manifest = EDI:::inference_hierarchy_migration_manifest_as_list()
	groups = EDI:::inference_no_likelihood_concrete_groups(manifest)
	no_likelihood = names(Filter(function(record) {
		!isTRUE(record$current_abstract) &&
			identical(record$current_likelihood_tier, "none")
	}, manifest))

	expect_named(groups, c("exact", "randomization", "jackknife_asymptotic", "pure_estimator"))
	expect_setequal(unlist(groups, use.names = FALSE), no_likelihood)
	expect_equal(anyDuplicated(unlist(groups, use.names = FALSE)), 0L)
	expect_true(all(c(
		"InferenceIncidExactBinomial",
		"InferenceIncidExactFisher",
		"InferenceIncidenceExactZhang"
	) %in% groups$exact))
})

test_that("exact incidence behavior manifest records current surface and target components", {
	EDI:::populate_inference_class_registry()
	manifest = EDI:::exact_incidence_behavior_manifest()

	expect_named(manifest, c(
		"InferenceIncidExactBinomial",
		"InferenceIncidExactFisher",
		"InferenceIncidenceExactZhang"
	))
	for (name in names(manifest)) {
		record = manifest[[name]]
		expect_identical(record$name, name)
		expect_true(is.character(record$current_public_optional_methods))
		expect_identical(record$intentional_capabilities, "exact_test")
		expect_setequal(record$intentional_public_methods, c(
			"compute_exact_confidence_interval",
			"compute_exact_two_sided_pval_for_treatment_effect"
		))
		expect_true(all(record$legacy_optional_surface %in% record$current_public_optional_methods))
		expect_true("ExactTest" %in% record$target_components)
		expect_true(nzchar(record$notes))
	}
	expect_true("ExactBinomialIncidence" %in% manifest$InferenceIncidExactBinomial$target_components)
	expect_true("ExactFisherIncidence" %in% manifest$InferenceIncidExactFisher$target_components)
	expect_true("ExactZhangIncidence" %in% manifest$InferenceIncidenceExactZhang$target_components)
	expect_true(length(manifest$InferenceIncidenceExactZhang$legacy_optional_surface) > 0L)
})

test_that("CMH is classified as blocked-incidence Wald behavior, not exact behavior", {
	classification = EDI:::get_no_likelihood_special_classification("InferenceIncidCMH")
	expect_type(classification, "list")
	expect_identical(classification$classification, "wald_blocked_incidence")
	expect_false(isTRUE(classification$exact_test_component))
	expect_true(grepl("Wald", classification$notes))
})

test_that("custom randomization has no remaining direct InferenceRand hosts", {
	EDI:::populate_inference_class_registry()
	hosts = EDI:::custom_randomization_direct_host_names()
	target_names = names(EDI:::EDI_CUSTOM_RANDOMIZATION_TARGETS)

	expect_identical(hosts, character())
	expect_identical(target_names, "InferenceCustomRand")
})

test_that("custom randomization hosts are split by migration role", {
	EDI:::populate_inference_class_registry()
	groups = EDI:::custom_randomization_host_groups()

	expect_named(groups, c("extension_bases", "concrete_package_estimators"))
	expect_identical(groups$extension_bases, "InferenceCustomRand")
	expect_identical(groups$concrete_package_estimators, character())
	expect_setequal(
		unlist(groups, use.names = FALSE),
		EDI:::custom_randomization_host_names()
	)
})

test_that("custom randomization behavior manifest records target metadata", {
	EDI:::populate_inference_class_registry()
	manifest = EDI:::custom_randomization_behavior_manifest()

	expect_named(manifest, "InferenceCustomRand")
	record = manifest$InferenceCustomRand
	expect_identical(record$name, "InferenceCustomRand")
	expect_identical(record$host_kind, "extension_base")
	expect_identical(record$current_parent, "Inference")
	expect_identical(record$target_parent, "Inference")
	expect_identical(record$target_components, "RandomizationTest")
	expect_identical(record$class_owned_capabilities, character())
	expect_identical(record$intentional_capabilities, "randomization_test")
	expect_identical(record$migration_status, "migrated")
	expect_setequal(record$migration_evidence, c("method_snapshot", "golden_randomization"))
	expect_identical(EDI:::get_direct_components("InferenceCustomRand"), "RandomizationTest")
	expect_identical(EDI:::get_effective_components("InferenceCustomRand"), "RandomizationTest")
	expect_identical(EDI:::get_effective_capabilities("InferenceCustomRand"), "randomization_test")
	expect_setequal(record$intentional_public_methods, c(
		"approximate_randomization_distribution_beta_hat_T",
		"compute_rand_two_sided_pval"
	))
	expect_equal(intersect(record$legacy_optional_surface, c(
		"compute_rand_confidence_interval",
		"approximate_rand_bootstrap_distribution_beta_hat_T",
		"compute_rand_bootstrap_two_sided_pval",
		"compute_rand_bootstrap_confidence_interval"
	)), character())
	expect_true(nzchar(record$notes))
})

test_that("reviewed custom randomization hosts are marked migrated through the migration gate", {
	EDI:::populate_inference_class_registry()

	expect_silent(EDI:::mark_custom_randomization_classes_migrated("InferenceCustomRand"))
	expect_identical(
		EDI:::get_inference_hierarchy_migration_record("InferenceCustomRand")$migration_status,
		"migrated"
	)
})

test_that("migration order lists leaf concrete classes before concrete parents", {
	manifest = list(
		InferenceTemporaryConcreteParent = list(
			name = "InferenceTemporaryConcreteParent",
			current_parent = "Inference",
			current_abstract = FALSE,
			migration_status = "pending"
		),
		InferenceTemporaryConcreteChild = list(
			name = "InferenceTemporaryConcreteChild",
			current_parent = "InferenceTemporaryConcreteParent",
			current_abstract = FALSE,
			migration_status = "pending"
		),
		InferenceTemporaryConcreteGrandchild = list(
			name = "InferenceTemporaryConcreteGrandchild",
			current_parent = "InferenceTemporaryConcreteChild",
			current_abstract = FALSE,
			migration_status = "pending"
		),
		InferenceTemporaryAbstractChild = list(
			name = "InferenceTemporaryAbstractChild",
			current_parent = "InferenceTemporaryConcreteParent",
			current_abstract = TRUE,
			migration_status = "infrastructure"
		)
	)

	order = EDI:::inference_hierarchy_migration_order(manifest)
	expect_true(match("InferenceTemporaryConcreteGrandchild", order) < match("InferenceTemporaryConcreteChild", order))
	expect_true(match("InferenceTemporaryConcreteChild", order) < match("InferenceTemporaryConcreteParent", order))
	expect_false("InferenceTemporaryAbstractChild" %in% order)
})

test_that("strict shallow hierarchy flag fails when concrete classes are pending", {
	old_value = Sys.getenv("EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY", unset = NA_character_)
	on.exit({
		if (is.na(old_value)) {
			Sys.unsetenv("EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY")
		} else {
			Sys.setenv(EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY = old_value)
		}
	}, add = TRUE)

	manifest = list(
		InferenceTemporaryPending = list(
			name = "InferenceTemporaryPending",
			current_abstract = FALSE,
			migration_status = "pending"
		)
	)

	Sys.unsetenv("EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY")
	expect_false(EDI:::edi_require_shallow_inference_hierarchy())
	expect_silent(EDI:::assert_shallow_inference_hierarchy_complete(manifest))

	Sys.setenv(EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY = "true")
	expect_true(EDI:::edi_require_shallow_inference_hierarchy())
	expect_error(
		EDI:::assert_shallow_inference_hierarchy_complete(manifest),
		"still pending migration"
	)

	manifest$InferenceTemporaryPending$migration_status = "migrated"
	expect_silent(EDI:::assert_shallow_inference_hierarchy_complete(manifest))
})

test_that("strict shallow hierarchy flag gates the current migration manifest", {
	EDI:::populate_inference_class_registry()
	expect_silent(EDI:::assert_shallow_inference_hierarchy_complete())
})

test_that("resolve_inference_components expands dependencies deterministically", {
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)
	EDI:::register_inference_class(
		name = "InferenceTemporaryComponentDependencyHost",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "KKCompound"
	)

	expect_identical(
		EDI:::resolve_inference_components("InferenceTemporaryComponentDependencyHost"),
		c("KKPassThrough", "KKCompound")
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
		direct_components = "KKPassThrough"
	)
	EDI:::register_inference_class(
		name = "InferenceTemporaryChildComponentHost",
		parent = "InferenceTemporaryParentComponentHost",
		metadata = temporary_inference_metadata(likelihood_tier = "full"),
		direct_components = "OffOptimumLikelihoodEval"
	)

	expect_identical(
		EDI:::get_direct_components("InferenceTemporaryChildComponentHost"),
		"OffOptimumLikelihoodEval"
	)
	expect_identical(
		EDI:::resolve_inference_components("InferenceTemporaryChildComponentHost"),
		c("KKPassThrough", "OffOptimumLikelihoodEval")
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
		direct_components = "KKPassThrough"
	)
	EDI:::register_inference_class(
		name = "InferenceTemporaryDependencyChild",
		parent = "InferenceTemporaryDependencyParent",
		metadata = temporary_inference_metadata(),
		direct_components = "KKCompound"
	)

	expect_identical(
		EDI:::resolve_inference_components("InferenceTemporaryDependencyChild"),
		c("KKPassThrough", "KKCompound")
	)
})

test_that("component resolution rejects inherited re-listing and duplicate transitive components", {
	on.exit(EDI:::populate_inference_class_registry(), add = TRUE)
	EDI:::register_inference_class(
		name = "InferenceTemporaryRelistParent",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = "KKPassThrough"
	)
	EDI:::register_inference_class(
		name = "InferenceTemporaryRelistChild",
		parent = "InferenceTemporaryRelistParent",
		metadata = temporary_inference_metadata(),
		direct_components = "KKPassThrough"
	)
	expect_error(
		EDI:::resolve_inference_components("InferenceTemporaryRelistChild"),
		"re-lists inherited component"
	)

	EDI:::register_inference_class(
		name = "InferenceTemporaryDuplicateTransitiveHost",
		parent = "Inference",
		metadata = temporary_inference_metadata(),
		direct_components = c("KKPassThrough", "KKCompound")
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
