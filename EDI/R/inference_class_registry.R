#' Inference Class Metadata Registry
#'
#' The registry is the source of truth for structural metadata used by the
#' future shallow inference hierarchy. It is intentionally separated from R6
#' generator environments so metadata can be validated without instantiating
#' inference classes.
#'
#' @keywords internal
#' @noRd
EDI_INFERENCE_CLASS_REGISTRY = new.env(parent = emptyenv())

EDI_INFERENCE_ALLOWED_LIKELIHOOD_TIERS = c("none", "quasi", "partial", "full")

EDI_EXACT_INCIDENCE_CLASS_NAMES = c(
	"InferenceIncidExactBinomial",
	"InferenceIncidExactFisher",
	"InferenceIncidenceExactZhang"
)

EDI_EXACT_INCIDENCE_TARGETS = list(
	InferenceIncidExactBinomial = list(
		intentional_capabilities = "exact_test",
		target_components = c("ExactTest", "ExactBinomialIncidence"),
		notes = "Matched-pair exact binomial inference should retain exact APIs only; inherited bootstrap, randomization, Bayesian-bootstrap, and jackknife APIs are treated as legacy surface."
	),
	InferenceIncidExactFisher = list(
		intentional_capabilities = "exact_test",
		target_components = c("ExactTest", "ExactFisherIncidence"),
		notes = "Fisher/mantelhaen exact incidence inference should retain exact APIs only; inherited bootstrap, randomization, Bayesian-bootstrap, and jackknife APIs are treated as legacy surface."
	),
	InferenceIncidenceExactZhang = list(
		intentional_capabilities = "exact_test",
		target_components = c("ExactTest", "ExactZhangIncidence"),
		notes = "Zhang exact incidence inference should retain exact APIs only; generic resampling-weight methods are legacy surface and are not valid exact Zhang capabilities."
	)
)

EDI_NO_LIKELIHOOD_SPECIAL_CLASSIFICATIONS = list(
	InferenceIncidCMH = list(
		classification = "wald_blocked_incidence",
		exact_test_component = FALSE,
		notes = "CMH is a blocked-incidence Wald/asymptotic class with randomization-based standard-error computation, not an exact-test class."
	)
)

EDI_INFERENCE_ALGORITHM_COMPATIBILITY_BASES = c(
	"InferenceAsymp",
	"InferenceJackknife",
	"InferenceNonParamBootstrap",
	"InferenceBayesianBootstrap",
	"InferenceRand",
	"InferenceRandCI",
	"InferenceRandBootstrap",
	"InferenceRandBootstrapCI",
	"InferenceExact",
	"InferenceAsympLik",
	"InferenceParamBootstrap",
	"InferenceAsympLikStdModCache",
	"InferenceAsympLikStdModCacheNoParamBootstrap",
	"InferenceCountLikelihood",
	"InferenceCountLikelihoodNoParamBootstrap",
	"InferenceCountCompositeLikelihood",
	"InferenceKKPassThroughCompound",
	"InferenceKKPassThroughCompoundNoParamBootstrap",
	"InferenceMLEorKMSummaryTable"
)

EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST = new.env(parent = emptyenv())

EDI_INFERENCE_ABSTRACT_CLASS_NAMES = c(
	"Inference",
	"InferenceAsymp",
	"InferenceJackknife",
	"InferenceNonParamBootstrap",
	"InferenceBayesianBootstrap",
	"InferenceRand",
	"InferenceRandCI",
	"InferenceRandBootstrap",
	"InferenceRandBootstrapCI",
	"InferenceExact",
	"InferenceAsympLik",
	"InferenceParamBootstrap",
	"InferenceKKPassThroughCompound",
	"InferenceKKPassThroughCompoundNoParamBootstrap",
	"InferenceAsympLikStdModCache",
	"InferenceAsympLikStdModCacheNoParamBootstrap",
	"InferenceCountLikelihood",
	"InferenceCountLikelihoodNoParamBootstrap",
	"InferenceCountCompositeLikelihood",
	"InferenceMLEorKMSummaryTable"
)

is_inference_r6_generator = function(obj) {
	if (!inherits(obj, "R6ClassGenerator")) return(FALSE)
	anc = obj
	while (!is.null(anc)) {
		if (identical(anc$classname, "Inference")) return(TRUE)
		anc = anc$get_inherit()
	}
	FALSE
}

infer_inference_response_types = function(name) {
	if (identical(name, "Inference") || grepl("^InferenceAll", name)) {
		return(c("continuous", "incidence", "count", "proportion", "survival", "ordinal"))
	}
	if (grepl("^InferenceContin", name) || grepl("^InferenceBai", name)) return("continuous")
	if (grepl("^InferenceCount", name)) return("count")
	if (grepl("^InferenceIncid", name) || grepl("^InferenceIncidence", name)) return("incidence")
	if (grepl("^InferenceOrdinal", name)) return("ordinal")
	if (grepl("^InferenceProp", name)) return("proportion")
	if (grepl("^InferenceSurvival", name)) return("survival")
	character()
}

infer_inference_likelihood_tier = function(name) {
	if (grepl("GEE|Quasi|Robust|Composite", name)) return("quasi")
	if (grepl("Cox|CondLogit|CondAdjCat|LWA", name)) return("partial")
	if (grepl(
		"OLS|Lin|LogRegr|LogBinomial|Probit|Poisson|NegBin|Hurdle|ZeroInflated|Beta|Weibull|GLMM|CLMM|PropOdds|Cauchit|Cloglog|Stereotype|ContRatio|AdjCat|OrderedProbit|BinomialIdentity|FractionalLogit|CountLikelihood|AsympLik|ParamBootstrap|OneLik|Copula|Frailty",
		name
	)) {
		return("full")
	}
	"none"
}

infer_inference_abstract = function(name) {
	identical(name, "Inference") ||
		name %in% EDI_INFERENCE_ABSTRACT_CLASS_NAMES ||
		grepl("Abstract", name, fixed = TRUE)
}

infer_inference_direct_components = function(name) {
	switch(
		name,
		InferenceRand = "RandomizationTest",
		InferenceRandCI = "RandomizationCI",
		InferenceNonParamBootstrap = "NonparametricBootstrap",
		InferenceRandBootstrap = "RandomizationBootstrap",
		InferenceBayesianBootstrap = "BayesianBootstrap",
		InferenceJackknife = "Jackknife",
		InferenceAsymp = "Wald",
		InferenceExact = "ExactTest",
		InferenceAsympLik = "LikelihoodTests",
		InferenceParamBootstrap = "ParametricLikelihoodBootstrap",
		InferenceAsympLikStdModCache = "StandardModelCache",
		InferenceAsympLikStdModCacheNoParamBootstrap = "StandardModelCache",
		InferenceCountLikelihood = "CountLikelihoodPlumbing",
		InferenceCountLikelihoodNoParamBootstrap = "CountLikelihoodPlumbing",
		InferenceCountCompositeLikelihood = "CountLikelihoodPlumbing",
		InferenceKKPassThroughCompound = "KKCompound",
		InferenceKKPassThroughCompoundNoParamBootstrap = "KKCompound",
		character()
	)
}

always_compatible_inference_metadata = function(design_metadata) {
	TRUE
}

validate_inference_class_metadata = function(metadata) {
	required = c(
		"name", "parent", "abstract", "exported", "response_types",
		"design_families", "compatibility", "likelihood_tier",
		"direct_components", "required_packages"
	)
	missing = setdiff(required, names(metadata))
	if (length(missing) > 0L) {
		stop(sprintf(
			"Inference metadata for %s is missing required field(s): %s",
			metadata$name %||% "<unknown>",
			paste(missing, collapse = ", ")
		), call. = FALSE)
	}
	if ("alias_of" %in% names(metadata)) {
		stop(sprintf(
			"Inference metadata for %s contains forbidden alias metadata.",
			metadata$name
		), call. = FALSE)
	}
	if (!is.character(metadata$name) || length(metadata$name) != 1L || !nzchar(metadata$name)) {
		stop("Inference metadata field `name` must be a non-empty character scalar.", call. = FALSE)
	}
	if (!(is.null(metadata$parent) || (is.character(metadata$parent) && length(metadata$parent) == 1L && nzchar(metadata$parent)))) {
		stop(sprintf("Inference metadata for %s has invalid `parent`.", metadata$name), call. = FALSE)
	}
	if (!is.logical(metadata$abstract) || length(metadata$abstract) != 1L || is.na(metadata$abstract)) {
		stop(sprintf("Inference metadata for %s has invalid `abstract`.", metadata$name), call. = FALSE)
	}
	if (!is.logical(metadata$exported) || length(metadata$exported) != 1L || is.na(metadata$exported)) {
		stop(sprintf("Inference metadata for %s has invalid `exported`.", metadata$name), call. = FALSE)
	}
	if (!is.character(metadata$response_types)) {
		stop(sprintf("Inference metadata for %s has invalid `response_types`.", metadata$name), call. = FALSE)
	}
	if (!is.character(metadata$design_families)) {
		stop(sprintf("Inference metadata for %s has invalid `design_families`.", metadata$name), call. = FALSE)
	}
	if (!is.function(metadata$compatibility)) {
		stop(sprintf("Inference metadata for %s has invalid `compatibility`.", metadata$name), call. = FALSE)
	}
	if (!is.character(metadata$likelihood_tier) ||
			length(metadata$likelihood_tier) != 1L ||
			!(metadata$likelihood_tier %in% EDI_INFERENCE_ALLOWED_LIKELIHOOD_TIERS)) {
		stop(sprintf(
			"Inference metadata for %s has invalid `likelihood_tier`: %s",
			metadata$name,
			paste(metadata$likelihood_tier, collapse = ", ")
		), call. = FALSE)
	}
	if (!is.character(metadata$direct_components)) {
		stop(sprintf("Inference metadata for %s has invalid `direct_components`.", metadata$name), call. = FALSE)
	}
	if (!is.character(metadata$required_packages)) {
		stop(sprintf("Inference metadata for %s has invalid `required_packages`.", metadata$name), call. = FALSE)
	}
	invisible(TRUE)
}

register_inference_class = function(name, parent = NULL, metadata = list(), direct_components = character()) {
	record = utils::modifyList(
		list(
			name = name,
			parent = parent,
			abstract = FALSE,
			exported = FALSE,
			response_types = character(),
			design_families = "all",
			compatibility = always_compatible_inference_metadata,
			likelihood_tier = "none",
			direct_components = direct_components,
			required_packages = character(),
			capabilities = character()
		),
		metadata
	)
	record$name = name
	record["parent"] = list(parent)
	record$direct_components = direct_components
	validate_inference_class_metadata(record)
	if (exists(name, envir = EDI_INFERENCE_CLASS_REGISTRY, inherits = FALSE)) {
		stop(sprintf("Inference class metadata already registered for %s.", name), call. = FALSE)
	}
	assign(name, record, envir = EDI_INFERENCE_CLASS_REGISTRY)
	invisible(record)
}

inference_class_registry_as_list = function() {
	mget(ls(EDI_INFERENCE_CLASS_REGISTRY), envir = EDI_INFERENCE_CLASS_REGISTRY, inherits = FALSE)
}

validate_inference_class_registry = function(registry = inference_class_registry_as_list()) {
	if (!is.list(registry)) {
		stop("Inference class registry must be a list of metadata records.", call. = FALSE)
	}
	for (name in names(registry)) {
		validate_inference_class_metadata(registry[[name]])
		if (!identical(registry[[name]]$name, name)) {
			stop(sprintf("Inference metadata name mismatch for registry key %s.", name), call. = FALSE)
		}
		parent = registry[[name]]$parent
		if (!is.null(parent) && !(parent %in% names(registry))) {
			stop(sprintf("Inference metadata for %s has unregistered parent %s.", name, parent), call. = FALSE)
		}
	}
	invisible(TRUE)
}

get_inference_class_metadata = function(name) {
	if (!exists(name, envir = EDI_INFERENCE_CLASS_REGISTRY, inherits = FALSE)) {
		stop(sprintf("No inference class metadata registered for %s.", name), call. = FALSE)
	}
	get(name, envir = EDI_INFERENCE_CLASS_REGISTRY, inherits = FALSE)
}

get_direct_components = function(name) {
	get_inference_class_metadata(name)$direct_components
}

resolve_inference_components = function(name) {
	metadata = get_inference_class_metadata(name)
	parent_components = if (is.null(metadata$parent)) {
		character()
	} else {
		resolve_inference_components(metadata$parent)
	}
	duplicate_components = intersect(parent_components, metadata$direct_components)
	if (length(duplicate_components) > 0L) {
		stop(sprintf(
			"%s re-lists inherited component(s): %s",
			name,
			paste(duplicate_components, collapse = ", ")
		), call. = FALSE)
	}
	direct_components = tryCatch(
		resolve_component_dependencies(
			metadata$direct_components,
			satisfied_components = parent_components
		),
		error = function(e) {
			stop(sprintf("%s: %s", name, conditionMessage(e)), call. = FALSE)
		}
	)
	duplicate_transitive = intersect(parent_components, direct_components)
	if (length(duplicate_transitive) > 0L) {
		stop(sprintf(
			"%s duplicates inherited transitive component(s): %s",
			name,
			paste(duplicate_transitive, collapse = ", ")
		), call. = FALSE)
	}
	c(parent_components, direct_components)
}

get_effective_components = function(name) {
	resolve_inference_components(name)
}

get_effective_capabilities = function(name) {
	metadata = get_inference_class_metadata(name)
	component_capabilities = as.character(unlist(lapply(get_effective_components(name), function(component_name) {
		get_inference_component(component_name)$provides_capabilities
	}), use.names = FALSE))
	unique(c(component_capabilities, metadata$capabilities %||% character()))
}

inference_class_ancestor_names = function(name, registry = inference_class_registry_as_list()) {
	ancestors = character()
	current = registry[[name]]
	while (!is.null(current) && !is.null(current$parent)) {
		ancestors = c(ancestors, current$parent)
		current = registry[[current$parent]]
	}
	ancestors
}

target_inference_components = function(name) {
	if (identical(name, "Inference")) return(character())
	metadata = get_inference_class_metadata(name)
	if (isTRUE(metadata$abstract)) return(get_effective_components(name))
	get_effective_components(name)
}

target_inference_parent = function(name) {
	metadata = get_inference_class_metadata(name)
	if (identical(name, "Inference")) return(NULL)
	if (isTRUE(metadata$abstract)) return(metadata$parent)
	"Inference"
}

build_inference_hierarchy_migration_record = function(name) {
	metadata = get_inference_class_metadata(name)
	ancestors = inference_class_ancestor_names(name)
	algorithmic_ancestors = intersect(ancestors, EDI_INFERENCE_ALGORITHM_COMPATIBILITY_BASES)
	target_components = target_inference_components(name)
	target_capabilities = unique(c(
		as.character(unlist(lapply(target_components, function(component_name) {
			get_inference_component(component_name)$provides_capabilities
		}), use.names = FALSE)),
		metadata$capabilities %||% character()
	))
	migration_status = if (identical(name, "Inference")) {
		"root"
	} else if (isTRUE(metadata$abstract)) {
		"infrastructure"
	} else if (length(algorithmic_ancestors) == 0L && identical(metadata$parent, "Inference")) {
		"migrated"
	} else {
		"pending"
	}
	list(
		name = name,
		current_parent = metadata$parent,
		current_ancestors = ancestors,
		current_abstract = metadata$abstract,
		current_exported = metadata$exported,
		current_response_types = metadata$response_types,
		current_design_families = metadata$design_families,
		current_likelihood_tier = metadata$likelihood_tier,
		algorithmic_compatibility_ancestors = algorithmic_ancestors,
		target_parent = target_inference_parent(name),
		target_components = target_components,
		target_capabilities = target_capabilities,
		target_likelihood_tier = metadata$likelihood_tier,
		target_response_types = metadata$response_types,
		target_design_families = metadata$design_families,
		migration_status = migration_status
	)
}

validate_inference_hierarchy_migration_record = function(record) {
	required = c(
		"name", "current_parent", "current_ancestors", "current_abstract",
		"current_exported", "current_response_types", "current_design_families",
		"current_likelihood_tier", "algorithmic_compatibility_ancestors",
		"target_parent", "target_components", "target_capabilities",
		"target_likelihood_tier", "target_response_types", "target_design_families",
		"migration_status"
	)
	missing = setdiff(required, names(record))
	if (length(missing) > 0L) {
		stop(sprintf(
			"Inference hierarchy migration record for %s is missing field(s): %s",
			record$name %||% "<unknown>",
			paste(missing, collapse = ", ")
		), call. = FALSE)
	}
	if (!is.character(record$name) || length(record$name) != 1L || !nzchar(record$name)) {
		stop("Migration record `name` must be a non-empty character scalar.", call. = FALSE)
	}
	if (!(record$migration_status %in% c("root", "infrastructure", "pending", "migrated"))) {
		stop(sprintf("Invalid migration status for %s.", record$name), call. = FALSE)
	}
	if (!(record$target_likelihood_tier %in% EDI_INFERENCE_ALLOWED_LIKELIHOOD_TIERS)) {
		stop(sprintf("Invalid target likelihood tier for %s.", record$name), call. = FALSE)
	}
	unknown_components = setdiff(record$target_components, ls(EDI_INFERENCE_COMPONENTS))
	if (length(unknown_components) > 0L) {
		stop(sprintf(
			"Migration record for %s has unknown target component(s): %s",
			record$name,
			paste(unknown_components, collapse = ", ")
		), call. = FALSE)
	}
	if (any(grepl("^InferenceMixin", record$target_components))) {
		stop(sprintf("Migration record for %s contains legacy mixin component names.", record$name), call. = FALSE)
	}
	invisible(TRUE)
}

clear_inference_hierarchy_migration_manifest = function() {
	rm(list = ls(EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST), envir = EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
	invisible(TRUE)
}

populate_inference_hierarchy_migration_manifest = function() {
	clear_inference_hierarchy_migration_manifest()
	for (name in ls(EDI_INFERENCE_CLASS_REGISTRY)) {
		record = build_inference_hierarchy_migration_record(name)
		validate_inference_hierarchy_migration_record(record)
		assign(name, record, envir = EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
	}
	invisible(EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
}

inference_hierarchy_migration_manifest_as_list = function() {
	mget(ls(EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST), envir = EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST, inherits = FALSE)
}

get_inference_hierarchy_migration_record = function(name) {
	if (!exists(name, envir = EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST, inherits = FALSE)) {
		stop(sprintf("No inference hierarchy migration record registered for %s.", name), call. = FALSE)
	}
	get(name, envir = EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST, inherits = FALSE)
}

inference_public_method_names = function(name, generator = NULL) {
	if (is.null(generator)) {
		ns = environment(populate_inference_class_registry)
		if (!exists(name, envir = ns, inherits = FALSE)) {
			stop(sprintf("No R6 generator found for inference class %s.", name), call. = FALSE)
		}
		generator = get(name, envir = ns, inherits = FALSE)
	}
	if (!inherits(generator, "R6ClassGenerator")) {
		stop(sprintf("Object for %s is not an R6 generator.", name), call. = FALSE)
	}
	public_names = character()
	current = generator
	while (!is.null(current)) {
		public_names = c(public_names, names(current$public_methods) %||% character())
		current = current$get_inherit()
	}
	unique(public_names)
}

public_methods_required_for_capabilities = function(capabilities) {
	unique(as.character(unlist(
		public_methods_for_capability[intersect(names(public_methods_for_capability), capabilities)],
		use.names = FALSE
	)))
}

inference_optional_method_names_for_capabilities = function(capabilities) {
	sort(unique(public_methods_required_for_capabilities(capabilities)))
}

validate_inference_public_optional_method_presence = function(name, capabilities = get_effective_capabilities(name), public_method_names = NULL) {
	required_public_methods = public_methods_required_for_capabilities(capabilities)
	if (length(required_public_methods) == 0L) return(invisible(TRUE))
	if (is.null(public_method_names)) {
		public_method_names = inference_public_method_names(name)
	}
	missing_public_methods = setdiff(required_public_methods, public_method_names)
	if (length(missing_public_methods) > 0L) {
		stop(sprintf(
			"%s is missing public method(s) required by effective capabilities: %s",
			name,
			paste(missing_public_methods, collapse = ", ")
		), call. = FALSE)
	}
	invisible(TRUE)
}

mark_inference_class_migrated = function(name, public_method_names = NULL) {
	record = get_inference_hierarchy_migration_record(name)
	metadata = get_inference_class_metadata(name)
	if (identical(record$migration_status, "root") || identical(record$migration_status, "infrastructure") || isTRUE(metadata$abstract)) {
		stop(sprintf("%s is not a concrete inference class and cannot be marked migrated.", name), call. = FALSE)
	}
	if (!identical(metadata$parent, record$target_parent)) {
		stop(sprintf(
			"%s cannot be marked migrated: current parent `%s` does not match target parent `%s`.",
			name,
			metadata$parent %||% "<none>",
			record$target_parent %||% "<none>"
		), call. = FALSE)
	}
	current_ancestors = inference_class_ancestor_names(name)
	algorithmic_ancestors = intersect(current_ancestors, EDI_INFERENCE_ALGORITHM_COMPATIBILITY_BASES)
	if (length(algorithmic_ancestors) > 0L) {
		stop(sprintf(
			"%s cannot be marked migrated: it still descends from algorithmic compatibility base(s): %s",
			name,
			paste(algorithmic_ancestors, collapse = ", ")
		), call. = FALSE)
	}
	effective_components = get_effective_components(name)
	if (!identical(effective_components, record$target_components)) {
		stop(sprintf(
			"%s cannot be marked migrated: effective components do not match target components.",
			name
		), call. = FALSE)
	}
	effective_capabilities = get_effective_capabilities(name)
	if (!identical(effective_capabilities, record$target_capabilities)) {
		stop(sprintf(
			"%s cannot be marked migrated: effective capabilities do not match target capabilities.",
			name
		), call. = FALSE)
	}
	validate_inference_public_optional_method_presence(
		name = name,
		capabilities = effective_capabilities,
		public_method_names = public_method_names
	)
	record$migration_status = "migrated"
	validate_inference_hierarchy_migration_record(record)
	assign(name, record, envir = EDI_INFERENCE_HIERARCHY_MIGRATION_MANIFEST)
	invisible(record)
}

is_kk_inference_migration_record = function(record) {
	grepl("KK|IVWC", record$name) ||
		any(grepl("KK|IVWC", record$current_ancestors)) ||
		any(grepl("^KK", record$target_components))
}

count_inference_migration_records_by = function(records, values_fn, field_name) {
	values = as.character(unlist(lapply(records, values_fn), use.names = FALSE))
	if (length(values) == 0L) {
		return(data.frame(
			metric = character(),
			value = character(),
			n = integer(),
			row.names = NULL,
			stringsAsFactors = FALSE
		))
	}
	counts = sort(table(values), decreasing = TRUE)
	data.frame(
		metric = field_name,
		value = names(counts),
		n = as.integer(counts),
		row.names = NULL,
		stringsAsFactors = FALSE
	)
}

inference_hierarchy_migration_summary = function(manifest = inference_hierarchy_migration_manifest_as_list()) {
	pending = Filter(function(record) {
		identical(record$migration_status, "pending") && !isTRUE(record$current_abstract)
	}, manifest)
	list(
		by_current_parent = count_inference_migration_records_by(
			pending,
			function(record) record$current_parent %||% "<none>",
			"current_parent"
		),
		by_likelihood_tier = count_inference_migration_records_by(
			pending,
			function(record) record$current_likelihood_tier,
			"likelihood_tier"
		),
		by_response_type = count_inference_migration_records_by(
			pending,
			function(record) record$current_response_types,
			"response_type"
		),
		by_kk_status = count_inference_migration_records_by(
			pending,
			function(record) if (is_kk_inference_migration_record(record)) "KK_or_IVWC" else "non_KK",
			"kk_status"
		)
	)
}

inference_hierarchy_migration_order = function(
	manifest = inference_hierarchy_migration_manifest_as_list(),
	statuses = c("pending", "migrated")
) {
	records = Filter(function(record) {
		!isTRUE(record$current_abstract) && record$migration_status %in% statuses
	}, manifest)
	if (length(records) == 0L) return(character())
	class_names = names(records)
	parent_map = stats::setNames(vapply(records, function(record) record$current_parent %||% NA_character_, character(1L)), class_names)
	children_for = function(parent_name) class_names[!is.na(parent_map) & parent_map == parent_name]
	depth_cache = new.env(parent = emptyenv())
	descendant_depth = function(name) {
		if (exists(name, envir = depth_cache, inherits = FALSE)) {
			return(get(name, envir = depth_cache, inherits = FALSE))
		}
		children = children_for(name)
		depth = if (length(children) == 0L) 0L else 1L + max(vapply(children, descendant_depth, integer(1L)))
		assign(name, depth, envir = depth_cache)
		depth
	}
	depths = vapply(class_names, descendant_depth, integer(1L))
	class_names[order(depths, class_names)]
}

edi_require_shallow_inference_hierarchy = function() {
	tolower(Sys.getenv("EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY", unset = "false")) %in% c("1", "true", "yes")
}

assert_shallow_inference_hierarchy_complete = function(manifest = inference_hierarchy_migration_manifest_as_list(), require = edi_require_shallow_inference_hierarchy()) {
	if (!isTRUE(require)) return(invisible(TRUE))
	pending = Filter(function(record) {
		identical(record$migration_status, "pending") && !isTRUE(record$current_abstract)
	}, manifest)
	if (length(pending) > 0L) {
		stop(sprintf(
			"EDI_REQUIRE_SHALLOW_INFERENCE_HIERARCHY is enabled, but %d concrete inference class(es) are still pending migration: %s",
			length(pending),
			paste(utils::head(sort(names(pending)), 10L), collapse = ", ")
		), call. = FALSE)
	}
	invisible(TRUE)
}

inference_no_likelihood_group = function(record) {
	capabilities = record$target_capabilities %||% character()
	ancestors = record$current_ancestors %||% character()
	if ("exact_test" %in% capabilities || "InferenceExact" %in% ancestors || grepl("Exact", record$name)) {
		return("exact")
	}
	if ("randomization_test" %in% capabilities || "InferenceRand" %in% ancestors) {
		return("randomization")
	}
	if (any(c("jackknife", "wald") %in% capabilities) ||
			any(c("InferenceJackknife", "InferenceAsymp") %in% ancestors)) {
		return("jackknife_asymptotic")
	}
	"pure_estimator"
}

inference_no_likelihood_concrete_groups = function(manifest = inference_hierarchy_migration_manifest_as_list()) {
	no_likelihood = Filter(function(record) {
		!isTRUE(record$current_abstract) &&
			identical(record$current_likelihood_tier, "none")
	}, manifest)
	groups = list(
		exact = character(),
		randomization = character(),
		jackknife_asymptotic = character(),
		pure_estimator = character()
	)
	for (record in no_likelihood) {
		group = inference_no_likelihood_group(record)
		groups[[group]] = c(groups[[group]], record$name)
	}
	lapply(groups, sort)
}

get_no_likelihood_special_classification = function(name) {
	EDI_NO_LIKELIHOOD_SPECIAL_CLASSIFICATIONS[[name]] %||% NULL
}

build_exact_incidence_behavior_record = function(name) {
	if (!(name %in% EDI_EXACT_INCIDENCE_CLASS_NAMES)) {
		stop(sprintf("%s is not a registered exact incidence class.", name), call. = FALSE)
	}
	target = EDI_EXACT_INCIDENCE_TARGETS[[name]]
	current_public_methods = inference_public_method_names(name)
	current_optional_methods = sort(unique(intersect(
		current_public_methods,
		inference_optional_method_names_for_capabilities(c(
			"exact_test",
			"nonparametric_bootstrap",
			"randomization_test",
			"randomization_ci",
			"randomization_bootstrap",
			"bayesian_bootstrap",
			"jackknife"
		))
	)))
	intentional_methods = inference_optional_method_names_for_capabilities(target$intentional_capabilities)
	list(
		name = name,
		current_public_optional_methods = current_optional_methods,
		intentional_capabilities = target$intentional_capabilities,
		intentional_public_methods = intentional_methods,
		legacy_optional_surface = setdiff(current_optional_methods, intentional_methods),
		target_components = target$target_components,
		notes = target$notes
	)
}

exact_incidence_behavior_manifest = function() {
	stats::setNames(lapply(EDI_EXACT_INCIDENCE_CLASS_NAMES, build_exact_incidence_behavior_record), EDI_EXACT_INCIDENCE_CLASS_NAMES)
}

clear_inference_class_registry = function() {
	rm(list = ls(EDI_INFERENCE_CLASS_REGISTRY), envir = EDI_INFERENCE_CLASS_REGISTRY)
	invisible(TRUE)
}

populate_inference_class_registry = function(ns = environment(populate_inference_class_registry)) {
	clear_inference_class_registry()
	exports = tryCatch(getNamespaceExports("EDI"), error = function(e) character())
	all_names = sort(ls(ns))
	for (name in all_names) {
		obj = get(name, envir = ns)
		if (!is_inference_r6_generator(obj)) next
		if (!identical(obj$classname, name)) next
		parent_gen = obj$get_inherit()
		parent_name = if (is.null(parent_gen)) NULL else parent_gen$classname
		register_inference_class(
			name = name,
			parent = parent_name,
			metadata = list(
				abstract = infer_inference_abstract(name),
				exported = name %in% exports,
				response_types = infer_inference_response_types(name),
				design_families = "all",
				compatibility = always_compatible_inference_metadata,
				likelihood_tier = infer_inference_likelihood_tier(name),
				required_packages = character(),
				capabilities = character()
			),
			direct_components = infer_inference_direct_components(name)
		)
	}
	validate_inference_class_registry()
	populate_inference_hierarchy_migration_manifest()
	invisible(EDI_INFERENCE_CLASS_REGISTRY)
}

populate_inference_component_registry()
populate_inference_class_registry()
