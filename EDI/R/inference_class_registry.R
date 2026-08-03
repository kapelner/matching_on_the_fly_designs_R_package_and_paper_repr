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
			direct_components = character()
		)
	}
	validate_inference_class_registry()
	invisible(EDI_INFERENCE_CLASS_REGISTRY)
}

populate_inference_component_registry()
populate_inference_class_registry()
