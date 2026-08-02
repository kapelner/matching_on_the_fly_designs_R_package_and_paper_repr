#' Inference Suite
#'
#' Bundles a design object with a set of inference classes and their
#' constructor arguments.  On initialization the suite automatically
#' introspects the package namespace to discover every concrete
#' \code{Inference} subclass that is compatible with the supplied design
#' object, so the applicable list never goes stale as new classes are added.
#'
#' Discovery uses package metadata rather than constructor probing: only
#' exported, non-abstract \code{Inference} generators are considered, and
#' compatibility is evaluated from the completed design's response and design
#' metadata. Optional-package failures and constructor side effects therefore
#' cannot affect the discovered class list.
#'
#' @export
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "continuous")
#' for (i in 1:20) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(20))
#'
#' suite = InferenceSuite$new(seq_des)
#' suite$applicable_design_classes
#' }
InferenceSuite = R6::R6Class("InferenceSuite",
	lock_objects = FALSE,
	public = list(
		#' @field applicable_design_classes Character vector of applicable inference
		#'   class names derived during initialization.
		applicable_design_classes = NULL,
		#' @description Initialize an \code{\link[EDI:InferenceSuite]{InferenceSuite}}
		#'   wrapper that coordinates multiple inference objects for the same
		#'   completed design and exposes their estimates, p-values, and confidence
		#'   intervals through one object.
		#' @param des_obj A completed \code{Design} object.
		#' @param inference_params A named list of lists supplying additional
		#'   constructor arguments for specific inference classes.  Each name
		#'   must be the name of a concrete \code{Inference} subclass; the
		#'   corresponding list contains keyword arguments (beyond
		#'   \code{des_obj}) forwarded to that class's \code{initialize}.
		#'   Defaults to an empty list (no extra arguments for any class).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		initialize = function(des_obj, model_formula = NULL, inference_params = list()) {
			if (should_run_asserts()) {
				if (!is(des_obj, "Design")) {
					stop("InferenceSuite: des_obj must be a Design object.")
				}
				if (!is.list(inference_params)) {
					stop("InferenceSuite: inference_params must be a list.")
				}
			}
			if (length(inference_params) > 0L &&
					(is.null(names(inference_params)) ||
					 any(nchar(names(inference_params)) == 0L))) {
				stop("InferenceSuite: inference_params must be a fully named list.")
			}
			# ── 1. Discover applicable classes ────────────────────────────────
			self$applicable_design_classes = private$.discover_applicable_design_classes(des_obj)
			# ── 2. Validate inference_params ──────────────────────────────────
			for (cls_name in names(inference_params)) {
				# Must be applicable for this design
				if (should_run_asserts()) {
					if (!(cls_name %in% self$applicable_design_classes)) {
						stop(sprintf(
							"InferenceSuite: '%s' is not applicable for this design/response_type combination.",
							cls_name))
					}
				}
				# All supplied param names must be formals of initialize
				params = inference_params[[cls_name]]
				if (should_run_asserts()) {
					if (!is.list(params)) {
						stop(sprintf(
							"InferenceSuite: params for '%s' must be a list.", cls_name))
					}
					if (length(params) > 0L) {
						cls        = get(cls_name, envir = getNamespace("EDI"))
						init_fn    = cls$public_methods$initialize
						valid_args = setdiff(names(formals(init_fn)), c("des_obj", "..."))
						unknown    = setdiff(names(params), valid_args)
						if (length(unknown) > 0L) {
							stop(sprintf(
								"InferenceSuite: unknown argument(s) for '%s': %s\n  Valid: %s",
								cls_name,
								paste(unknown,    collapse = ", "),
								paste(valid_args, collapse = ", ")))
						}
					}
				}
			}
			private$des_obj          = des_obj
			private$inference_params = inference_params
		}
	),
	private = list(
		des_obj          = NULL,
		inference_params = NULL,
		# Known abstract / infrastructure classes. Discovery also treats any
		# generator whose name contains "Abstract" as abstract.
		.abstract_class_names = c(
			"Inference",
			"InferenceAsymp", "InferenceJackknife", "InferenceNonParamBootstrap", "InferenceBayesianBootstrap", "InferenceRand",
			"InferenceRandCI", "InferenceRandBootstrap", "InferenceRandBootstrapCI", "InferenceExact",
			"InferenceAsympLik", "InferenceParamBootstrap",
			"InferenceKKPassThroughCompound", "InferenceKKPassThroughCompoundNoParamBootstrap",
			"InferenceAsympLikStdModCache", "InferenceAsympLikStdModCacheNoParamBootstrap",
			"InferenceCountLikelihood", "InferenceCountLikelihoodNoParamBootstrap",
			"InferenceCountCompositeLikelihood",
			"InferenceMLEorKMSummaryTable"
		),
		.is_inference_subclass = function(obj) {
			if (!inherits(obj, "R6ClassGenerator")) return(FALSE)
			anc = obj
			while (!is.null(anc)) {
				if (identical(anc$classname, "Inference")) return(TRUE)
				anc = anc$get_inherit()
			}
			FALSE
		},
		.is_exported_inference_class = function(nm) {
			nm %in% getNamespaceExports("EDI")
		},
		.is_abstract_inference_class = function(nm) {
			nm %in% private$.abstract_class_names || grepl("Abstract", nm, fixed = TRUE)
		},
		.design_metadata = function(des_obj) {
			list(
				response_type = des_obj$get_response_type(),
				is_kk = isTRUE(des_obj$is_a_kk_matching_capable()),
				is_blocking = inherits(des_obj, "DesignBlocking") || inherits(des_obj, "DesignFixedBlocking"),
				any_censoring = isTRUE(des_obj$any_censoring())
			)
		},
		.class_compatibility_metadata = function(nm) {
			if (grepl("^InferenceContin", nm) || grepl("^InferenceBai", nm)) {
				response_types = "continuous"
			} else if (grepl("^InferenceCount", nm)) {
				response_types = "count"
			} else if (grepl("^InferenceIncid", nm) || grepl("^InferenceIncidence", nm)) {
				response_types = "incidence"
			} else if (grepl("^InferenceOrdinal", nm)) {
				response_types = "ordinal"
			} else if (grepl("^InferenceProp", nm)) {
				response_types = "proportion"
			} else if (grepl("^InferenceSurvival", nm)) {
				response_types = "survival"
			} else if (grepl("^InferenceAll", nm)) {
				response_types = c("continuous", "incidence", "count", "proportion", "survival", "ordinal")
			} else {
				response_types = character()
			}
			list(
				abstract = private$.is_abstract_inference_class(nm),
				exported = private$.is_exported_inference_class(nm),
				response_types = response_types,
				requires_kk = grepl("KK", nm, fixed = TRUE),
				requires_blocking = nm %in% c("InferenceIncidCMH", "InferenceIncidExtendedRobins"),
				requires_uncensored = nm %in% c(
					"InferenceSurvivalKMDiff",
					"InferenceSurvivalLogRank",
					"InferenceSurvivalRestrictedMeanDiff",
					"InferenceSurvivalGehanWilcox"
				)
			)
			},
			.is_compatible_with_design_metadata = function(nm, design_meta) {
				class_meta = private$.class_compatibility_metadata(nm)
				if (isTRUE(class_meta$abstract) || !isTRUE(class_meta$exported)) return(FALSE)
				if (length(class_meta$response_types) == 0L) return(FALSE)
				if (!(design_meta$response_type %in% class_meta$response_types)) {
					return(FALSE)
				}
				if (isTRUE(class_meta$requires_kk) && !isTRUE(design_meta$is_kk)) return(FALSE)
				if (isTRUE(class_meta$requires_blocking) && !isTRUE(design_meta$is_blocking)) return(FALSE)
				if (isTRUE(class_meta$requires_uncensored) && isTRUE(design_meta$any_censoring)) return(FALSE)
				TRUE
		},
		.discover_applicable_design_classes = function(des_obj) {
			ns        = getNamespace("EDI")
			all_names = ls(ns)
			design_meta = private$.design_metadata(des_obj)
			# Filter to exported, concrete Inference subclasses. Aliases are removed
			# by requiring the namespace binding name to match the generator classname.
			candidates = Filter(function(nm) {
				obj = get(nm, envir = ns)
				private$.is_inference_subclass(obj) &&
					identical(obj$classname, nm) &&
					private$.is_exported_inference_class(nm) &&
					!private$.is_abstract_inference_class(nm)
			}, all_names)
			applicable = Filter(function(nm) {
				private$.is_compatible_with_design_metadata(nm, design_meta)
			}, candidates)
			sort(applicable)
		}
	)
)
