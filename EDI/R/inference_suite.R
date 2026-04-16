#' Inference Suite
#'
#' Bundles a design object with a set of inference classes and their
#' constructor arguments.  On initialization the suite automatically
#' introspects the package namespace to discover every concrete
#' \code{Inference} subclass that is compatible with the supplied design
#' object, so the applicable list never goes stale as new classes are added.
#'
#' Compatibility is determined by attempting to construct each candidate class
#' with \code{des_obj}.  A class is deemed \emph{not} applicable only when its
#' initializer explicitly rejects the combination via a response-type or
#' design-type mismatch error; any other outcome (successful construction
#' \emph{or} an unrelated error such as a missing optional package) is treated
#' as applicable.
#'
#' @export
InferenceSuite = R6::R6Class("InferenceSuite",
	lock_objects = FALSE,
	public = list(

		#' @field applicable_design_classes Character vector of applicable inference
		#'   class names derived during initialization.
		applicable_design_classes = NULL,

		#' @description
		#' Initialize an \code{InferenceSuite}.
		#' @param des_obj A completed \code{Design} object.
		#' @param inference_params A named list of lists supplying additional
		#'   constructor arguments for specific inference classes.  Each name
		#'   must be the name of a concrete \code{Inference} subclass; the
		#'   corresponding list contains keyword arguments (beyond
		#'   \code{des_obj}) forwarded to that class's \code{initialize}.
		#'   Defaults to an empty list (no extra arguments for any class).
		initialize = function(des_obj, inference_params = list()) {
			if (!is(des_obj, "Design")) {
				stop("InferenceSuite: des_obj must be a Design object.")
			}
			if (!is.list(inference_params)) {
				stop("InferenceSuite: inference_params must be a list.")
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
				if (!(cls_name %in% self$applicable_design_classes)) {
					stop(sprintf(
						"InferenceSuite: '%s' is not applicable for this design/response_type combination.",
						cls_name))
				}

				# All supplied param names must be formals of initialize
				params = inference_params[[cls_name]]
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

			private$des_obj          = des_obj
			private$inference_params = inference_params
		}
	),

	private = list(
		des_obj          = NULL,
		inference_params = NULL,

		# Known non-instantiable base / infrastructure classes
		.base_class_names = c(
			"Inference",
			"InferenceAsymp", "InferenceBoot", "InferenceRand",
			"InferenceRandCI", "InferenceExact",
			"InferenceKKPassThrough", "InferenceKKPassThroughCompound",
			"InferenceMLEorKMforGLMs", "InferenceMLEorKMSummaryTable"
		),

		# Error patterns that signal a genuine incompatibility
		.incompat_pattern =
			"only available for|requires a .* design|only available for uncensored|not recommended for",

		.is_inference_subclass = function(obj) {
			if (!inherits(obj, "R6ClassGenerator")) return(FALSE)
			anc = obj
			while (!is.null(anc)) {
				if (identical(anc$classname, "Inference")) return(TRUE)
				anc = anc$get_inherit()
			}
			FALSE
		},

		.discover_applicable_design_classes = function(des_obj) {
			ns        = getNamespace("EDI")
			all_names = ls(ns)

			# Filter to concrete Inference subclasses
			candidates = Filter(function(nm) {
				if (nm %in% private$.base_class_names)   return(FALSE)
				if (grepl("Abstract", nm, fixed = TRUE)) return(FALSE)
				private$.is_inference_subclass(get(nm, envir = ns))
			}, all_names)

			# Probe compatibility via try-construction (stdout suppressed to silence cat() calls)
			applicable = Filter(function(nm) {
				cls = get(nm, envir = ns)
				result = NULL
				utils::capture.output(
					result <- tryCatch(
						suppressWarnings(suppressMessages(cls$new(des_obj))),
						error = function(e) e
					)
				)
				if (inherits(result, "error")) {
					msg = conditionMessage(result)
					return(!grepl(private$.incompat_pattern, msg, ignore.case = TRUE))
				}
				TRUE
			}, candidates)

			sort(applicable)
		}
	)
)
