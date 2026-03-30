#' Exact Inference API
#'
#' @name InferenceExact
#' @description Internal method.
#' Thin interface for exact inference methods.
#'
#' @keywords internal
InferenceExact = R6::R6Class("InferenceExact",
	lock_objects = FALSE,
	inherit = Inference,
	public = list(
		#' @description
		#' Computes an exact confidence interval for the requested exact method.
		#' @param alpha            Significance level.
		#' @param pval_epsilon     Reserved for future exact methods that invert p-values.
		#' @param type             Exact inference type. Concrete subclasses provide a default.
		#' @param args_for_type    Structured per-type arguments as a list.
		compute_exact_confidence_interval = function(alpha = 0.05, pval_epsilon = 0.005, type = NULL, args_for_type = NULL){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			exact_type = private$resolve_exact_type(type)
			exact_args = private$normalize_exact_inference_args(exact_type, args_for_type = args_for_type)
			private$compute_exact_confidence_interval_by_type(exact_type, alpha, exact_args)
		},

		#' @description
		#' Computes an exact two-sided p-value for the requested exact method.
		#' @param delta            Null treatment effect on the log-odds-ratio scale.
		#' @param type             Exact inference type. Concrete subclasses provide a default.
		#' @param args_for_type    Structured per-type arguments as a list.
		compute_exact_two_sided_pval_for_treatment_effect = function(delta = 0, type = NULL, args_for_type = NULL){
			assertNumeric(delta, len = 1)
			exact_type = private$resolve_exact_type(type)
			exact_args = private$normalize_exact_inference_args(exact_type, args_for_type = args_for_type)
			private$compute_exact_two_sided_pval_for_treatment_effect_by_type(exact_type, delta, exact_args)
		}
	),

	private = list(
		default_exact_type = NULL,

		resolve_exact_type = function(type){
			if (is.null(type)) type = private$default_exact_type
			assertString(type, min.chars = 1)
			type
		},

		normalize_exact_inference_args = function(type, args_for_type = NULL){
			assertList(args_for_type, null.ok = TRUE)
			utils::modifyList(setNames(list(list()), type), if (is.null(args_for_type)) list() else args_for_type)
		},

		assert_exact_inference_params = function(type, args_for_type){
			assertString(type, min.chars = 1)
			assertList(args_for_type)
			if (!(type %in% names(args_for_type))) stop("args_for_type must contain a list for ", type)
			args = args_for_type[[type]]
			assertList(args)
			invisible(args)
		},

		compute_exact_confidence_interval_by_type = function(type, alpha, args_for_type){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$assert_exact_inference_params(type, args_for_type)
			stop("Exact confidence intervals are not implemented for this exact inference class.")
		},

		compute_exact_two_sided_pval_for_treatment_effect_by_type = function(type, delta, args_for_type){
			assertNumeric(delta, len = 1)
			private$assert_exact_inference_params(type, args_for_type)
			stop("Exact p-values are not implemented for this exact inference class.")
		}
	)
)
