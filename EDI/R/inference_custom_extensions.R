#' Internal base for user-defined asymptotic inference extensions
#'
#' \code{InferenceCustomAsymp} is intentionally not exported. Extension packages
#' may retrieve it with \code{getFromNamespace("InferenceCustomAsymp", "EDI")}
#' while this API is experimental.
#'
#' Subclasses implement a public \code{fit(estimate_only = FALSE)} method and
#' return a named list with the custom-fit result contract:
#' \describe{
#'   \item{\code{estimate}}{Required numeric scalar treatment-effect estimate.}
#'   \item{\code{se}}{Optional numeric scalar standard error. Required for Wald
#'     confidence intervals and asymptotic p-values unless \code{estimate_only}
#'     is \code{TRUE}.}
#'   \item{\code{df}}{Optional numeric scalar degrees of freedom. Use
#'     \code{NA_real_} for z inference.}
#'   \item{\code{model}}{Optional fitted model object retained for
#'     \code{get_mod()} and \code{get_summary()}.}
#'   \item{\code{nonestimable_reason}}{Optional character scalar. When supplied
#'     with a non-finite estimate or standard error, EDI records the result as
#'     explicitly non-estimable.}
#' }
#'
#' Subclasses should use public accessors such as \code{get_analysis_data()},
#' \code{get_response()}, \code{get_treatment()}, and \code{get_covariates()}
#' rather than EDI private fields.
#'
#' @keywords internal
InferenceCustomAsymp = R6::R6Class("InferenceCustomAsymp",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description User-defined fit method.
		#' @param estimate_only If TRUE, skip variance calculations.
		#' @return A list with fit results.
		fit = function(estimate_only = FALSE){
			stop("Custom inference subclasses must implement public$fit(estimate_only = FALSE).")
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only If TRUE, skip variance calculations.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			private$run_custom_fit(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute asymptotic confidence interval.
		#' @param alpha Significance level.
		#' @return Confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$run_custom_fit(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute asymptotic p-value.
		#' @param delta Null treatment effect.
		#' @return P-value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$run_custom_fit(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		run_custom_fit = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			result = self$fit(estimate_only = estimate_only)
			private$cache_custom_fit_result(result, estimate_only = estimate_only)
		},
		cache_custom_fit_result = function(result, estimate_only = FALSE){
			if (!is.list(result)) {
				stop("Custom inference fit() must return a named list.", call. = FALSE)
			}
			if (is.null(result$estimate) || length(result$estimate) != 1L) {
				stop("Custom inference fit() result must include numeric scalar 'estimate'.", call. = FALSE)
			}
			estimate = as.numeric(result$estimate)[1L]
			reason = result$nonestimable_reason
			if (!is.null(reason)) reason = as.character(reason)[1L]
			private$cached_values$beta_hat_T = estimate
			private$cached_mod = result$model
			if (!is.finite(estimate)) {
				private$cache_nonestimable_estimate(if (is.null(reason)) "custom_estimate_unavailable" else reason)
				return(invisible(NULL))
			}
			if (isTRUE(estimate_only)) {
				return(invisible(NULL))
			}
			se = if (is.null(result$se)) NA_real_ else as.numeric(result$se)[1L]
			df = if (is.null(result$df)) NA_real_ else as.numeric(result$df)[1L]
			private$cached_values$s_beta_hat_T = se
			private$cached_values$df = df
			if (!is.finite(se) || se <= 0) {
				private$cache_nonestimable_se(if (is.null(reason)) "custom_standard_error_unavailable" else reason)
			}
			invisible(NULL)
		},
		get_standard_error = function(){
			private$run_custom_fit(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			private$run_custom_fit(estimate_only = FALSE)
			private$cached_values$df
		}
	)
)
#' Internal base for user-defined randomization inference extensions
#'
#' This class uses the same \code{fit()} result contract as
#' \code{InferenceCustomAsymp}, but only promises estimate/randomization
#' behavior.
#'
#' @keywords internal
InferenceCustomRand = R6::R6Class("InferenceCustomRand",
	lock_objects = FALSE,
	inherit = InferenceRand,
	public = list(
		#' @description User-defined fit method.
		#' @param estimate_only If TRUE, skip variance calculations.
		#' @return A list with fit results.
		fit = function(estimate_only = FALSE){
			stop("Custom inference subclasses must implement public$fit(estimate_only = FALSE).")
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only If TRUE, skip variance calculations.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(private$cached_values$beta_hat_T)
			result = self$fit(estimate_only = estimate_only)
			if (!is.list(result) || is.null(result$estimate) || length(result$estimate) != 1L) {
				stop("Custom inference fit() result must include numeric scalar 'estimate'.", call. = FALSE)
			}
			private$cached_values$beta_hat_T = as.numeric(result$estimate)[1L]
			private$cached_mod = result$model
			if (!is.finite(private$cached_values$beta_hat_T)) {
				reason = if (is.null(result$nonestimable_reason)) "custom_estimate_unavailable" else as.character(result$nonestimable_reason)[1L]
				private$cache_nonestimable_estimate(reason)
			}
			private$cached_values$beta_hat_T
		}
	)
)
#' Internal base for user-defined bootstrap inference extensions
#'
#' This class uses the same \code{fit()} result contract as
#' \code{InferenceCustomAsymp}, but only promises estimate/bootstrap behavior.
#'
#' @keywords internal
InferenceCustomBoot = R6::R6Class("InferenceCustomBoot",
	lock_objects = FALSE,
	inherit = InferenceNonParamBootstrap,
	public = list(
		#' @description User-defined fit method.
		#' @param estimate_only If TRUE, skip variance calculations.
		#' @return A list with fit results.
		fit = function(estimate_only = FALSE){
			stop("Custom inference subclasses must implement public$fit(estimate_only = FALSE).")
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only If TRUE, skip variance calculations.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(private$cached_values$beta_hat_T)
			result = self$fit(estimate_only = estimate_only)
			if (!is.list(result) || is.null(result$estimate) || length(result$estimate) != 1L) {
				stop("Custom inference fit() result must include numeric scalar 'estimate'.", call. = FALSE)
			}
			private$cached_values$beta_hat_T = as.numeric(result$estimate)[1L]
			private$cached_mod = result$model
			if (!is.finite(private$cached_values$beta_hat_T)) {
				reason = if (is.null(result$nonestimable_reason)) "custom_estimate_unavailable" else as.character(result$nonestimable_reason)[1L]
				private$cache_nonestimable_estimate(reason)
			}
			private$cached_values$beta_hat_T
		}
	)
)
