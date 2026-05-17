#' Simple Mean Difference Inference with Pooled Variance
#'
#' Unadjusted mean-difference inference using the simple treated-minus-control
#' difference with pooled equal-variance t inference.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceAllSimpleMeanDiffPooledVar$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceAllSimpleMeanDiffPooledVar = R6::R6Class("InferenceAllSimpleMeanDiffPooledVar",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description Initialize simple pooled-variance inference.
		#' @param des_obj A completed design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		#' @return A new \code{InferenceAllSimpleMeanDiffPooledVar} object.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = TRUE){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			super$compute_estimate(estimate_only = estimate_only)
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			super$compute_asymp_confidence_interval(alpha = alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			super$compute_asymp_two_sided_pval(delta = delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		}
	),
	private = list(
		get_standard_error = function(){
			if (is.null(private$cached_values$simple_mean_diff_pooled_se)) {
				private$compute_simple_mean_diff_pooled_components()
			}
			private$cached_values$simple_mean_diff_pooled_se
		},
		get_degrees_of_freedom = function(){
			if (is.null(private$cached_values$simple_mean_diff_pooled_df)) {
				private$compute_simple_mean_diff_pooled_components()
			}
			private$cached_values$simple_mean_diff_pooled_df
		},
		compute_simple_mean_diff_pooled_components = function(){
			if (is.null(private$cached_values$beta_hat_T)) {
				self$compute_estimate()
			}
			y_t = private$cached_values$yTs
			y_c = private$cached_values$yCs
			n_t = length(y_t)
			n_c = length(y_c)
			if (n_t <= 1L || n_c <= 1L) {
				private$cached_values$simple_mean_diff_pooled_se = NA_real_
				private$cached_values$simple_mean_diff_pooled_df = NA_real_
				return(invisible(NULL))
			}
			s2_t = stats::var(y_t)
			s2_c = stats::var(y_c)
			df = n_t + n_c - 2L
			s2_pooled = ((n_t - 1L) * s2_t + (n_c - 1L) * s2_c) / df
			var_hat = s2_pooled * (1 / n_t + 1 / n_c)
			private$cached_values$simple_mean_diff_pooled_se =
				if (is.finite(var_hat) && var_hat >= 0) sqrt(var_hat) else NA_real_
			private$cached_values$simple_mean_diff_pooled_df =
				if (is.finite(var_hat) && is.finite(df) && df > 0) {
					as.numeric(df)
				} else {
					NA_real_
				}
			invisible(NULL)
		}
	)
)
