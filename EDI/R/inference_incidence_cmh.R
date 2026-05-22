#' CMH Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a block-stratified standard error.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneRandomBlockSize$new(n = 20, response_type = 'incidence', strata_cols = 'x1')
#' for (i in 1:20) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = factor(rep(1:2, 10)[i], levels=1:2)))
#' }
#' seq_des$add_all_subject_responses(rbinom(20, 1, 0.5))
#' inf = InferenceIncidCMH$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceIncidCMH = R6::R6Class("InferenceIncidCMH",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description Computes an approximate confidence interval.
		#' @param alpha Numeric. Significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$get_standard_error()
			super$compute_asymp_confidence_interval(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Numeric. Null treatment effect value (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			private$get_standard_error()
			super$compute_asymp_two_sided_pval(delta)
		},
		#' @description Initialize CMH incidence inference.
		#' @param des_obj A completed design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param se_est_num_vectors For non-block designs, the number of randomization vectors
		#'   drawn from the design to estimate the standard error. Default \code{1000L}.
		#' @param verbose Logical. Whether to print progress messages.
		#' @return A new \code{InferenceIncidCMH} object.
		initialize = function(des_obj, model_formula = NULL, se_est_num_vectors = 1000L, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
				assertCount(se_est_num_vectors, positive = TRUE)
				if (des_obj$is_blocking_design()) {
					des_obj$assert_even_allocation()
					des_obj$assert_equal_block_sizes()
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$se_est_num_vectors = as.integer(se_est_num_vectors)
		}
	),
	private = list(
		se_est_num_vectors = 1000L,
		get_standard_error = function(){
			if (!is.null(private$cached_values$cmh_s_beta_hat_T)) {
				return(private$cached_values$cmh_s_beta_hat_T)
			}
			if (private$des_obj$is_blocking_design()) {
				private$cached_values$cmh_s_beta_hat_T = compute_cmh_block_se_cpp(
					private$des_obj_priv_int$y,
					private$des_obj$get_block_ids(),
					private$des_obj_priv_int$n
				)
			} else {
				precomp = private$des_obj$get_cmh_se_w_mat()
				w_mat = if (!is.null(precomp)) precomp else private$des_obj$draw_ws_according_to_design(private$se_est_num_vectors)
				ytw      = drop(private$y %*% w_mat)
				# E[y·W] = (n/2)·ȳ is known exactly — use it rather than mean(ytw).
				eytw_sq  = ((private$n / 2L) * mean(private$y))^2
				private$cached_values$cmh_s_beta_hat_T = 4 / private$n * sqrt(sum(ytw^2) / (length(ytw) - 1L) - eytw_sq)
			}
			private$cached_values$s_beta_hat_T = private$cached_values$cmh_s_beta_hat_T
			private$cached_values$df = NA_real_
			private$cached_values$cmh_s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			NA_real_
		}
	)
)
