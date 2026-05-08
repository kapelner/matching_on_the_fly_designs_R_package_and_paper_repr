#' Extended Robins Blocked Incidence Inference
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
#' inf = InferenceIncidExtendedRobins$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceIncidExtendedRobins = R6::R6Class("InferenceIncidExtendedRobins",
	lock_objects = FALSE,
	inherit = InferenceIncidCMH,

	private = list(
		get_standard_error = function(){
			if (!is.null(private$cached_values$robins_s_beta_hat_T)) {
				return(private$cached_values$robins_s_beta_hat_T)
			}
			private$cached_values$robins_s_beta_hat_T = compute_extended_robins_block_se_cpp(
				private$des_obj_priv_int$y,
				private$des_obj_priv_int$w,
				private$des_obj$get_block_ids(),
				private$des_obj_priv_int$n
			)
			private$cached_values$robins_s_beta_hat_T
		}
	)
)
