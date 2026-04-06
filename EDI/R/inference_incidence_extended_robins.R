#' Extended Robins Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a block-stratified standard error.
#'
#' @export
InferenceIncidExtendedRobins = R6::R6Class("InferenceIncidExtendedRobins",
	lock_objects = FALSE,
	inherit = InferenceIncidAzriel,

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
