#' Azriel Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a block-stratified standard error.
#'
#' @export
InferenceIncidAzriel = R6::R6Class("InferenceIncidAzriel",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description
		#' Initialize Azriel blocked-design incidence inference.
		#' @param des_obj A completed design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceIncidAzriel} object.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			des_obj$assert_even_allocation()
			des_obj$assert_equal_block_sizes()
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
			private$cached_values$is_z = TRUE
		}
	),

	private = list(
		get_standard_error = function(){
			if (!is.null(private$cached_values$azriel_s_beta_hat_T)) {
				return(private$cached_values$azriel_s_beta_hat_T)
			}
			private$cached_values$azriel_s_beta_hat_T = compute_azriel_block_se_cpp(
				private$des_obj_priv_int$y,
				private$des_obj$get_block_ids(),
				private$des_obj_priv_int$n
			)
			private$cached_values$azriel_s_beta_hat_T
		}
	)
)
