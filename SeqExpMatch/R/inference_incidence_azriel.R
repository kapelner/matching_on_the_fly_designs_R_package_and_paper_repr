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
			if (!private$design_is_supported_blocking(des_obj)) {
				stop(class(self)[1], " requires a blocking design.")
			}
			private$m = des_obj$get_block_ids()
			if (is.null(private$m)){
				stop(class(self)[1], " requires a blocking design.")
			}
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		m = NULL,

		design_is_supported_blocking = function(des_obj){
			is(des_obj, "DesignSeqOneByOneKK14") ||
				is(des_obj, "FixedDesignBlocking") ||
				is(des_obj, "DesignSeqOneByOneSPBR") ||
				is(des_obj, "DesignSeqOneByOneRandomBlockSize")
		},

		get_standard_error = function(){
			if (!is.null(private$cached_values$azriel_s_beta_hat_T)) {
				return(private$cached_values$azriel_s_beta_hat_T)
			}
			private$cached_values$azriel_s_beta_hat_T = compute_azriel_block_se_cpp(
				private$des_obj_priv_int$y,
				private$m,
				private$des_obj_priv_int$n
			)
			private$cached_values$azriel_s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			NA_real_
		}
	)
)
