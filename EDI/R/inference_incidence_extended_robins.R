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
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		#' @description Initialize Extended Robins blocked-design incidence inference.
		#' @param des_obj A completed design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceIncidExtendedRobins} object.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
				des_obj$assert_even_allocation()
				des_obj$assert_equal_block_sizes()
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),
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
		},
		get_degrees_of_freedom = function(){
			NA_real_
		}
	)
)
