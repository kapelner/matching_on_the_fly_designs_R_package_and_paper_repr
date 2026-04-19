#' Cox Proportional Hazards Regression Inference for Survival Responses
#'
#' Fits a Cox proportional hazards regression for survival responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceSurvivalCoxPHRegr = R6::R6Class("InferenceSurvivalCoxPHRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a Cox PH inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose)
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
		},

		#' @description
		#' Computes the Cox PH estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			super$compute_treatment_estimate(estimate_only = estimate_only)
		}
	),

	private = list(
		include_covariates = NULL,

		generate_mod = function(estimate_only = FALSE){
			surv_obj = survival::Surv(private$y, private$dead)
			X_rhs = if (private$include_covariates) {
				cbind(private$w, private$get_X())
			} else {
				private$w
			}
			tryCatch({
				coxph_mod = suppressWarnings(survival::coxph(surv_obj ~ X_rhs))

				if (estimate_only) {
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = NA_real_
					)
				} else {
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = as.numeric(stats::vcov(coxph_mod))[1]
					)
				}
			}, error = function(e){
				list(
					b = c(NA_real_, NA_real_),
					ssq_b_2 = NA_real_
				)
			})
		}
	)
)
