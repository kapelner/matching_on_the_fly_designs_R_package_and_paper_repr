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
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the optimized Rcpp
		#'   implementation. If \code{FALSE}, use \pkg{survival::coxph}.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, use_rcpp = TRUE, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFlag(include_covariates, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			super$initialize(des_obj, verbose)
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
			private$use_rcpp = use_rcpp
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
		use_rcpp = TRUE,

		generate_mod = function(estimate_only = FALSE){
			X_fit = if (private$include_covariates) {
				cbind(treatment = private$w, private$get_X())
			} else {
				matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			}

			if (private$use_rcpp) {
				fit = tryCatch(
					fast_coxph_regression(X_fit, private$y, private$dead, use_rcpp = TRUE, estimate_only = estimate_only),
					error = function(e) NULL
				)
				if (is.null(fit)) {
					return(list(b = rep(NA_real_, ncol(X_fit)), vcov = matrix(NA_real_, ncol(X_fit), ncol(X_fit))))
				}
				# fast_coxph_regression returns coefficients with intercept=FALSE for Cox.
				# InferenceMLEorKMforGLMs expects intercept in first position if it exists.
				# But Cox has no intercept. So we prefix 0.
				return(list(
					b = c(0, fit$b),
					ssq_b_2 = if (estimate_only) NA_real_ else fit$vcov[1, 1],
					vcov = if (estimate_only) NULL else {
						v = matrix(0, ncol(X_fit) + 1, ncol(X_fit) + 1)
						v[2:(ncol(X_fit) + 1), 2:(ncol(X_fit) + 1)] = fit$vcov
						v
					}
				))
			}

			surv_obj = survival::Surv(private$y, private$dead)
			tryCatch({
				coxph_mod = suppressWarnings(survival::coxph(surv_obj ~ X_fit))

				if (estimate_only) {
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = NA_real_,
						vcov = NULL
					)
				} else {
					vcov_mat = stats::vcov(coxph_mod)
					v = matrix(0, ncol(X_fit) + 1, ncol(X_fit) + 1)
					v[2:(ncol(X_fit) + 1), 2:(ncol(X_fit) + 1)] = vcov_mat
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = as.numeric(vcov_mat[1, 1]),
						vcov = v
					)
				}
			}, error = function(e){
				list(
					b = rep(NA_real_, ncol(X_fit) + 1),
					vcov = matrix(NA_real_, ncol(X_fit) + 1, ncol(X_fit) + 1)
				)
			})
		}
	)
)
