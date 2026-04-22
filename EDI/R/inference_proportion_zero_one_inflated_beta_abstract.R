#' Zero/One-Inflated Beta Inference for Proportion Responses
#'
#' Internal base class for non-KK zero/one-inflated beta regression models. The
#' response is modeled as a three-component mixture with point masses at 0 and 1
#' plus a beta-distributed interior component on \eqn{(0, 1)}. The reported
#' treatment effect is the treatment coefficient from the beta mean submodel, on
#' the logit scale.
#'
#' @details
#' The inflation probabilities for exact 0 and exact 1 are intercept-only in this
#' first implementation. The beta mean submodel uses treatment alone in the
#' univariate class and treatment plus covariates in the multivariate class.
#'
#' @keywords internal
#' @noRd
InferencePropZeroOneInflatedBetaAbstract = R6::R6Class("InferencePropZeroOneInflatedBetaAbstract",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a zero-one-inflated beta regression inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),

	private = list(
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			Xmm = if (length(Xmm_cols) == 0L) {
				cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				cbind(`(Intercept)` = 1, treatment = private$w, X_data[, Xmm_cols, drop = FALSE])
			}
			
			res = tryCatch(fast_zero_one_inflated_beta_cpp(Xmm, private$y), error = function(e) NULL)
			if (is.null(res) || length(res$b) < 2 || !is.finite(res$b[2])) return(NA_real_)
			as.numeric(res$b[2])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		generate_mod = function(estimate_only = FALSE){
			# Use the common GLM fitting pattern
			attempt = private$fit_with_hardened_qr_column_dropping(
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					res = fast_zero_one_inflated_beta_cpp(X_fit, private$y)
					res$j_treat = j_treat
					res
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), c("(Intercept)", "treatment"))
			}
			attempt$fit
		}
	)
)
