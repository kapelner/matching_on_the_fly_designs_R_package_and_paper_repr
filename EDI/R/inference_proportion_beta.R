#' Beta Regression Inference for Proportion Responses
#'
#' Fits a beta regression for proportion responses (constrained to (0, 1)) using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferencePropBetaRegr = R6::R6Class("InferencePropBetaRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
		#' @description
		#' Initialize a beta-regression inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
				assertFormula(model_formula, null.ok = TRUE)
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

			if (length(Xmm_cols) == 0L){
				Xmm = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				Xmm = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}

			res = fast_beta_regression_cpp(X = Xmm, y = as.numeric(private$y))
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			as.numeric(res$b[2])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit){
					if (estimate_only) {
						res = fast_beta_regression_cpp(X_fit, private$y)
						list(b = res$b, ssq_b_2 = NA_real_)
					} else {
						fast_beta_regression_with_var_cpp(X_fit, private$y)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
			}
			attempt$fit
		},

		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		}
	)
)
