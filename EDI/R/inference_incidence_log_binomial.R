#' Log-Binomial Regression Inference for Incidence Responses
#'
#' Fits a log-binomial regression for binary (incidence) responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceIncidLogBinomial = R6::R6Class("InferenceIncidLogBinomial",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
				
		#' @description
		#' Initialize a log-binomial regression inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
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
				Xmm = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				Xmm = cbind(1, treatment = private$w, X_cov)
			}

			res = tryCatch(fast_log_binomial_regression_cpp(X = Xmm, y = as.numeric(private$y)), error = function(e) NULL)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
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
					if (estimate_only) {
						res = fast_log_binomial_regression_cpp(X = X_fit, y = private$y)
						list(b = res$b, ssq_b_j = NA_real_, j_treat = j_treat)
					} else {
						res = fast_log_binomial_regression_with_var_cpp(Xmm = X_fit, y = private$y, j = j_treat)
						res$j_treat = j_treat
						res
					}
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
