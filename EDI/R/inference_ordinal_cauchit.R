#' Cumulative Cauchit Inference for Ordinal Responses
#'
#' Cumulative Cauchit model inference for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceOrdinalCauchitRegr = R6::R6Class("InferenceOrdinalCauchitRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
						#' @description
		#' Initialize a cumulative-cauchit inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),

	private = list(
		cauchit_polr_fallback = function(){
			if (!check_package_installed("MASS")) return(NULL)
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)
			mod = tryCatch(
				MASS::polr(y ~ w, data = dat, method = "cauchit", Hess = TRUE),
				error = function(e) NULL
			)
			if (is.null(mod) || !"w" %in% names(stats::coef(mod))) return(NULL)
			coef_w = as.numeric(stats::coef(mod)["w"])
			var_w = tryCatch(vcov(mod)["w", "w"], error = function(e) NA_real_)
			ssq = if (is.finite(var_w) && var_w > 0) var_w else NA_real_
			list(b = c(NA, coef_w), ssq_b_2 = ssq)
		},

		generate_mod = function(estimate_only = FALSE){
			# Use the common GLM fitting pattern
			attempt = private$try_fit_full_and_harden(
				fit_fun = function(X_fit, keep){
					res = fast_ordinal_cauchit_regression_with_var_cpp(X = X_fit, y = as.numeric(private$y))
					# Add intercept placeholder to b for MLEorKM logic
					res$b = c(NA_real_, res$b)
					res
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (!is.null(mod$converged) && !mod$converged) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			attempt$fit
		}
	)
)
