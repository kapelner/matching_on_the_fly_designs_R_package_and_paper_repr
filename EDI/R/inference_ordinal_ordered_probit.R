#' Ordered Probit Inference for Ordinal Responses
#'
#' Ordinal probit (ordered probit) model inference for ordinal responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceOrdinalOrderedProbitRegr = R6::R6Class("InferenceOrdinalOrderedProbitRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
		#' @description
		#' Initialize an ordered-probit inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE, harden = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose, harden)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
		}
	),

	private = list(
		include_covariates = NULL,

		probit_polr_fallback = function(){
			if (!check_package_installed("MASS")) return(NULL)
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)
			mod = tryCatch(
				MASS::polr(y ~ w, data = dat, method = "probit", Hess = TRUE),
				error = function(e) NULL
			)
			if (is.null(mod) || !"w" %in% names(stats::coef(mod))) return(NULL)
			coef_w = as.numeric(stats::coef(mod)["w"])
			var_w = tryCatch(vcov(mod)["w", "w"], error = function(e) NA_real_)
			ssq = if (is.finite(var_w) && var_w > 0) var_w else NA_real_
			list(b = c(NA, coef_w), ssq_b_2 = ssq)
		},

		generate_mod = function(estimate_only = FALSE){
			if (private$include_covariates) {
				X_full = cbind(treatment = private$w, private$get_X())
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = 1L,
					fit_fun = function(X_fit){
						res = fast_ordinal_probit_regression_with_var_cpp(X = X_fit, y = as.numeric(private$y))
						b1 = tryCatch(res$b[1], error = function(e) NA_real_)
						if (is.finite(b1)){
							if (estimate_only) return(list(b = c(NA, res$b), ssq_b_2 = NA_real_, converged = res$converged))
							if (is.finite(res$ssq_b_2) && res$ssq_b_2 > 0) return(list(b = c(NA, res$b), ssq_b_2 = res$ssq_b_2, converged = res$converged))
						}
						fallback = private$probit_polr_fallback()
						if (!is.null(fallback)){
							if (estimate_only) return(list(b = fallback$b, ssq_b_2 = NA_real_, converged = TRUE))
							return(c(fallback, list(converged = TRUE)))
						}
						list(b = c(NA, res$b), ssq_b_2 = if (estimate_only) NA_real_ else res$ssq_b_2, converged = res$converged)
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
						if (!is.null(mod$converged) && !mod$converged) return(FALSE)
						if (estimate_only) return(TRUE)
						is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
					}
				)
				return(attempt$fit)
			} else {
				Xmm = matrix(private$w, ncol = 1)
				colnames(Xmm) = "treatment"

				res = fast_ordinal_probit_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))
				b1 = tryCatch(res$b[1], error = function(e) NA_real_)
				if (is.finite(b1)){
					if (estimate_only) return(list(b = c(NA, b1), ssq_b_2 = NA_real_))
					if (is.finite(res$ssq_b_2) && res$ssq_b_2 > 0) return(list(b = c(NA, b1), ssq_b_2 = res$ssq_b_2))
				}
				fallback = private$probit_polr_fallback()
				if (!is.null(fallback)){
					if (estimate_only) return(list(b = fallback$b, ssq_b_2 = NA_real_))
					return(fallback)
				}
				return(list(b = c(NA, b1), ssq_b_2 = if (estimate_only) NA_real_ else res$ssq_b_2))
			}
		}
	)
)
