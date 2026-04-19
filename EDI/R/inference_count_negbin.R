#' Negative Binomial Regression Inference for Count Responses
#'
#' Fits a negative binomial regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountNegBin = R6::R6Class("InferenceCountNegBin",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a negative binomial regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose)
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
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()

			if (length(Xmm_cols) == 0L){
				Xmm = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				Xmm = cbind(1, treatment = private$w, X_cov)
			}

			res = fast_negbin_regression(Xmm = Xmm, y = private$y)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			as.numeric(res$b[2])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		generate_mod = function(estimate_only = FALSE){
			if (private$include_covariates) {
				if (estimate_only && !is.null(private$best_Xmm_colnames)) {
					X_data = private$get_X()
					Xmm_cols = private$best_Xmm_colnames
					if (length(Xmm_cols) == 0L) {
						Xmm = cbind(1, private$w)
					} else {
						X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
						Xmm = cbind(1, treatment = private$w, X_cov)
					}
					res = fast_negbin_regression(Xmm = Xmm, y = private$y)
					if (is.null(res) || !is.finite(res$b[2])) return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_))
					return(list(b = res$b, ssq_b_2 = NA_real_))
				}

				X_full = private$create_design_matrix()
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = match("treatment", colnames(X_full)),
					fit_fun = function(X_fit, j_treat){
						if (estimate_only) {
							res = fast_negbin_regression(Xmm = X_fit, y = private$y)
							list(b = res$b, ssq_b_2 = NA_real_, j_treat = j_treat)
						} else {
							res = fast_negbin_regression_with_var(Xmm = X_fit, y = private$y)
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
			} else {
				Xmm = cbind(1, private$w)
				full_names = c("(Intercept)", "treatment")
				colnames(Xmm) = full_names[seq_len(ncol(Xmm))]

				if (estimate_only) {
					res = fast_negbin_regression(Xmm = Xmm, y = private$y)
					list(b = res$b, ssq_b_2 = NA_real_)
				} else {
					fast_negbin_regression_with_var(Xmm = Xmm, y = private$y)
				}
			}
		}
	)
)
