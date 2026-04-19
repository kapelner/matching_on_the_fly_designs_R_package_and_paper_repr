#' Modified Poisson Regression Inference for Incidence Responses
#'
#' Fits a modified Poisson regression (Zou 2004) for binary (incidence) responses
#' using the treatment indicator and, optionally, all recorded covariates as
#' predictors. This model provides an alternative to log-binomial regression for
#' estimating risk ratios.
#'
#' @export
InferenceIncidModifiedPoisson = R6::R6Class("InferenceIncidModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a modified Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
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

			res = fast_poisson_regression_cpp(X = Xmm, y = as.numeric(private$y))
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
					res = fast_poisson_regression_cpp(X = Xmm, y = as.numeric(private$y))
					if (is.null(res) || !is.finite(res$b[2])) return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_))
					return(list(b = res$b, ssq_b_2 = NA_real_))
				}

				X_full = private$create_design_matrix()
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = 1L,
					fit_fun = function(X_fit, j_treat){
						if (estimate_only) {
							res = fast_poisson_regression_cpp(X_fit, private$y)
							list(b = res$b, ssq_b_2 = NA_real_, j_treat = j_treat)
						} else {
							res = fast_poisson_regression_with_var_cpp(X_fit, private$y, j = j_treat)
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
					res = fast_poisson_regression_cpp(Xmm, private$y)
					list(b = res$b, ssq_b_2 = NA_real_)
				} else {
					fast_poisson_regression_with_var_cpp(Xmm, private$y)
				}
			}
		}
	)
)
