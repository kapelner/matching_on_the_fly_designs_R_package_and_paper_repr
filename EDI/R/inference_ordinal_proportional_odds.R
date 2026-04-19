#' Proportional Odds Inference for Ordinal Responses
#'
#' Fits a proportional odds model (cumulative logit) for ordinal responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceOrdinalPropOddsRegr = R6::R6Class("InferenceOrdinalPropOddsRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a proportional-odds inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
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
			
			Xmm = if (length(Xmm_cols) == 0L){
				# Univariate case
				matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				cbind(treatment = private$w, X_cov)
			}

			res = fast_ordinal_regression_cpp(X = Xmm, y = as.numeric(private$y))
			if (is.null(res) || !is.finite(res$b[1])){
				return(NA_real_)
			}
			as.numeric(res$b[1])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		generate_mod = function(estimate_only = FALSE){
			# Proportional odds model logit(P(Y <= k)) = alpha_k - beta * w
			# We use fast_ordinal_regression_with_var_cpp
			# Xmm should NOT include intercept
			
			if (private$include_covariates) {
				if (estimate_only && !is.null(private$best_Xmm_colnames)) {
					X_data = private$get_X()
					Xmm_cols = private$best_Xmm_colnames
					Xmm = if (length(Xmm_cols) == 0L) {
						matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
					} else {
						X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
						cbind(treatment = private$w, X_cov)
					}
					res = fast_ordinal_regression_cpp(X = Xmm, y = as.numeric(private$y))
					if (is.null(res) || !is.finite(res$b[1])) return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_))
					return(list(b = c(NA, res$b[1]), ssq_b_2 = NA_real_))
				}

				X_full = private$get_X()
				if (ncol(X_full) == 0L) {
					Xmm = matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
					res = if (estimate_only) {
						r = fast_ordinal_regression_cpp(Xmm, as.numeric(private$y))
						list(b = c(NA, r$b[1]), ssq_b_2 = NA_real_)
					} else {
						r = fast_ordinal_regression_with_var_cpp(Xmm, as.numeric(private$y))
						list(b = c(NA, r$b[1]), ssq_b_2 = r$ssq_b_2)
					}
					private$best_Xmm_colnames = character(0)
					return(res)
				}

				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = cbind(treatment = private$w, X_full),
					required_cols = 1L, # treatment is first
					fit_fun = function(X_fit, j_treat){
						if (estimate_only) {
							res = fast_ordinal_regression_cpp(X_fit, as.numeric(private$y))
							list(b = c(NA, res$b[j_treat]), ssq_b_2 = NA_real_, j_treat = j_treat)
						} else {
							res = fast_ordinal_regression_with_var_cpp(X_fit, as.numeric(private$y))
							list(b = c(NA, res$b[j_treat]), ssq_b_2 = res$ssq_b_2, j_treat = j_treat)
						}
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod) || !is.finite(mod$b[2])) return(FALSE)
						if (estimate_only) return(TRUE)
						is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
					}
				)

				if (!is.null(attempt$fit)){
					private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), "treatment")
				}
				attempt$fit
			} else {
				Xmm = matrix(private$w, ncol = 1)
				colnames(Xmm) = "treatment"

				if (estimate_only) {
					res = fast_ordinal_regression_cpp(Xmm, as.numeric(private$y))
					list(b = c(NA, res$b[1]), ssq_b_2 = NA_real_)
				} else {
					res = fast_ordinal_regression_with_var_cpp(Xmm, as.numeric(private$y))
					list(b = c(NA, res$b[1]), ssq_b_2 = res$ssq_b_2)
				}
			}
		}
	)
)
