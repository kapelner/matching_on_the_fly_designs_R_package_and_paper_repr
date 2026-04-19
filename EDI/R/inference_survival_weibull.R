#' Weibull Regression Inference for Survival Responses
#'
#' Fits a Weibull regression for survival responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceSurvivalWeibullRegr = R6::R6Class("InferenceSurvivalWeibullRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMSummaryTable,
	public = list(

		#' @description
		#' Initialize a Weibull regression inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the optimized Rcpp
		#'   implementation. If \code{FALSE}, use \pkg{survival::survreg}.
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
		}
	),

	private = list(
		include_covariates = NULL,
		best_Xmm_colnames = NULL,
		use_rcpp = TRUE,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
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

			mod = tryCatch(private$weibull_generate_mod_from_X(Xmm, estimate_only = estimate_only), error = function(e) NULL)
			if (is.null(mod) || !("treatment" %in% names(mod$coefficients))){
				return(NA_real_)
			}
			as.numeric(mod$coefficients["treatment"])
		},

		generate_mod = function(estimate_only = FALSE){
			if (private$include_covariates) {
				X_cov_orig = private$get_X()
				if (ncol(X_cov_orig) == 0L) {
					full_X_matrix = matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
					mod = tryCatch(private$weibull_generate_mod_from_X(full_X_matrix, estimate_only = estimate_only), error = function(e) NULL)
					if (should_run_asserts()) {
						if (is.null(mod)) stop("Weibull regression failed to converge (univariate).")
					}
					private$best_Xmm_colnames = character(0)
					return(mod)
				}

				# Multivariate fallback loop for near-collinearity
				thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
				full_X_matrix_last = NULL
				for (thresh in thresholds) {
					if (is.finite(thresh)) {
						X_cov = drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M
					} else {
						X_cov = X_cov_orig
					}
					full_X_matrix = cbind(treatment = private$w, X_cov)
					
					mod = tryCatch(private$weibull_generate_mod_from_X(full_X_matrix, estimate_only = estimate_only), error = function(e) NULL)
					if (!is.null(mod)) {
						private$best_Xmm_colnames = colnames(X_cov)
						return(mod)
					}
				}
				stop("Weibull regression failed to converge even after robust retries.")
			} else {
				# Univariate logic
				full_X_matrix = matrix(private$w, ncol = 1)
				colnames(full_X_matrix) = "treatment"
				mod = tryCatch(private$weibull_generate_mod_from_X(full_X_matrix, estimate_only = estimate_only), error = function(e) NULL)
				if (!is.null(mod)) {
					private$best_Xmm_colnames = character(0)
					return(mod)
				}
				stop("Weibull regression failed to converge.")
			}
		},

		weibull_generate_mod_from_X = function(full_X_matrix, estimate_only = FALSE){
			weibull_regr_mod = fast_weibull_regression(
				private$y,
				private$dead,
				as.matrix(full_X_matrix),
				use_rcpp = private$use_rcpp,
				estimate_only = estimate_only
			)
			if (should_run_asserts()) {
				if (is.null(weibull_regr_mod$coefficients) || (!estimate_only && (is.null(weibull_regr_mod$vcov) || !is.matrix(weibull_regr_mod$vcov)))){
					stop("fast_weibull_regression failed to return valid coefficients or vcov.")
				}
			}
			full_coefficients = c(weibull_regr_mod$coefficients, "log(scale)" = weibull_regr_mod$log_sigma)
			if (estimate_only) {
				return(list(coefficients = full_coefficients, vcov = NULL))
			}
			full_vcov = weibull_regr_mod$vcov
			colnames(full_vcov) = rownames(full_vcov) = names(full_coefficients)
			list(coefficients = full_coefficients, vcov = full_vcov)
		}
	)
)
