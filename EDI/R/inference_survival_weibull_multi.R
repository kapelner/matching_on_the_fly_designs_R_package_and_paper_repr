#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- InferenceSurvivalMultiWeibullRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalMultiWeibullRegr = R6::R6Class("InferenceSurvivalMultiWeibullRegr",
	lock_objects = FALSE,
	inherit = InferenceSurvivalUniWeibullRegr,
	public = list(
		#' @description
		#' Initialize the Inference object.
		#'
		#' @param des_obj The design object.
		#' @param verbose If TRUE, print additional information.
		initialize = function(des_obj, verbose = FALSE) {
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose)
		}

	),

	private = list(
		generate_mod = function(estimate_only = FALSE){
			# Multivariate: covariates + treatment (mirrors InferenceSurvivalMultiCoxPHRegr)
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

			# We already have a complex fallback loop here, but we need to populate best_Xmm_colnames
			# for the optimized randomization/bootstrap path.
			
			thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			full_X_matrix_last = NULL
			for (thresh in thresholds) {
				if (is.finite(thresh)) {
					X_cov = drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M
				} else {
					X_cov = X_cov_orig
				}
				full_X_matrix = cbind(X_cov, treatment = private$w)
				full_X_matrix_last = full_X_matrix
				mod = tryCatch(private$weibull_generate_mod_from_X(full_X_matrix, estimate_only = estimate_only), error = function(e) NULL)
				if (!is.null(mod)) {
					private$best_Xmm_colnames = setdiff(colnames(full_X_matrix), "treatment")
					return(mod)
				}
			}
			# All fast C++ paths failed: one robust survreg fallback on the most-reduced matrix
			surv_mod = robust_survreg_with_surv_object(survival::Surv(private$y, private$dead), full_X_matrix_last)
			if (!is.null(surv_mod)) {
				private$best_Xmm_colnames = setdiff(colnames(full_X_matrix_last), "treatment")
				full_coefficients = c(surv_mod$coefficients, "log(scale)" = log(surv_mod$scale))
				if (estimate_only) {
					return(list(coefficients = full_coefficients, vcov = NULL))
				}
				full_vcov = surv_mod$var
				if (!is.null(full_vcov) && is.matrix(full_vcov) && all(is.finite(diag(full_vcov)))) {
					colnames(full_vcov) = rownames(full_vcov) = names(full_coefficients)
					return(list(coefficients = full_coefficients, vcov = full_vcov))
				}
			}
			# Progressive QR-ordered column dropping
			X_cov_last = full_X_matrix_last[, seq_len(ncol(full_X_matrix_last) - 1L), drop = FALSE]
			if (ncol(X_cov_last) > 0L) {
				keep_js = qr(X_cov_last)$pivot
				while (length(keep_js) > 0L) {
					keep_js = keep_js[-length(keep_js)]
					if (length(keep_js) > 0L) {
						X_try = cbind(X_cov_last[, keep_js, drop = FALSE], treatment = private$w)
					} else {
						X_try = matrix(private$w, ncol = 1L, dimnames = list(NULL, "treatment"))
					}
					mod = tryCatch(private$weibull_generate_mod_from_X(X_try, estimate_only = estimate_only), error = function(e) NULL)
					if (!is.null(mod)) {
						private$best_Xmm_colnames = setdiff(colnames(X_try), "treatment")
						return(mod)
					}
				}
			}
			stop("Weibull regression failed to converge even after progressive reduction.")
		}
	)
)
