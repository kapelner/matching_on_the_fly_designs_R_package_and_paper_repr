#' Poisson Regression Inference for Count Responses
#'
#' Fits a Poisson log-link regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceCountPoisson = R6::R6Class("InferenceCountPoisson",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a Poisson regression inference object.
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
		},

		#' @description
		#' Computes the Poisson-regression estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
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

		fit_poisson_with_var = function(Xmm, j_treat, estimate_only = FALSE){
			mod = tryCatch(
				if (estimate_only) {
					fast_poisson_regression_cpp(Xmm, private$y)
				} else {
					fast_poisson_regression_with_var_cpp(Xmm, private$y, j = j_treat)
				},
				error = function(e) NULL
			)
			if (is.null(mod) || (is.logical(mod$converged) && !isTRUE(mod$converged))){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			coef_hat = as.numeric(mod$b)
			ssq_b_j = if (estimate_only) NA_real_ else as.numeric(mod$ssq_b_j)
			
			list(b = coef_hat, ssq_b_2 = ssq_b_j)
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
					required_cols = match("treatment", colnames(X_full)),
					fit_fun = function(X_fit, j_treat){
						res = private$fit_poisson_with_var(X_fit, j_treat = j_treat, estimate_only = estimate_only)
						res$j_treat = j_treat
						res
					},
					fit_ok = function(mod, X_fit, keep){
						j_treat = mod$j_treat
						if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
						if (estimate_only) return(TRUE)
						is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
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
				private$fit_poisson_with_var(Xmm, j_treat = 2L, estimate_only = estimate_only)
			}
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			if (is.null(w_mat)) return(NULL)
			X_covars = if (private$include_covariates) private$get_X() else NULL
			log_transform = transform_responses == "log"
			compute_poisson_distr_parallel_cpp(as.numeric(y), X_covars, w_mat, as.numeric(delta), log_transform, private$n_cpp_threads(ncol(w_mat)))
		}
	)
)
