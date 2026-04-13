#' Univariate Robust Poisson Regression Inference for Count Responses
#'
#' Fits a Poisson log-link regression for count responses using only the treatment
#' indicator. The treatment effect is reported on the log-rate scale and inference
#' uses a Huber-White sandwich variance rather than the model-based Poisson variance.
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
#'   response_type = "count",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- InferenceCountUnivRobustPoissonRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceCountUnivRobustPoissonRegr = R6::R6Class("InferenceCountUnivRobustPoissonRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "count")
			super$initialize(des_obj, verbose)
		},

		#' @description
		#' Computes the appropriate estimate.
		#'
		#' @return	The log-rate treatment-effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = TRUE)
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		build_design_matrix = function(){
			Xmm = cbind(1, private$w)
			full_names = c("(Intercept)", "treatment")
			colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
			Xmm
		},

		fit_count_model_with_var = function(Xmm, estimate_only = FALSE){
			reduced = private$reduce_design_matrix_preserving_treatment_fixed_covariates(Xmm)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			mod = tryCatch(
				fast_poisson_regression_cpp(X = X_fit, y = as.numeric(private$y)),
				error = function(e) NULL
			)
			if (is.null(mod)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			coef_hat = as.numeric(mod$b)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			if (estimate_only){
				b_full = rep(NA_real_, ncol(Xmm))
				b_full[reduced$keep] = coef_hat
				names(b_full) = colnames(Xmm)
				return(list(b = b_full, ssq_b_2 = NA_real_))
			}

			mu_hat = as.numeric(exp(X_fit %*% coef_hat))
			if (length(mu_hat) != nrow(X_fit) || any(!is.finite(mu_hat)) || any(mu_hat <= 0)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			bread = NULL
			var_keep = seq_len(ncol(X_fit))
			X_var = X_fit
			cross_mat = tryCatch(
				crossprod(X_fit, X_fit * mu_hat),
				error = function(e) NULL
			)
			if (!is.null(cross_mat)) {
				bread = tryCatch(solve(cross_mat), error = function(e) NULL)
			}
			if (is.null(bread)) {
				sqrt_mu = sqrt(mu_hat)
				X_weighted = X_fit * sqrt_mu
				reduced_weighted = private$reduce_design_matrix_preserving_treatment(X_weighted)
				keep_sub = as.integer(reduced_weighted$keep)
				if (length(keep_sub) == 0L) {
					return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
				}
				X_var_candidate = as.matrix(X_fit[, keep_sub, drop = FALSE])
				cross_mat = tryCatch(
					crossprod(X_var_candidate, X_var_candidate * mu_hat),
					error = function(e) NULL
				)
				if (!is.null(cross_mat)) {
					bread = tryCatch(solve(cross_mat), error = function(e) NULL)
					if (!is.null(bread)) {
						var_keep = keep_sub
						X_var = X_var_candidate
					}
				}
			}
			if (is.null(bread)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			resid = as.numeric(private$y) - mu_hat
			meat = crossprod(X_var, X_var * (resid^2))
			vcov_trim = bread %*% meat %*% bread
			vcov_full = matrix(NA_real_, ncol(X_fit), ncol(X_fit))
			if (length(var_keep) == ncol(X_fit)) {
				vcov_full = vcov_trim
			} else {
				for (i in seq_along(var_keep)){
					for (j in seq_along(var_keep)){
						vcov_full[var_keep[i], var_keep[j]] = vcov_trim[i, j]
					}
				}
			}
			vcov_robust = (vcov_full + t(vcov_full)) / 2

			ssq_b_2 = as.numeric(vcov_robust[j_treat, j_treat])
			if (!is.finite(ssq_b_2) || ssq_b_2 < 0){
				ssq_b_2 = NA_real_
			}

			b_full = rep(NA_real_, ncol(Xmm))
			b_full[reduced$keep] = coef_hat
			names(b_full) = colnames(Xmm)
			list(b = b_full, ssq_b_2 = ssq_b_2)
		},

		generate_mod = function(estimate_only = FALSE){
			private$fit_count_model_with_var(private$build_design_matrix(), estimate_only = estimate_only)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			has_cached_se = !is.null(private$cached_values$s_beta_hat_T) &&
				length(private$cached_values$s_beta_hat_T) == 1L &&
				is.finite(private$cached_values$s_beta_hat_T)
			if (!is.null(private$cached_values$is_z) && (estimate_only || has_cached_se)) return(invisible(NULL))
			model_output = private$generate_mod(estimate_only = estimate_only)
			private$cached_values$beta_hat_T = model_output$b[2]
			if (estimate_only) return(invisible(NULL))

			ssq = model_output$ssq_b_2
			if (!is.null(ssq) && !is.na(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
			private$cached_values$df = model_output$df %||% NA_real_
		}
	)
)
