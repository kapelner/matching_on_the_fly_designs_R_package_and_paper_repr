#' Univariate Modified Poisson Inference for Binary Responses
#'
#' Fits the classic modified Poisson model for binary (incidence) responses using
#' a Poisson log-link working model with Huber-White
#' sandwich variance. The treatment effect is reported on the log-risk-ratio
#' scale.
#'
#' @details
#' This is the practical risk-ratio approach popularized by Zou (2004). The
#' working model is Poisson with log link, but inference uses a robust sandwich
#' covariance rather than the naive Poisson model-based variance.
#'
#' @export
InferenceIncidUnivModifiedPoisson = R6::R6Class("InferenceIncidUnivModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a modified-Poisson inference object for a completed
		#' design with a binary response.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an
		#'   incidence response.
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "incidence")
		#' for (i in 1:20) {
		#' 	x_i = data.frame(x1 = rnorm(1), x2 = rnorm(1))
		#' 	w_i = seq_des$add_one_subject_to_experiment_and_assign(x_i)
		#' 	p_i = pmin(0.95, exp(-1 + 0.4 * w_i))
		#' 	seq_des$add_one_subject_response(i, rbinom(1, 1, p_i))
		#' }
		#' seq_des_inf = InferenceIncidUnivModifiedPoisson$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Computes the modified-Poisson estimate of the treatment effect on the
		#' log-risk-ratio scale.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = TRUE)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval using the sandwich
		#' standard error.
		#' @param alpha The confidence level in the computed confidence interval is
		#'   1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect on the log-risk-ratio scale.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		best_Xmm_colnames = NULL,

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
			
			if (length(Xmm_cols) == 0L){
				# Univariate case
				return(private$fit_univariate_modified_poisson(estimate_only = TRUE)$beta_hat)
			}

			# Multivariate case
			X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
			Xmm = cbind("(Intercept)" = 1, treatment = private$w, X_cov)

			fit = private$fit_modified_poisson(Xmm, j_treat = 2L, estimate_only = estimate_only)
			if (is.null(fit)) return(NA_real_)
			as.numeric(fit$beta_hat)
		},

		max_abs_reasonable_coef = 1e4,

		build_design_matrix = function(){
			cbind(1, private$w)
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		create_bootstrap_worker_state = function(){
			worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			worker$num_cores = 1L
			worker_priv = worker$.__enclos_env__$private
			X_fixed = private$get_X()
			worker_priv$X = X_fixed
			list(
				worker = worker,
				worker_priv = worker_priv,
				base_w = as.numeric(private$w),
				base_y = as.numeric(private$y),
				base_dead = as.numeric(private$dead),
				n = private$n
			)
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			indices = as.integer(indices)
			w_priv = worker_state$worker_priv
			w_priv$w[] = worker_state$base_w[indices]
			w_priv$y[] = worker_state$base_y[indices]
			w_priv$dead[] = worker_state$base_dead[indices]
			w_priv$y_temp = w_priv$y
			w_priv$n = worker_state$n
			w_priv$cached_values = list()
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			as.numeric(worker_state$worker$compute_treatment_estimate(estimate_only = TRUE))[1L]
		},

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},

		set_failed_fit_cache = function(){
			private$cache_nonestimable_estimate("modified_poisson_fit_unavailable")
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},

		coefficients_are_usable = function(coef_hat){
			length(coef_hat) > 0L &&
				all(is.finite(coef_hat)) &&
				max(abs(coef_hat), na.rm = TRUE) <= private$max_abs_reasonable_coef
		},

		fit_modified_poisson = function(X_fit, j_treat, estimate_only = FALSE){
			mod = tryCatch(
				fast_poisson_regression_cpp(X = X_fit, y = as.numeric(private$y)),
				error = function(e) NULL
			)
			if (is.null(mod)){
				return(NULL)
			}

			coef_hat = as.numeric(mod$b)
			if (length(coef_hat) != ncol(X_fit) || !private$coefficients_are_usable(coef_hat)){
				return(NULL)
			}

			if (estimate_only){
				return(list(
					beta_hat = coef_hat[j_treat],
					se = NA_real_,
					coefficients = coef_hat,
					vcov = NULL,
					summary_table = NULL
				))
			}

			mu_hat = as.numeric(exp(X_fit %*% coef_hat))
			if (length(mu_hat) != nrow(X_fit) || any(!is.finite(mu_hat)) || any(mu_hat <= 0)){
				return(NULL)
			}

			post_fit = tryCatch(
				glm_sandwich_post_fit_cpp(
					X_fit = X_fit,
					y = as.numeric(private$y),
					coef_hat = coef_hat,
					mu_hat = mu_hat,
					working_weights = mu_hat,
					j_treat = j_treat
				),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				return(NULL)
			}

			coef_names = colnames(X_fit)
			names(coef_hat) = coef_names
			vcov_robust = post_fit$vcov
			colnames(vcov_robust) = rownames(vcov_robust) = coef_names
			std_err = post_fit$std_err
			names(std_err) = coef_names
			z_vals = post_fit$z_vals
			names(z_vals) = coef_names

			list(
				beta_hat = post_fit$beta_hat,
				se = post_fit$se,
				coefficients = coef_hat,
				vcov = vcov_robust,
				summary_table = cbind(
					Value = coef_hat,
					`Std. Error` = std_err,
					`z value` = z_vals,
					`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
				)
			)
		},

		fit_univariate_modified_poisson = function(estimate_only = FALSE){
			X_uni = cbind("(Intercept)" = 1, treatment = as.numeric(private$w))
			private$fit_modified_poisson(X_uni, 2L, estimate_only = estimate_only)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T) && (estimate_only || !is.null(private$cached_values$summary_table))) return(invisible(NULL))

			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			if (is.null(colnames(X_full))) {
				colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)
			}

			reduced = private$reduce_design_matrix_preserving_treatment_fixed_covariates(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			fit = NULL
			fallback_to_univariate = FALSE
			if (nrow(X_fit) > ncol(X_fit)){
				fit = private$fit_modified_poisson(X_fit, j_treat, estimate_only = estimate_only)
			}
			if (is.null(fit)){
				fallback_fit = private$fit_univariate_modified_poisson(estimate_only = estimate_only)
				if (is.null(fallback_fit)){
					private$set_failed_fit_cache()
					return(invisible(NULL))
				}
				fit = fallback_fit
				fallback_to_univariate = TRUE
			}

			private$best_Xmm_colnames = if (fallback_to_univariate) character(0) else setdiff(colnames(X_fit), c("(Intercept)", "treatment"))
			private$cached_values$beta_hat_T = fit$beta_hat
			if (estimate_only) return(invisible(NULL))

			local_df = if (fallback_to_univariate) {
				nrow(cbind("(Intercept)" = 1, treatment = as.numeric(private$w))) - 2L
			} else {
				nrow(X_fit) - ncol(X_fit)
			}
			if (!is.finite(fit$se) || fit$se <= 0){
				if (!fallback_to_univariate){
					fallback_fit = private$fit_univariate_modified_poisson(estimate_only = FALSE)
					if (is.null(fallback_fit) || !is.finite(fallback_fit$se) || fallback_fit$se <= 0){
						private$set_failed_fit_cache()
						return(invisible(NULL))
					}
					fit = fallback_fit
					fallback_to_univariate = TRUE
				}
				local_df = nrow(cbind("(Intercept)" = 1, treatment = as.numeric(private$w))) - 2L
			}

			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = local_df
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$full_vcov = fit$vcov
			private$cached_values$summary_table = fit$summary_table
		}
	)
)
