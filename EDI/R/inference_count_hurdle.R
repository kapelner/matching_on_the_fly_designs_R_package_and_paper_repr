#' Hurdle Poisson Regression Inference for Count Responses
#'
#' Fits a hurdle Poisson regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountHurdleNegBin$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountHurdlePoisson = R6::R6Class("InferenceCountHurdlePoisson",
	lock_objects = FALSE,
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	public = list(
	),
	private = list(
		za_family = function() glmmTMB::truncated_poisson(link = "log"),
		za_description = function() "Hurdle Poisson"
	)
)

#' Hurdle Negative Binomial Regression Inference for Count Responses
#'
#' Fits a hurdle negative binomial regression for count responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountHurdleNegBin$new(seq_des, model_formula = ~ x1)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountHurdleNegBin = R6::R6Class("InferenceCountHurdleNegBin",
	lock_objects = FALSE,
	inherit = InferenceCountLikelihood,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute asymp two sided pval for treatment effect
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Compute gradient / likelihood-based alternatives
		#' @param delta The null treatment effect (default 0).
		compute_gradient_two_sided_pval = function(delta = 0){
			self$compute_asymp_two_sided_pval(delta = delta)
		},

		#' @description
		#' Compute likelihood-based confidence interval
		#' @param alpha The significance level (default 0.05).
		compute_gradient_confidence_interval = function(alpha = 0.05){
			self$compute_asymp_confidence_interval(alpha = alpha)
		}
	),
	private = list(
		hurdle_description = function() "Hurdle Negative Binomial",
		cached_mod = NULL,

		predictors_df = function(){
			data.frame(w = private$w)
		},

		try_hurdle_negbin_fit = function(X_f, j_t, estimate_only = FALSE){
			mod = tryCatch(
				if (estimate_only) {
					fast_hurdle_negbin_cpp(X_f, private$y)
				} else {
					fast_hurdle_negbin_with_var_cpp(X_f, private$y, j = j_t)
				},
				error = function(e) NULL
			)
			if (is.null(mod)) return(NULL)
			b = as.numeric(mod$b)
			ssq = if (estimate_only) NA_real_ else as.numeric(mod$ssq_b_j)
			if (length(b) != ncol(X_f) || any(!is.finite(b))) return(NULL)
			if (!estimate_only && (!is.finite(ssq) || ssq < 0)) return(NULL)
			list(mod = mod, b = b, ssq = ssq, j = j_t, X = X_f)
		},

		supports_likelihood_tests = function(){
			TRUE
		},

		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$count_likelihood_context
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			opt_alg = private$optimization_alg %||% "lbfgs"
			full_fit = tryCatch(
				fast_truncated_negbin_count_cpp(
					X = X_fit,
					y = y,
					start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L),
					estimate_only = FALSE,
					optimization_alg = opt_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_fit) || !isTRUE(full_fit$converged)) return(NULL)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					start_params = start %||% private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L)
					fast_truncated_negbin_count_cpp(
						X = X_fit,
						y = y,
						start_params = start_params,
						estimate_only = FALSE,
						optimization_alg = opt_alg,
						fixed_idx = j_treat,
						fixed_values = delta
					)
				},
				extract_start = function(fit){
					as.numeric(fit$params %||% c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))
				},
				score = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))
					as.numeric(fit$score %||% get_hurdle_negbin_count_score_cpp(X_fit, y, params))
				},
				observed_information = function(fit){
					as.matrix(fit$information %||% fit$observed_information)
				},
				fisher_information = function(fit){
					as.matrix(fit$information %||% fit$observed_information)
				},
				information = function(fit){
					as.matrix(fit$information %||% fit$observed_information)
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			X_full = cbind(1, private$w, as.matrix(private$predictors_df()[, setdiff(colnames(private$predictors_df()), "w"), drop = FALSE]))
			colnames(X_full)[1:2] = c("(Intercept)", "treatment")
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$cache_nonestimable_estimate("hurdle_negbin_design_unusable")
				return(invisible(NULL))
			}

			fit = private$try_hurdle_negbin_fit(X_fit, j_treat, estimate_only = estimate_only)
			if (private$harden && is.null(fit) && ncol(X_fit) > 2L){
				X_treat_only = X_fit[, 1:2, drop = FALSE]
				fit = private$try_hurdle_negbin_fit(X_treat_only, 2L, estimate_only = estimate_only)
				reduced = list(X = X_treat_only, keep = 1:2, j_treat = 2L)
			}
			if (is.null(fit)){
				private$cache_nonestimable_estimate("hurdle_negbin_fit_unavailable")
				return(invisible(NULL))
			}

			private$cached_mod = fit$mod
			private$cached_values$count_likelihood_context = list(X = fit$X, j_treat = fit$j)
			private$set_fit_warm_start(c(as.numeric(fit$b), log(as.numeric(fit$mod$theta_hat))), "params")
			b_full = rep(NA_real_, ncol(X_full))
			b_full[reduced$keep] = fit$b
			names(b_full) = colnames(X_full)

			private$cached_values$beta_hat_T = fit$b[fit$j]
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(fit$ssq)
			private$cached_values$full_coefficients = b_full
			private$cached_values$theta_hat = as.numeric(fit$mod$theta_hat)
			private$cached_values$hurdle_coefficients = fit$mod$hurdle_b
		},

		get_standard_error = function(){
			private$shared()
			as.numeric(private$cached_values$s_beta_hat_T)
		},

		get_degrees_of_freedom = function() Inf,

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)
)
