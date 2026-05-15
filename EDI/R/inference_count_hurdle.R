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
		#' @description Initialize a hurdle Poisson inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula Optional formula for the count submodel.
		#' @param model_formula_hurdle Formula for the hurdle submodel. Defaults to
		#'   \code{~ .}, meaning treatment plus all available covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the internal Rcpp
		#'   implementation. If \code{FALSE}, use \pkg{glmmTMB}.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, model_formula_hurdle = ~ ., use_rcpp = TRUE, verbose = FALSE, smart_default = TRUE, optimization_alg = NULL){
			super$initialize(des_obj, model_formula = model_formula, model_formula_zero = model_formula_hurdle, use_rcpp = use_rcpp, verbose = verbose, smart_default = smart_default, optimization_alg = optimization_alg)
		}
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
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param model_formula_hurdle Formula for the hurdle submodel. Defaults to
		#'   \code{~ .}, meaning treatment plus all available covariates.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj, model_formula = NULL, model_formula_hurdle = ~ ., verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFormula(model_formula_hurdle, null.ok = FALSE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$model_formula_hurdle = model_formula_hurdle
		},
		#' @description Compute asymp confidence interval
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
		#' @description Compute asymp two sided pval for treatment effect
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
		#' @description Compute gradient / likelihood-based alternatives
		#' @param delta The null treatment effect (default 0).
		compute_gradient_two_sided_pval = function(delta = 0){
			private$compute_likelihood_test_two_sided_pval(delta = delta, testing_type = "gradient")
		},
		#' @description Compute likelihood-based confidence interval
		#' @param alpha The significance level (default 0.05).
		compute_gradient_confidence_interval = function(alpha = 0.05){
			private$invert_test_pval_confidence_interval(alpha)
		}
	),
	private = list(
		hurdle_description = function() "Hurdle Negative Binomial",
		cached_mod = NULL,
		model_formula_hurdle = NULL,
		best_hurdle_X_colnames = NULL,
		build_component_matrix = function(model_formula, selected_colnames = NULL, treatment_name = "treatment"){
			if (is.null(selected_colnames)) {
				if (identical(model_formula, ~ .)) {
					X_cov = private$get_X()
				} else {
					X_imp = private$des_obj$get_X_imp()
					X_cov = if (is.null(X_imp)) matrix(NA_real_, nrow = private$n, ncol = 0) else create_model_matrix_from_features(model_formula, X_imp)
				}
			} else {
				if (identical(model_formula, ~ .)) {
					X_cov_all = private$get_X()
				} else {
					X_imp = private$des_obj$get_X_imp()
					X_cov_all = if (is.null(X_imp)) matrix(NA_real_, nrow = private$n, ncol = 0) else create_model_matrix_from_features(model_formula, X_imp)
				}
				X_cov = if (is.null(X_cov_all) || !length(selected_colnames)) matrix(NA_real_, nrow = private$n, ncol = 0) else as.matrix(X_cov_all[, intersect(selected_colnames, colnames(X_cov_all)), drop = FALSE])
			}
			if (is.null(X_cov) || ncol(as.matrix(X_cov)) == 0L) {
				X_fit = cbind(1, private$w)
				colnames(X_fit) = c("(Intercept)", treatment_name)
				return(X_fit)
			}
			X_fit = cbind(1, private$w, as.matrix(X_cov))
			colnames(X_fit)[1:2] = c("(Intercept)", treatment_name)
			X_fit
		},
		predictors_df = function(){
			data.frame(w = private$w)
		},
		try_hurdle_negbin_fit = function(X_f, X_hurdle, j_t, estimate_only = FALSE){
			n_params = ncol(X_f) + 1L
			warm_start_params = private$get_fit_warm_start_for_length("params", n_params)
			warm_fisher = private$get_fit_warm_start_fisher(n_params)
			mod = tryCatch(
				if (estimate_only) {
					fast_hurdle_negbin_cpp(
						X_f, private$y, X_hurdle = X_hurdle,
						warm_start_params = warm_start_params,
						smart_start = private$smart_default,
						warm_start_fisher_info = warm_fisher
					)
				} else {
					fast_hurdle_negbin_with_var_cpp(
						X_f, private$y, X_hurdle = X_hurdle, j = j_t,
						warm_start_params = warm_start_params,
						smart_start = private$smart_default,
						warm_start_fisher_info = warm_fisher
					)
				},
				error = function(e) NULL
			)
			if (is.null(mod)) return(NULL)
			b = as.numeric(mod$b)
			ssq = if (estimate_only) NA_real_ else as.numeric(mod$ssq_b_j)
			if (length(b) != ncol(X_f) || any(!is.finite(b))) return(NULL)
			if (!estimate_only && (!is.finite(ssq) || ssq < 0)) return(NULL)
			list(mod = mod, b = b, ssq = ssq, j = j_t, X = X_f, X_hurdle = X_hurdle)
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$count_likelihood_context
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			X_hurdle = ctx$X_hurdle
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			opt_alg = private$optimization_alg %||% "lbfgs"
			n_params = ncol(X_fit) + 1L
			full_fit = tryCatch(
				fast_truncated_negbin_count_cpp(
					X = X_fit,
					y = y,
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					estimate_only = FALSE,
					optimization_alg = opt_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_fit) || !isTRUE(full_fit$converged)) return(NULL)
				list(
					X = X_fit,
					X_hurdle = X_hurdle,
					y = y,
					j = j_treat,
					full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					warm_start_params = start %||% private$get_fit_warm_start_for_length("params", n_params)
					warm_fisher = private$get_fit_warm_start_fisher(n_params)
					fast_truncated_negbin_count_cpp(
						X = X_fit,
						y = y,
						warm_start_params = warm_start_params,
						warm_start_fisher_info = warm_fisher,
						smart_start = private$smart_default,
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
					as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% -get_hurdle_negbin_count_hessian_cpp(X_fit, y, as.numeric(fit$params %||% c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))))
				},
				fisher_information = function(fit){
					as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% -get_hurdle_negbin_count_hessian_cpp(X_fit, y, as.numeric(fit$params %||% c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))))
				},
				information = function(fit){
					as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% -get_hurdle_negbin_count_hessian_cpp(X_fit, y, as.numeric(fit$params %||% c(as.numeric(fit$b), log(as.numeric(fit$theta_hat))))))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_component_matrix(private$model_formula, treatment_name = "treatment")
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$cache_nonestimable_estimate("hurdle_negbin_design_unusable")
				return(NULL)
			}
			colnames(X_fit) = colnames(X_full)[reduced$keep]
			if (is.null(private$best_hurdle_X_colnames)) {
				X_hurdle_full = private$build_component_matrix(private$model_formula_hurdle, treatment_name = "treatment")
				reduced_hurdle = private$reduce_design_matrix_preserving_treatment(X_hurdle_full)
				X_hurdle = reduced_hurdle$X
				if (is.null(X_hurdle)) {
					private$cache_nonestimable_estimate("hurdle_negbin_aux_design_unusable")
					return(NULL)
				}
				colnames(X_hurdle) = colnames(X_hurdle_full)[reduced_hurdle$keep]
				private$best_hurdle_X_colnames = setdiff(colnames(X_hurdle), c("(Intercept)", "treatment"))
			} else {
				X_hurdle = private$build_component_matrix(private$model_formula_hurdle, private$best_hurdle_X_colnames, treatment_name = "treatment")
			}
			fit = private$try_hurdle_negbin_fit(X_fit, X_hurdle, j_treat, estimate_only = estimate_only)
			if (private$harden && is.null(fit) && ncol(X_fit) > 2L){
				X_treat_only = X_fit[, 1:2, drop = FALSE]
				X_hurdle_treat_only = X_hurdle[, 1:2, drop = FALSE]
				fit = private$try_hurdle_negbin_fit(X_treat_only, X_hurdle_treat_only, 2L, estimate_only = estimate_only)
				reduced = list(X = X_treat_only, keep = 1:2, j_treat = 2L)
			}
			if (is.null(fit)){
				private$cache_nonestimable_estimate("hurdle_negbin_fit_unavailable")
				return(NULL)
			}
			
			private$cached_values$count_likelihood_context = list(X = fit$X, X_hurdle = fit$X_hurdle, j_treat = fit$j)
			
			out = list(
				mod = fit$mod,
				beta_hat_T = fit$b[fit$j],
				ssq_b_j = fit$ssq,
				params = c(as.numeric(fit$b), log(as.numeric(fit$mod$theta_hat))),
				fisher_information = fit$mod$fisher_information
			)
			
			# Cache extra values for this specific class
			b_full = rep(NA_real_, ncol(X_full))
			b_full[reduced$keep] = fit$b
			names(b_full) = colnames(X_full)
			private$cached_values$full_coefficients = b_full
			private$cached_values$theta_hat = as.numeric(fit$mod$theta_hat)
			private$cached_values$hurdle_coefficients = fit$mod$hurdle_b
			
			out
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
