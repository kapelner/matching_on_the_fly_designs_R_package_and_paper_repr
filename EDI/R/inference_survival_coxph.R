#' Cox Proportional Hazards Regression Inference for Survival Responses
#'
#' Fits a Cox proportional hazards regression for survival responses using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalCoxPHRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalCoxPHRegr = R6::R6Class("InferenceSurvivalCoxPHRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
				
		#' @description Initialize a Cox PH inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use internal Rcpp.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values by default.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
			
			
			private$use_rcpp = use_rcpp
		},
		#' @description Computes the Cox PH estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			super$compute_estimate(estimate_only = estimate_only)
		},
		#' @description Recomputes the Cox PH treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) return(beta_hat_T)
			}
			X_cov = private$get_X()
			X_fit = if (!is.null(X_cov) && ncol(X_cov) > 0) cbind(treatment = private$w, X_cov) else matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			fit = weighted_cox_bootstrap_surrogate_fit(private$y, private$dead, X_fit, row_weights)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		}
	),
	private = list(
		use_rcpp = TRUE,
		cox_data_cache = NULL,
		cox_w_cache = NULL,
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$b)
			if (!all(is.finite(b_null))) return(NULL)
			X_fit = spec$X
			y_obs = as.numeric(private$y)
			dead_obs = as.numeric(private$dead)
			breslow = .breslow_hazard(y_obs, dead_obs, X_fit, b_null)
			if (length(breslow$times) == 0L) return(NULL)
			sim = .cox_simulate_from_breslow(breslow, y_obs, dead_obs, X_fit, b_null)
			y_sim = sim$y_sim; dead_sim = sim$dead_sim
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			j = spec$j
			full_res = tryCatch(
				fast_coxph_regression_cpp(
					X = X_fit, y = y_sim, dead = dead_sim,
					estimate_only = FALSE
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$coefficients[j])) return(NULL)
			full_fit_boot = list(b = as.numeric(full_res$coefficients), neg_loglik = as.numeric(full_res$neg_ll))
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = tryCatch(
						fast_coxph_regression_cpp(
							X = X_fit, y = y_sim, dead = dead_sim,
							warm_start_beta = start %||% as.numeric(full_res$coefficients),
							fixed_idx = j, fixed_values = d,
							estimate_only = FALSE
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), neg_loglik = as.numeric(res$neg_ll))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X = ctx$X
			y = as.numeric(private$y)
			dead = as.numeric(private$dead)
			full_b_cox = as.numeric(private$cached_mod$b[-1])  # strip 0-prefix
			full_fit = list(b = full_b_cox, neg_loglik = ctx$full_neg_loglik)
			list(
				X = X, y = y, j = 1L,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					res = tryCatch(
						fast_coxph_regression_cpp(
							X, y, dead,
							warm_start_beta = start %||% private$get_fit_warm_start_for_length("beta", ncol(X)),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X)),
							fixed_idx = 1L,
							fixed_values = delta,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = TRUE
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), neg_loglik = as.numeric(res$neg_ll), fisher_information = res$fisher_information)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_coxph_score_cpp(X, y, dead, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_coxph_hessian_cpp(X, y, dead, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_coxph_hessian_cpp(X, y, dead, as.numeric(fit$b))
				},
				information = function(fit){
					-get_coxph_hessian_cpp(X, y, dead, as.numeric(fit$b))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_cov = private$get_X()
			X_fit = if (!is.null(X_cov) && ncol(X_cov) > 0){
				cbind(treatment = private$w, X_cov)
			} else {
				matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			}
			if (private$use_rcpp) {
				if (is.null(private$cox_data_cache) || !identical(private$w, private$cox_w_cache)) {
					private$cox_data_cache = build_cox_data_cache_cpp(X_fit, private$y, private$dead)
					private$cox_w_cache = private$w
				}
				optimization_alg = .normalize_optimizer_algorithm(NULL, allow_irls = FALSE, default = "newton_raphson")
				raw_res = tryCatch(
					fast_coxph_regression_prebuilt_cpp(
						private$cox_data_cache,
						estimate_only = estimate_only,
						optimization_alg = optimization_alg,
						warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
						smart_cold_start = private$smart_cold_start_default
					),
					error = function(e) NULL
				)
				fit = if (!is.null(raw_res) && isTRUE(raw_res$converged)) {
					b = as.numeric(raw_res$coefficients)
					names(b) = colnames(X_fit)
					list(b = b, coefficients = b, vcov = if (estimate_only) NULL else raw_res$vcov,
					     neg_log_lik = as.numeric(raw_res$neg_ll), fisher_information = raw_res$fisher_information)
				} else NULL
				fit = if (is.null(fit)) tryCatch(
					fast_coxph_regression(
						X_fit, private$y, private$dead,
						use_rcpp = TRUE,
						estimate_only = estimate_only,
						warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
						smart_cold_start = private$smart_cold_start_default
					),
					error = function(e) NULL
				) else fit
				if (is.null(fit)) {
					private$cached_values$likelihood_test_context = NULL
					return(list(b = rep(NA_real_, ncol(X_fit)), vcov = matrix(NA_real_, ncol(X_fit), ncol(X_fit))))
				}
				private$set_fit_warm_start(as.numeric(fit$b), "beta", fisher = fit$fisher_information)
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					full_neg_loglik = fit$neg_log_lik
				)
				# fast_coxph_regression returns coefficients with intercept=FALSE for Cox.
				# InferenceAsympLikStdModCache expects intercept in first position if it exists.
				# But Cox has no intercept. So we prefix 0.
				return(list(
					b = c(0, fit$b),
					ssq_b_2 = if (estimate_only) NA_real_ else fit$vcov[1, 1],
					vcov = if (estimate_only) NULL else {
						v = matrix(0, ncol(X_fit) + 1, ncol(X_fit) + 1)
						v[2:(ncol(X_fit) + 1), 2:(ncol(X_fit) + 1)] = fit$vcov
						v
					}
				))
			}
			surv_obj = survival::Surv(private$y, private$dead)
			tryCatch({
				coxph_mod = suppressWarnings(survival::coxph(surv_obj ~ X_fit))
				if (estimate_only) {
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = NA_real_,
						vcov = NULL
					)
				} else {
					vcov_mat = stats::vcov(coxph_mod)
					v = matrix(0, ncol(X_fit) + 1, ncol(X_fit) + 1)
					v[2:(ncol(X_fit) + 1), 2:(ncol(X_fit) + 1)] = vcov_mat
					list(
						b = c(0, stats::coef(coxph_mod)),
						ssq_b_2 = as.numeric(vcov_mat[1, 1]),
						vcov = v
					)
				}
			}, error = function(e){
				list(
					b = rep(NA_real_, ncol(X_fit) + 1),
					vcov = matrix(NA_real_, ncol(X_fit) + 1, ncol(X_fit) + 1)
				)
			})
		}
	)
)
