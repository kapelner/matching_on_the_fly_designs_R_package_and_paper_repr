#' Dependent-Censoring Transformation Inference for Survival Responses
#'
#' Fits a survival model accounting for dependent censoring via a transformation
#' approach using the treatment indicator and, optionally, all recorded covariates
#' as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalDepCensTransformRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalDepCensTransformRegr = R6::R6Class("InferenceSurvivalDepCensTransformRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a dependent-censoring transformation inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
		}
	),
	private = list(
		best_X_colnames = NULL,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_X_colnames
			X_data = private$get_X()
			if (length(X_cols) == 0L){
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			n_params = 2 * ncol(X) + 3L
			warm_start_params = private$get_fit_warm_start_for_length("params", n_params)
			if (is.null(warm_start_params)) warm_start_params = rep(0, n_params)
			warm_fisher = private$get_fit_warm_start_fisher(n_params)
			res = fast_dep_cens_transform_optim_cpp(
				y = private$y, dead = private$dead, X = X,
				warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
				warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
				smart_start = private$smart_default
			)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "params", fisher = res$fisher_information)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			dead = as.numeric(private$dead)
			j_treat = 2L  # treatment is second param = beta_event[treatment]
			n_params = 2L * ncol(X_fit) + 3L
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm_start_params = start %||% private$get_fit_warm_start_for_length("params", n_params)
					warm_fisher = private$get_fit_warm_start_fisher(n_params)
					tryCatch(fast_dep_cens_transform_optim_cpp(
						y = y, dead = dead, X = X_fit, warm_start_params = warm_start_params,
						warm_start_fisher_info = warm_fisher,
						smart_start = private$smart_default,
						fixed_idx = j_treat, fixed_values = delta
					), error = function(e) NULL)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_dep_cens_transform_score_cpp(X_fit, y, dead, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_dep_cens_transform_hessian_cpp(X_fit, y, dead, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_dep_cens_transform_hessian_cpp(X_fit, y, dead, as.numeric(fit$b))
				},
				information = function(fit){
					-get_dep_cens_transform_hessian_cpp(X_fit, y, dead, as.numeric(fit$b))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit){
					n_params = 2 * ncol(X_fit) + 3L
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params)
					warm_fisher = private$get_fit_warm_start_fisher(n_params)
					res = tryCatch(fast_dep_cens_transform_optim_cpp(
						y = private$y, dead = private$dead, X = X_fit, warm_start_params = warm_start_params,
						warm_start_fisher_info = warm_fisher,
						smart_start = private$smart_default
					), error = function(e) NULL)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(
						b = as.numeric(res$b),
						neg_loglik = as.numeric(res$neg_loglik),
						fisher_information = res$fisher_information,
						ssq_b_2 = if (estimate_only || is.null(res$vcov) || nrow(res$vcov) < 2L) NA_real_
						          else as.numeric(res$vcov[2L, 2L])
					)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$b, "params", fisher = attempt$fit$fisher_information)
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(X = attempt$X)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},
		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X = cbind(`(Intercept)` = 1, treatment = private$w, X_cov)
			}
			X
		}
	)
)
