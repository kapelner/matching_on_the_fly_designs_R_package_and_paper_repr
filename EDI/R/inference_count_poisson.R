#' Poisson Regression Inference for Count Responses
#'
#' Fits a Poisson log-link regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountPoisson$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountPoisson = R6::R6Class("InferenceCountPoisson",
	lock_objects = FALSE,
	inherit = InferenceCountLikelihood,
	public = list(
				
		#' @description Initialize a Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_default Whether to use smart optimizer start values by default.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			X_data = private$get_X()
			if (is.null(X_data) || ncol(X_data) == 0) {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					res = fast_poisson_regression_weighted_cpp(
						X = X_fit,
						y = private$y,
						weights = row_weights,
						warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
						smart_start = private$smart_default
					)
					list(b = res$b, XtWX = res$XtWX, ssq_b_j = NA_real_, j_treat = j_treat)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) && length(mod$b) >= mod$j_treat && is.finite(mod$b[mod$j_treat])
				}
			)
			private$cached_mod = attempt$fit
			if (is.null(attempt$fit) || is.null(attempt$fit$b)) {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(NA_real_)
			}
			j_treat = attempt$fit$j_treat %||% 2L
			private$cached_values$beta_hat_T = as.numeric(attempt$fit$b[j_treat])
			private$cached_values$s_beta_hat_T = NA_real_
			private$set_fit_warm_start(
				as.numeric(attempt$fit$b),
				"beta",
				fisher = attempt$fit$XtWX
			)
			private$cached_values$beta_hat_T
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
				X = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(1, treatment = private$w, X_cov)
			}
			res = tryCatch(
				fast_poisson_regression_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X)),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X)),
					smart_start = private$smart_default
				),
				error = function(e) NULL
			)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$XtWX)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		supports_fisher_information = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					fast_poisson_regression_with_var_cpp(
						X = X_fit,
						y = y,
						j = j_treat,
						warm_start_beta = start %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit)),
						fixed_idx = j_treat,
						fixed_values = delta,
						smart_start = private$smart_default
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					get_poisson_regression_score_cpp(X_fit, y, as.numeric(fit$b))
				},
				observed_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				fisher_information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				information = function(fit){
					-get_poisson_regression_hessian_cpp(X_fit, as.numeric(fit$b))
				},
				neg_loglik = function(fit){
					eta = as.numeric(X_fit %*% as.numeric(fit$b))
					-sum(y * eta - exp(eta) - lgamma(y + 1))
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			# Use the common GLM fitting pattern
			X_data = private$get_X()
			if (is.null(X_data) || ncol(X_data) == 0) {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w)
			} else {
				X_full = cbind(`(Intercept)` = 1, treatment = private$w, X_data)
			}
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit))
					warm_fisher = private$get_fit_warm_start_fisher(ncol(X_fit))
					if (estimate_only) {
						res = fast_poisson_regression_cpp(
							X = X_fit, y = private$y,
							warm_start_beta = warm_start_beta,
							warm_start_fisher_info = warm_fisher,
							smart_start = private$smart_default
						)
						list(b = res$b, XtWX = res$XtWX, ssq_b_j = NA_real_, j_treat = j_treat)
					} else {
						res = fast_poisson_regression_with_var_cpp(
							X = X_fit, y = private$y, j = j_treat,
							warm_start_beta = warm_start_beta,
							warm_start_fisher_info = warm_fisher,
							smart_start = private$smart_default
						)
						res$j_treat = j_treat
						res
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod)) return(FALSE)
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = which(attempt$keep == 2L)
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			if (is.null(w_mat)) return(NULL)
			X_covars = private$X
			log_transform = transform_responses == "log"
			compute_poisson_distr_parallel_cpp(X_covars, as.numeric(y), w_mat, as.numeric(delta), log_transform, private$n_cpp_threads(ncol(w_mat)))
		}
	)
)
