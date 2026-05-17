#' Robust Poisson Regression Inference for Count Responses
#'
#' Fits a Poisson log-link regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors. Inference
#' uses a Huber-White sandwich variance.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountRobustPoisson$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountRobustPoisson = R6::R6Class("InferenceCountRobustPoisson",
	lock_objects = FALSE,
	inherit = InferenceCountCompositeLikelihood,
	public = list(
				
		#' @description Initialize a robust Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart starting values for the optimizer.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
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
					smart_cold_start = private$smart_cold_start_default,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X))
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
		build_design_matrix = function(){
			private$create_design_matrix()
		},
		fit_count_model_with_var = function(X, estimate_only = FALSE){
			reduced = private$reduce_design_matrix_preserving_treatment_fixed_covariates(X)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(b = rep(NA_real_, ncol(X)), ssq_b_2 = NA_real_, X_fit = X_fit, j_treat = j_treat))
			}
			mod = tryCatch(
				fast_poisson_regression_cpp(
					X = X_fit, y = as.numeric(private$y),
					warm_start_beta = private$get_fit_warm_start_for_length("beta", ncol(X_fit)),
					smart_cold_start = private$smart_cold_start_default,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit))
				),
				error = function(e) NULL
			)

			if (is.null(mod)){
				return(list(b = rep(NA_real_, ncol(X)), ssq_b_2 = NA_real_, X_fit = X_fit, j_treat = j_treat))
			}
			coef_hat = as.numeric(mod$b)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				return(list(b = rep(NA_real_, ncol(X)), ssq_b_2 = NA_real_, X_fit = X_fit, j_treat = j_treat))
			}
			
			if (estimate_only){
				b_full = rep(NA_real_, ncol(X))
				b_full[reduced$keep] = coef_hat
				names(b_full) = colnames(X)
				return(list(b = b_full, ssq_b_2 = NA_real_, X_fit = X_fit, j_treat = j_treat, mod = mod, XtWX = mod$XtWX))
			}
			mu_hat = as.numeric(exp(X_fit %*% coef_hat))
			if (length(mu_hat) != nrow(X_fit) || any(!is.finite(mu_hat)) || any(mu_hat <= 0)){
				return(list(b = b_full, ssq_b_2 = NA_real_, X_fit = X_fit, j_treat = j_treat, mod = mod, XtWX = mod$XtWX))
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
					return(list(b = rep(NA_real_, ncol(X)), ssq_b_2 = NA_real_, X_fit = X_fit, j_treat = j_treat, mod = mod, XtWX = mod$XtWX))
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
				return(list(b = rep(NA_real_, ncol(X)), ssq_b_2 = NA_real_, X_fit = X_fit, j_treat = j_treat, mod = mod, XtWX = mod$XtWX))
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
			ssq_b_j = as.numeric(vcov_robust[j_treat, j_treat])
			if (!is.finite(ssq_b_j) || ssq_b_j < 0){
				ssq_b_j = NA_real_
			}
			b_full = rep(NA_real_, ncol(X))
			b_full[reduced$keep] = coef_hat
			names(b_full) = colnames(X)
			list(b = b_full, ssq_b_2 = ssq_b_j, mod = mod, XtWX = mod$XtWX, X_fit = X_fit, j_treat = j_treat)
		},
		generate_mod = function(estimate_only = FALSE){
			model_output = private$fit_count_model_with_var(private$build_design_matrix(), estimate_only = estimate_only)
			if (!is.null(model_output$b)) {
				private$cached_values$likelihood_test_context = list(
					X = model_output$X_fit %||% private$build_design_matrix(),
					j_treat = model_output$j_treat %||% 2L
				)
				private$best_X_colnames = setdiff(colnames(model_output$X_fit), c("(Intercept)", "treatment"))
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			model_output
		}
	)
)
