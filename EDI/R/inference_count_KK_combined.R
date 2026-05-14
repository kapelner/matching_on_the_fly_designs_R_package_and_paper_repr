#' GLMM Inference for KK Designs with Count Response
#'
#' Fits a Poisson GLMM for count responses under a KK matching-on-the-fly design.
#' The random intercept per matched pair is integrated out via Gauss-Hermite quadrature.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an internal
#' Rcpp routine. Set \code{use_rcpp = FALSE} to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountKKGLMM = R6::R6Class("InferenceCountKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
		#' @description Initialize a KK Poisson GLMM inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the internal Rcpp Poisson GLMM.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
			private$use_rcpp = use_rcpp
		}
	),
	private = list(
		use_rcpp = TRUE,
		glmm_response_type = function() "count",
		glmm_family        = function() stats::poisson(link = "log"),
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		shared = function(estimate_only = FALSE){
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				super$shared(estimate_only)
			}
		},
		shared_rcpp = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()
			private$cached_mod = NULL
			private$cached_values$likelihood_test_context = NULL
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			if (ncol(as.matrix(private$X)) > 0){
				X_fit = private$create_design_matrix()
			} else {
				X_fit = cbind(`(Intercept)` = 1, w = private$w)
			}
			X_fit = as.matrix(X_fit)
			j_T = 1L  # 0-based index of treatment column
			fit = tryCatch(
				fast_poisson_glmm_cpp(
					X        = X_fit,
					y        = as.numeric(private$y),
					group_id = as.integer(group_id),
					j_T      = j_T,
					start_par = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					smart_start = private$smart_default,
					estimate_only    = estimate_only,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)

			if (is.null(fit) || !isTRUE(fit$converged)) {
				# Rcpp failed; fall back to glmmTMB
				return(super$shared(estimate_only = estimate_only))
			}
			beta_hat_T = as.numeric(fit$b[j_T + 1L])
			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				return(super$shared(estimate_only = estimate_only))
			}
			private$cached_mod = fit
			private$set_fit_warm_start(as.numeric(c(fit$b, fit$log_sigma)), "params", fisher = fit$fisher_information)
			private$cached_values$likelihood_test_context = list(
				X = X_fit,
				y = as.numeric(private$y),
				group_id = as.integer(group_id),
				j_treat = j_T + 1L,
				start = as.numeric(c(fit$b, fit$log_sigma))
			)
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		},
		get_likelihood_test_spec = function(){
			if (!isTRUE(private$use_rcpp)) return(NULL)
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			group_id = as.integer(ctx$group_id)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				group_id = group_id,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					fast_poisson_glmm_cpp(
						X = X_fit,
						y = y,
						group_id = group_id,
						start_par = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start,
						warm_start_fisher_info = private$get_fit_warm_start_fisher(length(ctx$start)),
						smart_start = private$smart_default,
						j_T = j_treat - 1L,
						estimate_only = FALSE,
						n_gh = 20L,
						maxit = 300L,
						eps_g = 1e-6,
						fixed_idx = j_treat,
						fixed_values = delta,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){
					as.numeric(c(fit$b, fit$log_sigma))
				},
				score = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.numeric(get_poisson_glmm_score_cpp(X_fit, y, group_id, params))
				},
				observed_information = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.matrix(get_poisson_glmm_hessian_cpp(X_fit, y, group_id, params))
				},
				fisher_information = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.matrix(get_poisson_glmm_hessian_cpp(X_fit, y, group_id, params))
				},
				information = function(fit){
					params = as.numeric(c(fit$b, fit$log_sigma))
					as.matrix(get_poisson_glmm_hessian_cpp(X_fit, y, group_id, params))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		}
	)
)
#' KK Hurdle Poisson Combined-Likelihood Inference for Count Responses
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with count
#' responses using a joint likelihood over all subjects.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountKKHurdlePoissonOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountKKHurdlePoissonOneLik = R6::R6Class("InferenceCountKKHurdlePoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKHurdlePoissonOneLik,
	public = list(
	)
)
