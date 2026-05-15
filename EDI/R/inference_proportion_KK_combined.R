#' GEE Inference for KK Designs with Proportion Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for proportion (continuous values in (0, 1)) responses under a KK
#' matching-on-the-fly design using the treatment indicator and, optionally,
#' all recorded covariates as predictors.
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed in Suggests
#' and is not installed automatically with \pkg{EDI}.
#' Install \pkg{geepack} before using this class.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropKKGEE = R6::R6Class("InferencePropKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = c(InferenceMixinKKGEEShared$public, list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param use_rcpp Whether to use the internal Rcpp solver (TRUE) or fallback to geepack (FALSE).
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts() && !use_rcpp) {
				if (!check_package_installed("geepack")){
					stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
			private$init_kk_gee_shared(des_obj, use_rcpp = use_rcpp, model_formula = model_formula)
		}
	)),
	private = c(InferenceMixinKKGEEShared$private, list(
		gee_response_type = function() "proportion",
		gee_family        = function() stats::binomial(link = "logit"),
		shared_gee_dispatch = function(estimate_only = FALSE) private$shared_gee_default(estimate_only)
	))
)
#' GLMM Inference for KK Designs with Proportion Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) for proportion (continuous values
#' in (0, 1)) responses under a KK matching-on-the-fly design. The random intercept
#' per matched pair is integrated out via Gauss-Hermite quadrature.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an internal
#' Rcpp/L-BFGS routine that requires no external packages. Set \code{use_rcpp = FALSE}
#' to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropKKGLMM = R6::R6Class("InferencePropKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = c(InferenceMixinKKGLMMShared$public, list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_default = smart_default)
			private$init_kk_glmm_shared(des_obj)
			private$use_rcpp = use_rcpp
		}
	)),
	private = c(InferenceMixinKKGLMMShared$private, list(
		use_rcpp = TRUE,
		glmm_response_type  = function() "proportion",
		glmm_family         = function() stats::binomial(link = "logit"),
		supports_likelihood_tests = function() isTRUE(private$use_rcpp),
		shared = function(estimate_only = FALSE){
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				private$shared_glmm_tmb(estimate_only)
			}
		},
			shared_rcpp = function(estimate_only = FALSE){
				if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
				if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
				private$clear_nonestimable_state()
				private$cached_values$likelihood_test_context = NULL
				m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			# X WITH intercept; treatment is at 0-based index j_T = 1 (second column)
			if (ncol(as.matrix(private$X)) > 0){
				X_fit = private$create_design_matrix()
			} else {
				X_fit = cbind(`(Intercept)` = 1, w = private$w)
			}
			X_fit = as.matrix(X_fit)
			j_T = 1L  # 0-based index of treatment column
			y_vals = as.numeric(private$y)
			finite_y = y_vals[is.finite(y_vals)]
			if (length(finite_y) == 0L || any(finite_y < 0) || any(finite_y > 1)) {
				return(private$shared_glmm_tmb(estimate_only = estimate_only))
			}
			n_params = ncol(X_fit) + 1L
			fit = tryCatch(
				fast_logistic_glmm_cpp(
					X        = X_fit,
					y        = y_vals,
					group_id = as.integer(group_id),
					j_T      = j_T,
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					smart_start = private$smart_default,
					estimate_only    = estimate_only,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				return(private$shared_glmm_tmb(estimate_only = estimate_only))
			}
			beta_hat_T = as.numeric(fit$b[j_T + 1L])
				if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
					return(private$shared_glmm_tmb(estimate_only = estimate_only))
				}
				private$cached_mod = fit
				private$set_fit_warm_start(as.numeric(fit$params), "params", fisher = fit$fisher_information)
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					y = as.numeric(private$y),
					group_id = as.integer(group_id),
					j_T = j_T,
					j_treat = j_T + 1L,
					n_gh = 20L,
					start = as.numeric(fit$params)
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
					fast_logistic_glmm_cpp(
						X = X_fit,
						y = y,
						group_id = group_id,
						warm_start_params = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start,
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
					as.numeric(fit$params)
				},
				score = function(fit){
					as.numeric(get_logistic_glmm_score_cpp(X_fit, y, group_id, as.numeric(fit$params)))
				},
				observed_information = function(fit){
					as.matrix(get_logistic_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(fit$params)))
				},
				fisher_information = function(fit){
					as.matrix(get_logistic_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(fit$params)))
				},
				information = function(fit){
					as.matrix(get_logistic_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(fit$params)))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		}
	))
)
