#' Linear Mixed Model Inference for KK Designs with Continuous Response
#'
#' Fits a linear mixed model for continuous responses under a KK
#' matching-on-the-fly design. The matched-pair strata enter as a subject-level
#' random intercept \code{(1 | group_id)}, accounting for within-pair correlation.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an
#' internal Rcpp/L-BFGS routine that requires no external packages. Set
#' \code{use_rcpp = FALSE} to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'continuous')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rnorm(10))
#' inf = InferenceContinKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceContinKKGLMM = R6::R6Class("InferenceContinKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
		#' @description Initialize a KK GLMM inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the optimised Rcpp
		#'   Gaussian LMM implementation (no external package required). If \code{FALSE},
		#'   use \pkg{glmmTMB}.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg The optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			# If using Rcpp, skip glmmTMB package check in the parent initialize.
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
			private$use_rcpp = use_rcpp
		}
	),
	private = list(
		use_rcpp = TRUE,
		glmm_response_type = function() "continuous",
		glmm_family        = function() stats::gaussian(link = "identity"),
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		# glmmTMB path: univariate → treatment only
		glmm_predictors_df = function(){
			super$glmm_predictors_df()
		},
		# ── Dispatch ─────────────────────────────────────────────────────────
		shared = function(estimate_only = FALSE){
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				super$shared(estimate_only)
			}
		},
		# ── Rcpp Gaussian LMM path ────────────────────────────────────────────
		shared_rcpp = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			private$clear_nonestimable_state()
			private$cached_mod = NULL
			private$cached_values$likelihood_test_context = NULL
			# ── Group IDs (same logic as fit_glmm_on_data) ───────────────────
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			# ── Design matrix: (intercept, w, covariates) ────────────────────
			X_fit = private$create_design_matrix()
			# create_design_matrix uses "treatment"; Rcpp path expects "w" in the search for j_T_r
			if ("treatment" %in% colnames(X_fit))
				colnames(X_fit)[colnames(X_fit) == "treatment"] = "w"
			X_fit = as.matrix(X_fit)
			# Treatment column index (1-based in R, 0-based in C++ is handled via ssq_b_T)
			j_T_r = which(colnames(X_fit) == "w")
			if (length(j_T_r) == 0L) j_T_r = 2L   # fallback: second column
			fit = tryCatch(
				fast_gaussian_lmm_cpp(
					X             = X_fit,
					y             = as.numeric(private$y),
					group_id      = as.integer(group_id),
					estimate_only = estimate_only,
					maxit         = 300L,
					eps_g         = 1e-6,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				return(super$shared(estimate_only = estimate_only))
			}
			beta_hat_T = as.numeric(fit$b[j_T_r])
			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				return(super$shared(estimate_only = estimate_only))
			}
			private$cached_mod = fit
			private$set_fit_warm_start(as.numeric(fit$b), "params")
			private$cached_values$likelihood_test_context = list(
				X = X_fit,
				y = as.numeric(private$y),
				group_id = as.integer(group_id),
				j_treat = j_T_r,
				start = as.numeric(fit$b)
			)
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		}
		,
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
					start_par = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start
					fast_gaussian_lmm_cpp(
						X = X_fit,
						y = y,
						group_id = group_id,
						start_par = start_par,
						estimate_only = FALSE,
						maxit = 300L,
						eps_g = 1e-6,
						fixed_idx = j_treat,
						fixed_values = delta,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){
					as.numeric(fit$b)
				},
				score = function(fit){
					as.numeric(get_gaussian_lmm_score_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				observed_information = function(fit){
					as.matrix(get_gaussian_lmm_fisher_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				fisher_information = function(fit){
					as.matrix(get_gaussian_lmm_fisher_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				information = function(fit){
					as.matrix(get_gaussian_lmm_fisher_cpp(X_fit, y, group_id, as.numeric(fit$b)))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		}
	)
)
