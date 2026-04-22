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
#' @export
InferenceContinKKGLMM = R6::R6Class("InferenceContinKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
		#' @description
		#' Initialize a KK GLMM inference object.
		#' @param des_obj A completed \code{Design} object with a continuous response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the optimised Rcpp
		#'   Gaussian LMM implementation (no external package required). If \code{FALSE},
		#'   use \pkg{glmmTMB}.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg The optimization algorithm to use.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = "lbfgs"){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			# If using Rcpp, skip glmmTMB package check in the parent initialize.
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE, default = "lbfgs")
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)

			private$use_rcpp = use_rcpp
		}
	),
	private = list(
		use_rcpp = TRUE,

		glmm_response_type = function() "continuous",
		glmm_family        = function() stats::gaussian(link = "identity"),

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
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

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
					y             = as.numeric(private$y),
					X             = X_fit,
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

			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$is_z = TRUE
			private$cached_values$df   = Inf

			if (estimate_only) return(invisible(NULL))

			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		}
	)
)
