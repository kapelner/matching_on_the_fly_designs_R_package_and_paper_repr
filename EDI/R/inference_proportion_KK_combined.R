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
#' @export
InferencePropKKGEE = R6::R6Class("InferencePropKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				if (!check_package_installed("geepack")){
					stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	),
	private = list(
		gee_response_type = function() "proportion",
		gee_family        = function() stats::binomial(link = "logit")
	)
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
#' @export
InferencePropKKGLMM = R6::R6Class("InferencePropKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
			private$use_rcpp = use_rcpp
		}
	),
	private = list(
		use_rcpp = TRUE,
		glmm_response_type  = function() "proportion",
		glmm_family         = function() stats::binomial(link = "logit"),

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

			fit = tryCatch(
				fast_logistic_glmm_cpp(
					X        = X_fit,
					y        = as.numeric(private$y),
					group_id = as.integer(group_id),
					j_T      = j_T,
					estimate_only = estimate_only
				),
				error = function(e) NULL
			)

			if (is.null(fit) || !isTRUE(fit$converged)) {
				private$cache_nonestimable_estimate("kk_glmm_rcpp_failed")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			beta_hat_T = as.numeric(fit$b[j_T + 1L])

			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				private$cache_nonestimable_estimate("kk_glmm_rcpp_nonestimable")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
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
