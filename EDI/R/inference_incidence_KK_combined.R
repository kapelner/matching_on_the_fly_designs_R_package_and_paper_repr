#' GEE Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for binary (incidence) responses under a KK matching-on-the-fly design using
#' the treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKGEE = R6::R6Class("InferenceIncidKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param use_rcpp Whether to use the internal Rcpp solver (TRUE) or fallback to geepack (FALSE).
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			if (should_run_asserts() && !use_rcpp) {
				if (!check_package_installed("geepack")){
					stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, use_rcpp = use_rcpp)
		}
	),
	private = list(
		gee_response_type = function() "incidence",
		gee_family        = function() stats::binomial(link = "logit")
	)
)

#' GLMM Inference for KK Designs with Binary Response
#'
#' Fits a Generalized Linear Mixed Model (GLMM) for binary (incidence) responses
#' under a KK matching-on-the-fly design. The random intercept per matched pair is
#' integrated out via Gauss-Hermite quadrature.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an internal
#' Rcpp/L-BFGS routine that requires no external packages. Set \code{use_rcpp = FALSE}
#' to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKGLMM = R6::R6Class("InferenceIncidKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Whether to use the internal Rcpp solver.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			if (use_rcpp) private$skip_glmm_pkg_check = TRUE
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			private$use_rcpp = use_rcpp
		}
	),
	private = list(
		use_rcpp = TRUE,
		glmm_response_type  = function() "incidence",
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

			fit = tryCatch(
				fast_logistic_glmm_cpp(
					X        = X_fit,
					y        = as.numeric(private$y),
					group_id = as.integer(group_id),
					j_T      = j_T,
					estimate_only    = estimate_only,
					eps_g            = 1e-4,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)

				if (is.null(fit) || !isTRUE(fit$converged)) {
					return(super$shared(estimate_only = estimate_only))
				}

			# b[j_T+1] in R 1-based = treatment coefficient
			beta_hat_T = as.numeric(fit$b[j_T + 1L])

				if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
					return(super$shared(estimate_only = estimate_only))
				}

				private$cached_mod = fit
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					y = as.numeric(private$y),
					group_id = as.integer(group_id),
					j_T = j_T,
					j_treat = j_T + 1L,
					n_gh = 20L
				)
				private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf

			if (estimate_only) return(invisible(NULL))

			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		}
	)
)
