#' GEE Inference for KK Designs with Ordinal Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{multgee})
#' for ordinal responses under a KK matching-on-the-fly design using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKGEE$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKGEE = R6::R6Class("InferenceOrdinalKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				if (!check_package_installed("multgee")){
					stop("Package 'multgee' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	),
	private = list(
		gee_response_type = function() "ordinal",

		# Override: ordinal response requires ordLORgee, not geeglm.
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			pred_df = private$gee_predictors_df()
			dat = data.frame(y = factor(private$y, ordered = TRUE), pred_df, group_id = group_id)
			dat = dat[order(dat$group_id), ]
			id_sorted = dat$group_id
			
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			formula_gee = stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))

			mod = tryCatch({
				utils::capture.output(m <- suppressMessages(suppressWarnings(
					multgee::ordLORgee(
						formula_gee,
						data   = dat,
						id     = id_sorted,
						LORstr = "uniform",
						link   = "logit"
					)
				)))
				m
			}, error = function(e) NULL)

			if (is.null(mod)){
				private$cache_nonestimable_estimate("ordinal_kk_gee_fit_unavailable")
				return(invisible(NULL))
			}

			beta = stats::coef(mod)
			j_treat = private$gee_treatment_index(beta)
			private$cached_values$beta_hat_T = as.numeric(beta[j_treat])

			if (estimate_only) return(invisible(NULL))

			vcov_robust = tryCatch(stats::vcov(mod), error = function(e) NULL)
			if (is.null(vcov_robust)) {
				private$cached_values$s_beta_hat_T = NA_real_
			} else {
				private$cached_values$s_beta_hat_T = sqrt(as.numeric(vcov_robust[j_treat, j_treat]))
			}
			private$cached_values$df = Inf
			private$cached_values$summary_table = summary(mod)$coefficients
		}
	)
)

#' GLMM Inference for KK Designs with Ordinal Response
#'
#' Fits a cumulative-logit mixed model (proportional odds) for ordinal responses
#' under a KK matching-on-the-fly design. The random intercept per matched pair is
#' integrated out via Gauss-Hermite quadrature.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an internal
#' Rcpp/L-BFGS routine that requires no external packages. Set \code{use_rcpp = FALSE}
#' to fall back to \pkg{glmmTMB}.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKGLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKGLMM = R6::R6Class("InferenceOrdinalKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use internal Rcpp.
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
		glmm_response_type  = function() "ordinal",
		glmm_family         = function() glmmTMB::cumulative(link = "logit"),

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

			# X WITHOUT intercept (cutpoints serve as intercepts)
			if (ncol(as.matrix(private$X)) > 0){
				X_fit = as.matrix(private$glmm_predictors_df())  # [w, cov1, ...]
			} else {
				X_fit = matrix(private$w, ncol = 1L, dimnames = list(NULL, "w"))
			}

			# Convert y to 1-indexed integers
			y_levels = sort(unique(private$y))
			K = length(y_levels)
			y = as.integer(match(private$y, y_levels))
			n_alpha = K - 1L

			# Treatment is always the first column of X_fit (j_T = 0, 0-based)
			j_T = 0L

			# Warm start from fixed-effects ordinal MLE to avoid divergence
			start = tryCatch({
				nore = fast_ordinal_regression_cpp(X_fit, as.numeric(y) - 1L)
				alpha_direct = as.numeric(nore$alpha)  # K-1 direct cutpoints
				beta_nore    = as.numeric(nore$b)      # p betas
				# Convert direct alphas to log-diff parameterization
				alpha_par = numeric(n_alpha)
				if (n_alpha >= 1L) alpha_par[1L] = alpha_direct[1L]
				if (n_alpha >= 2L) {
					for (k in 2L:n_alpha) {
						diff_k = alpha_direct[k] - alpha_direct[k - 1L]
						alpha_par[k] = if (diff_k > 0) log(diff_k) else 0.0
					}
				}
				c(alpha_par, beta_nore, -3.0)  # log_sigma = -3 (small random effect)
			}, error = function(e) NULL)

			fit = tryCatch(
				fast_ordinal_glmm_cpp(
					X          = X_fit,
					y          = y,
					group_id   = as.integer(group_id),
					K          = K,
					j_T        = j_T,
					estimate_only = estimate_only,
					start      = start,
					eps_g      = 1e-3
				),
				error = function(e) NULL
			)

			if (is.null(fit) || !isTRUE(fit$converged)) {
				private$cache_nonestimable_estimate("kk_glmm_rcpp_failed")
				return(invisible(NULL))
			}

			# b is the beta vector (no cutpoints); treatment is at index j_T+1 (1-based R)
			beta_hat_T = as.numeric(fit$b[j_T + 1L])

			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				private$cache_nonestimable_estimate("kk_glmm_rcpp_nonestimable")
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf

			if (estimate_only) return(invisible(NULL))

			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		}
	)
)
