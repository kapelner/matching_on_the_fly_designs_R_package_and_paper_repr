#' Univariate Quantile Regression Inference for Continuous Responses
#'
#' Fits a quantile regression for continuous responses using only the treatment
#' indicator as a predictor. The treatment effect is reported on the response
#' scale at quantile \code{tau}; by default \code{tau = 0.5}, so this is median
#' regression.
#'
#' Standard errors use \pkg{quantreg}'s Powell-style \code{"nid"} estimator when
#' available, with fallback to \code{"iid"} if needed. Inference is based on the
#' resulting asymptotic normal approximation.
#'
#' This class requires the \pkg{quantreg} package, which is listed under
#' \code{Suggests} and is not installed automatically with \pkg{EDI}.
#' Install \pkg{quantreg} manually before use.
#'
#' @export
InferenceContinUnivQuantileRegr = R6::R6Class("InferenceContinUnivQuantileRegr",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a quantile-regression inference object for a completed design
		#' with a continuous response.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with a continuous response.
		#' @param tau The quantile level for regression, strictly between 0 and 1. The default
		#'   \code{tau = 0.5}
		#'   estimates the median treatment effect.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "continuous")
		#' for (t in 1:20) {
		#'   x_t = data.frame(x1 = rnorm(1))
		#'   w_t = seq_des$add_one_subject_to_experiment_and_assign(x_t)
		#'   seq_des$add_one_subject_response(t, x_t$x1 + 0.5 * w_t + rt(1, df = 3))
		#' }
		#' seq_des_inf = InferenceContinUnivQuantileRegr$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj, tau = 0.5,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "continuous")
			assertNumeric(tau, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			if (!requireNamespace("quantreg", quietly = TRUE)) {
				stop("Package 'quantreg' is required. Please install it with install.packages(\"quantreg\").")
			}
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
			private$tau = tau
		},

		#' @description
		#' Computes the quantile-regression estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an approximate confidence interval for the treatment effect.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes an approximate two-sided p-value for the treatment effect.
		#' @param delta The null difference to test against. Default is zero.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		tau = NULL,
		fit_warm_keep = NULL,

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},

		build_design_matrix = function(){
			X = cbind(1, private$w)
			colnames(X) = c("(Intercept)", "treatment")
			X
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		set_failed_fit_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
		},

		get_ci_fit_controls = function(){
			ctrl = private$randomization_mc_control
			list(
				warm_start = !is.null(ctrl) && isTRUE(ctrl$fit_warm_start_enable),
				reuse_factorizations = !is.null(ctrl) && isTRUE(ctrl$fit_reuse_factorizations)
			)
		},

		reduce_design_matrix_for_quantile = function(X_full, reuse_factorizations = FALSE){
			if (is.null(dim(X_full))) X_full = matrix(X_full, ncol = 2L)
			if (ncol(X_full) < 2L) return(list(X = NULL, j_treat = NA_integer_))

			if (reuse_factorizations && !is.null(private$fit_warm_keep) && length(private$fit_warm_keep) > 0L &&
				max(private$fit_warm_keep) <= ncol(X_full) && 2L %in% private$fit_warm_keep) {
				X_try = X_full[, private$fit_warm_keep, drop = FALSE]
				j_try = match(2L, private$fit_warm_keep)
				if (!is.na(j_try) && nrow(X_try) > ncol(X_try) && qr(X_try)$rank == ncol(X_try)) {
					return(list(X = X_try, j_treat = j_try))
				}
			}

			reduced = private$reduce_design_matrix_preserving_treatment_fixed_covariates(X_full)
			if (reuse_factorizations && !is.null(reduced$keep) && length(reduced$keep) > 0L && is.finite(reduced$j_treat)) {
				private$fit_warm_keep = reduced$keep
			}
			reduced
		},

		fit_quantile_model = function(X_fit, estimate_only = FALSE){
			if (estimate_only) {
				fit = tryCatch(
					suppressWarnings(getFromNamespace("rq.fit", "quantreg")(x = X_fit, y = private$y, tau = private$tau, method = "br")),
					error = function(e) NULL
				)
				if (is.null(fit) || is.null(fit$coefficients)) return(NULL)
				coef_vec = as.numeric(fit$coefficients)
				if (length(coef_vec) != ncol(X_fit)) return(NULL)
				names(coef_vec) = colnames(X_fit)
				fit$coefficients = coef_vec
				return(fit)
			}

			dat = as.data.frame(X_fit)
			dat$y__ = private$y
			tryCatch(
				suppressWarnings(quantreg::rq(y__ ~ . - 1, tau = private$tau, data = dat)),
				error = function(e) NULL
			)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			fit_controls = private$get_ci_fit_controls()
			X_full = private$build_design_matrix()
			reduced = private$reduce_design_matrix_for_quantile(X_full, reuse_factorizations = fit_controls$reuse_factorizations)
			X_fit = reduced$X
			j_treat = reduced$j_treat

			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			if (is.null(colnames(X_fit)) || length(colnames(X_fit)) != ncol(X_fit)) {
				if (ncol(X_fit) == 1L && isTRUE(j_treat == 1L)) {
					colnames(X_fit) = "treatment"
				} else if (ncol(X_fit) == 1L) {
					colnames(X_fit) = "(Intercept)"
				} else {
					colnames(X_fit) = c("(Intercept)", "treatment", if (ncol(X_fit) > 2L) paste0("x", seq_len(ncol(X_fit) - 2L)) else NULL)[seq_len(ncol(X_fit))]
				}
			}
			coef_names = colnames(X_fit)
			fit = private$fit_quantile_model(X_fit, estimate_only = estimate_only)
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			coef_vec = tryCatch(
				if (estimate_only) as.numeric(fit$coefficients) else as.numeric(stats::coef(fit)),
				error = function(e) NULL
			)
			if (!is.null(coef_vec) && length(coef_vec) == length(coef_names)){
				names(coef_vec) = coef_names
				private$cached_values$full_coefficients = coef_vec
			}

			beta = if (!is.null(coef_vec) && length(coef_vec) >= j_treat) as.numeric(coef_vec[j_treat]) else NA_real_
			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			if (estimate_only) return(invisible(NULL))

			se = .extract_se_from_rq_fit(fit, "treatment")
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
		}
	)
)
