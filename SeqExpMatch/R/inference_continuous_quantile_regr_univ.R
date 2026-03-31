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
#' \code{Suggests} and is not installed automatically with \pkg{SeqExpMatch}.
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
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "continuous")
		#' for (t in 1:20) {
		#'   x_t = data.frame(x1 = rnorm(1))
		#'   w_t = seq_des$add_subject_to_experiment_and_assign(x_t)
		#'   seq_des$add_subject_response(t, x_t$x1 + 0.5 * w_t + rt(1, df = 3))
		#' }
		#' seq_des_inf = InferenceContinUnivQuantileRegr$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj, tau = 0.5, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "continuous")
			assertNumeric(tau, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			if (!requireNamespace("quantreg", quietly = TRUE)) {
				stop("Package 'quantreg' is required. Please install it with install.packages(\"quantreg\").")
			}
			super$initialize(des_obj, num_cores, verbose)
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

		build_design_matrix = function(){
			cbind(1, private$w)
		},

		set_failed_fit_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}

			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat

			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			coef_names = c("(Intercept)", "treatment", if (ncol(X_fit) > 2L) paste0("x", seq_len(ncol(X_fit) - 2L)) else NULL)
			colnames(X_fit) = coef_names
			dat = as.data.frame(X_fit)
			dat$y__ = private$y

			fit = tryCatch(
				suppressWarnings(quantreg::rq(y__ ~ . - 1, tau = private$tau, data = dat)),
				error = function(e) NULL
			)
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			beta = tryCatch(as.numeric(coef(fit)[["treatment"]]), error = function(e) NA_real_)
			se = .extract_se_from_rq_fit(fit, "treatment")

			full_coef = tryCatch(as.numeric(coef(fit)), error = function(e) NULL)
			if (!is.null(full_coef) && length(full_coef) == length(coef_names)){
				names(full_coef) = coef_names
				private$cached_values$full_coefficients = full_coef
			}

			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
		}
	)
)
