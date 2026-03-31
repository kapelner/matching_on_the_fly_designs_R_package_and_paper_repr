#' Univariate Robust Regression Inference for Continuous Responses
#'
#' Fits a robust linear regression via \code{MASS::rlm} for continuous responses
#' using only the treatment indicator as a predictor (intercept + treatment).
#' This provides a Huber/MM-style robustness upgrade over ordinary least squares
#' when outcomes are heavy-tailed or outlier-prone. Inference is based on the
#' coefficient table returned by \code{summary.rlm()}.
#'
#' @details
#' The \code{method} argument is passed to \code{MASS::rlm} and may be either
#' \code{"M"} or \code{"MM"}. For \code{"M"}, the fit uses Huber's psi
#' function. Approximate confidence intervals and p-values use the reported
#' robust standard error with residual degrees of freedom \eqn{n - p}.
#'
#' @export
InferenceContinUnivRobustRegr = R6::R6Class("InferenceContinUnivRobustRegr",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a robust-regression inference object for a completed design
		#' with a continuous response.
		#' @param des_obj             A DesignSeqOneByOne object whose entire n subjects are assigned and
		#'   whose continuous response y is recorded.
		#' @param method                  Robust-regression fitting method for \code{MASS::rlm}; one
		#'   of \code{"M"} or \code{"MM"}. The default is \code{"MM"}.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#' 							and bootstrap resampling.
		#' @param verbose                 A flag indicating whether messages should be displayed to
		#'   the user. Default is \code{FALSE}.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 20, response_type = "continuous")
		#' for (i in 1 : 20){
		#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		#' }
		#' seq_des$add_all_subject_responses(rt(20, df = 3))
		#'
		#' seq_des_inf = InferenceContinUnivRobustRegr$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(des_obj, method = "MM", num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "continuous")
			assertChoice(method, c("M", "MM"))
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			private$rlm_method = method
		},

		#' @description
		#' Computes the robust-regression estimate of the treatment effect.
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
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		rlm_method = NULL,

		build_design_matrix = function(){
			cbind(1, private$w)
		},

		set_failed_fit_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = FALSE
			private$cached_values$df = NA_real_
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))) {
				X_full = matrix(X_full, ncol = 2L)
			}
			if (ncol(X_full) < 2L) {
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) colnames(private$get_X()) else NULL)
			qr_X = qr(X_full)
			keep = qr_X$pivot[seq_len(qr_X$rank)]
			if (!(2L %in% keep)) {
				keep[qr_X$rank] = 2L
			}
			keep = sort(unique(keep))
			X = X_full[, keep, drop = FALSE]

			df = nrow(X) - ncol(X)
			if (df <= 0L) {
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			mod = tryCatch({
				if (identical(private$rlm_method, "M")) {
					MASS::rlm(x = X, y = private$y, method = "M", psi = MASS::psi.huber)
				} else {
					MASS::rlm(x = X, y = private$y, method = private$rlm_method)
				}
			}, error = function(e) NULL)

			if (is.null(mod)) {
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if (is.null(coef_table) || !("treatment" %in% rownames(coef_table))) {
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}
			beta = as.numeric(coef_table["treatment", "Value"])
			se = as.numeric(coef_table["treatment", "Std. Error"])
			private$cached_values$full_coefficients = stats::coef(mod)
			private$cached_values$full_vcov = tryCatch(
				stats::vcov(mod),
				error = function(e) {
					mat = matrix(NA_real_, ncol(X), ncol(X))
					colnames(mat) = rownames(mat) = colnames(X)
					mat
				}
			)

			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = FALSE
			private$cached_values$df = df
		}
	)
)
