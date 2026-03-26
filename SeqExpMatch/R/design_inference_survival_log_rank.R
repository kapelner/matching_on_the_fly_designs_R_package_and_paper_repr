#' Log-Rank Inference for Survival Data with Censoring
#'
#' @description
#' Non-parametric all-subject inference for survival outcomes supporting right
#' censoring, based on the standard two-sample log-rank test. The treatment effect
#' estimate is the difference in mean martingale residuals between the treatment and
#' control groups under the pooled null hazard. The p-value uses the classic
#' log-rank score statistic with its hypergeometric tie-adjusted variance.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- DesignInferenceSurvivalLogRank$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceSurvivalLogRank = R6::R6Class("DesignInferenceSurvivalLogRank",
	inherit = DesignInference,
	public = list(

		#' @description
		#' Initialize a log-rank inference object for survival data.
		#' @param seq_des_obj A completed \code{SeqDesign} object with a survival response.
		#' @param num_cores The number of CPU cores to use for parallel processing.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			assertResponseType(seq_des_obj$get_response_type(), "survival")
		},

		#' @description
		#' Computes the treatment-effect estimate on the martingale-residual mean-difference scale.
		compute_treatment_estimate = function(){
			private$compute_shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a (1 - alpha)-level confidence interval based on the asymptotic
		#' normality of the martingale-residual mean-difference estimate.
		#' Falls back to bootstrap if the estimated standard error is unavailable.
		#' @param alpha Significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$compute_shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the standard two-sided log-rank p-value for a zero treatment effect.
		#' @param delta Null treatment effect to test against. Only \code{0} is supported.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$compute_shared()

			if (delta != 0){
				stop("Testing non-zero delta is not yet implemented for DesignInferenceSurvivalLogRank.")
			}

			if (!is.finite(private$cached_values$logrank_var) || private$cached_values$logrank_var <= 0){
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}

			chisq_stat = private$cached_values$logrank_score ^ 2 / private$cached_values$logrank_var
			stats::pchisq(chisq_stat, df = 1, lower.tail = FALSE)
		},

		#' @description
		#' Randomization confidence intervals are not supported for this class because
		#' the martingale-residual score scale is not commensurate with the transformed
		#' time-ratio null used by the randomization CI algorithm.
		#' @param alpha Unused.
		#' @param r Unused.
		#' @param pval_epsilon Unused.
		#' @param show_progress Unused.
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE){
			stop("Randomization confidence intervals are not supported for DesignInferenceSurvivalLogRank due to inconsistent estimator units on the log-rank score scale.")
		}
	),

	private = list(

		compute_shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			logrank_stats = tryCatch(
				fast_logrank_stats_cpp(
					time = as.numeric(private$y),
					dead = as.integer(private$dead),
					w = as.integer(private$w)
				),
				error = function(e) NULL
			)

			if (is.null(logrank_stats)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$logrank_score = NA_real_
				private$cached_values$logrank_var = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(logrank_stats$beta_hat)
			private$cached_values$s_beta_hat_T = as.numeric(logrank_stats$se_beta_hat)
			private$cached_values$logrank_score = as.numeric(logrank_stats$score)
			private$cached_values$logrank_var = as.numeric(logrank_stats$var_score)
			private$cached_values$is_z = TRUE
		}
	)
)
