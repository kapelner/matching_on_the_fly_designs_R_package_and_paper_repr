#' Gehan-Wilcoxon (Peto-Prentice) Inference for Survival Data with Censoring
#'
#' @description
#' Non-parametric inference for survival outcomes supporting censored data, using
#' the Peto-Prentice modification of the Gehan-Wilcoxon test. The treatment effect
#' estimate is the mean difference in Peto-Prentice weighted martingale residuals
#' between the treatment and control groups. Specifically, for each subject the
#' weighted residual is \eqn{M_i^w = \hat{S}(t_i^-) \cdot M_i}, where
#' \eqn{M_i = \delta_i - \hat\Lambda_0(t_i)} is the martingale residual and
#' \eqn{\hat{S}(t_i^-)} is the overall Kaplan-Meier survival estimate just before
#' time \eqn{t_i}. These weights downweight late events, analogously to the
#' Wilcoxon rank-sum test for uncensored data (which also weights early observations
#' more heavily via their larger rank denominator).
#'
#' The p-value uses \code{survival::survdiff(rho = 1)} (Peto-Prentice / Fleming-Harrington
#' p=1, q=0), which is distinct from the log-rank test (\code{rho = 0}) used in
#' \code{DesignInferenceSurvivalKMDiff}.
#'
#' @export
DesignInferenceSurvivalGehanWilcox = R6::R6Class("DesignInferenceSurvivalGehanWilcox",
	inherit = DesignInference,
	public = list(

		#' @description
		#' Initialize a Gehan-Wilcoxon (Peto-Prentice) inference object for survival data.
		#' @param seq_des_obj     A SeqDesign object whose entire n subjects are assigned and response
		#'   y is recorded within.
		#' @param	num_cores		The number of CPU cores to use for parallel processing. Default is 1.
		#' @param	verbose		Whether to print progress messages. Default is FALSE.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			assertResponseType(seq_des_obj$get_response_type(), "survival")
		},

		#' @description
		#' Returns the mean difference in Peto-Prentice weighted martingale residuals
		#' between the treatment and control groups. Positive values indicate that treatment
		#' subjects experienced fewer early events than expected.
		#'
		#' @return	A numeric scalar (the Peto-Prentice weighted score treatment effect estimate).
		#'
		#' @examples
		#' seq_des = SeqDesignBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2:10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = DesignInferenceSurvivalGehanWilcox$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		compute_treatment_estimate = function(){
			private$compute_shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a (1 - alpha)-level confidence interval based on the asymptotic normality
		#' of the Peto-Prentice weighted martingale residual mean difference. Falls back to
		#' bootstrap if the SE is unavailable.
		#'
		#' @param	alpha	Significance level. Default is 0.05.
		#'
		#' @return	A numeric vector of length 2: (lower, upper) confidence bounds.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2:10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = DesignInferenceSurvivalGehanWilcox$new(seq_des)
		#' seq_des_inf$compute_asymp_confidence_interval()
		#' }
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$compute_shared()
			if (is.na(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the Peto-Prentice (Gehan-Wilcoxon) two-sided p-value via
		#' \code{survival::survdiff(rho = 1)}, which puts greater weight on early events
		#' relative to the standard log-rank test (\code{rho = 0}).
		#' For delta != 0, not yet implemented.
		#'
		#' @param	delta	Null treatment effect to test against. Default is 0.
		#'
		#' @return	A p-value in [0, 1].
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2:10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2:10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = DesignInferenceSurvivalGehanWilcox$new(seq_des)
		#' seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect()
		#' }
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (delta == 0){
				surv_obj  = survival::Surv(private$y, private$dead)
				surv_diff = survival::survdiff(surv_obj ~ private$w, rho = 1)
				surv_diff$pvalue
			} else {
				stop("Testing non-zero delta is not yet implemented for DesignInferenceSurvivalGehanWilcox.")
			}
		},

		#' @description
		#' Randomization confidence intervals are not supported for this class because
		#' the Peto-Prentice weighted score scale is not commensurate with the time-ratio
		#' null used by the randomization CI bisection algorithm.
		#'
		#' @param	alpha			Unused.
		#' @param	nsim_exact_test	Unused.
		#' @param	pval_epsilon	Unused.
		#' @param	show_progress	Unused.
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.005, show_progress = TRUE){
			stop("Randomization confidence intervals are not supported for DesignInferenceSurvivalGehanWilcox due to inconsistent estimator units on the Peto-Prentice score scale.")
		}
	),

	private = list(

		# Computes the Peto-Prentice weighted martingale residual estimate and SE.
		# Results are cached in private$cached_values.
		compute_shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			y    = private$y
			dead = private$dead
			w    = private$w

			surv_obj = survival::Surv(y, dead)

			# Martingale residuals M_i = delta_i - Lambda_hat_0(t_i) from null Cox model
			cox_null = tryCatch(
				survival::coxph(surv_obj ~ 1),
				error = function(e) NULL
			)
			if (is.null(cox_null)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}
			M = tryCatch(
				as.numeric(residuals(cox_null, type = "martingale")),
				error = function(e) NULL
			)
			if (is.null(M) || length(M) != length(y) || !all(is.finite(M))){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}

			# Peto-Prentice weights: S_hat(t_i^-) from the overall KM estimate.
			# KM is right-continuous, so S(t^-) = km$surv at the largest event time < t.
			# findInterval with left.open=TRUE finds the largest km_time strictly < y[i].
			km_all = tryCatch(
				survival::survfit(surv_obj ~ 1),
				error = function(e) NULL
			)
			if (is.null(km_all)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}
			idx          = findInterval(y, km_all$time, left.open = TRUE)
			peto_weights = c(1.0, km_all$surv)[idx + 1L]

			# Weighted martingale residuals
			M_w  = peto_weights * M
			M_wT = M_w[w == 1]
			M_wC = M_w[w == 0]
			n_T  = length(M_wT)
			n_C  = length(M_wC)

			beta_hat = mean(M_wT) - mean(M_wC)

			# Welch-style SE treating Peto-Prentice weights as fixed
			v_T = if (n_T > 1L) var(M_wT) / n_T else 0
			v_C = if (n_C > 1L) var(M_wC) / n_C else 0
			se  = sqrt(v_T + v_C)

			private$cached_values$beta_hat_T   = beta_hat
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z         = TRUE
		}
	)
)
