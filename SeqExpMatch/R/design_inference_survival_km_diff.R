#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
DesignInferenceSurvivalKMDiff = R6::R6Class("DesignInferenceSurvivalKMDiff",
	inherit = DesignInferenceAsymp,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
			assertResponseType(des_obj$get_response_type(), "survival")
		},

		#' @description
		#' Computes the appropriate estimate for mean difference
		#'
		#' @return	The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @examples
		#' seq_des = SeqDesignBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = DesignInferenceSurvivalKMDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#'
		compute_treatment_estimate = function(){
			get_survival_stat_diff(
				private$y,
				private$dead,
				private$w,
				"median"
			)
		},

		#' @description
		#' Computes a (1 - alpha)-level confidence interval for the difference in Kaplan-Meier
		#' median survival times (treatment minus control).
		#'
		#' The Brookmeyer-Crowley confidence interval is obtained for each group's median
		#' separately via \code{survival::survfit} (using a log-log transformation of the
		#' survival function by default). The per-group SE is back-calculated from the CI
		#' half-width as \eqn{\hat\sigma_i = (\text{upper}_i - \text{lower}_i) / (2 z_{\alpha/2})}.
		#' The two groups are independent by design, so the SE of the difference is
		#' \eqn{\sqrt{\hat\sigma_T^2 + \hat\sigma_C^2}}, and the CI is
		#' \eqn{(\hat{m}_T - \hat{m}_C) \pm z_{\alpha/2} \cdot \sqrt{\hat\sigma_T^2 +
		#' \hat\sigma_C^2}}.
		#'
		#' Falls back to \code{compute_bootstrap_confidence_interval} when either group's
		#' median is not estimable (i.e., the Kaplan-Meier curve does not reach 0.5) or
		#' when the Brookmeyer-Crowley CI bounds are \code{NA}.
		#'
		#' @param alpha           The significance level; the confidence level is 1 - \code{alpha}.
		#'   Default is 0.05.
		#'
		#' @return	A numeric vector of length 2 giving the (lower, upper) confidence bounds
		#' 			for the difference in median survival times, on the original time scale.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = DesignInferenceSurvivalKMDiff$new(seq_des)
		#' seq_des_inf$compute_asymp_confidence_interval()
		#' }
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			y    = private$y
			dead = private$dead
			w    = private$w
			fit_T = tryCatch(survival::survfit(survival::Surv(y[w == 1], dead[w == 1]) ~ 1, conf.int = 1 - alpha), error = function(e) NULL)
			fit_C = tryCatch(survival::survfit(survival::Surv(y[w == 0], dead[w == 0]) ~ 1, conf.int = 1 - alpha), error = function(e) NULL)
			if (is.null(fit_T) || is.null(fit_C)){
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			q_T = quantile(fit_T, 0.5)
			q_C = quantile(fit_C, 0.5)
			med_T = as.numeric(q_T$quantile)
			lo_T  = as.numeric(q_T$lower)
			hi_T  = as.numeric(q_T$upper)
			med_C = as.numeric(q_C$quantile)
			lo_C  = as.numeric(q_C$lower)
			hi_C  = as.numeric(q_C$upper)
			# Fall back to bootstrap if either median or its CI is not estimable
			if (!is.finite(med_T) || !is.finite(med_C) || !is.finite(lo_T) || !is.finite(hi_T) || !is.finite(lo_C) || !is.finite(hi_C)){
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			# Back-calculate SE for each median from the Brookmeyer-Crowley CI,
			# then combine under independence for the difference
			z        = stats::qnorm(1 - alpha / 2)
			se_diff  = sqrt(((hi_T - lo_T) / (2 * z))^2 + ((hi_C - lo_C) / (2 * z))^2)
			diff     = med_T - med_C
			c(diff - z * se_diff, diff + z * se_diff)
		},

		#' @description
		#' Computes a 2-sided p-value via the log rank test
		#'
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		#'
		#' @return	The approximate frequentist p-value
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignBernoulli$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(
		#'   ys = c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43),
		#'   deads = c(1L, 0L, 1L, 1L, 0L, 1L)
		#' )
		#'
		#' seq_des_inf = DesignInferenceSurvivalKMDiff$new(seq_des)
		#' seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect()
		#' }
		#'
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)

			if (delta == 0){
				survival_obj = survival::Surv(private$y, private$dead)
				surv_diff = survival::survdiff(survival_obj ~ private$w)
				surv_diff$pvalue
			} else {
				stop("TO-DO")
			}
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#' @param	r		The number of randomization vectors. The default is 501.
		#' @param	pval_epsilon			The bisection algorithm tolerance. The default is 0.005.
		#' @param	show_progress		Show a text progress indicator.
		#' @return	A 1 - alpha sized frequentist confidence interval
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE){
			stop("Randomization confidence intervals are not supported for DesignInferenceSurvivalKMDiff due to inconsistent estimator units on the transformed scale.")
		}
	),

	private = list()
)
