# Abstract mixin: randomisation-based CI for binary (incidence) outcomes
#
# @description
# Sits between \code{SeqDesignInference} and incidence-focused inference
# classes. Implements the Zhang (2026) exact test-inversion CI and p-values.
#
# @keywords internal
SeqDesignInferenceIncidExactZhangAbstract = R6::R6Class("SeqDesignInferenceIncidExactZhangAbstract",
	inherit = SeqDesignInferenceAbstractZhangCombinedBase,
	public = list(

		# @description
		# Returns the estimated treatment effect (log-odds ratio).
		compute_treatment_estimate = function(){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		# @description
		# Returns the MLE-based confidence interval.
		compute_mle_confidence_interval = function(alpha = 0.05){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		# @description
		# Returns the MLE-based p-value.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		compute_confidence_interval_rand = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		compute_two_sided_pval_for_treatment_effect_rand = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		approximate_bootstrap_distribution_beta_hat_T = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		compute_bootstrap_confidence_interval = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		compute_bootstrap_two_sided_pval = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		# @description
		# Computes the exact randomisation-based CI via Zhang's combined method.
		# @param alpha            Significance level; CI covers 1-alpha.
		# @param pval_epsilon     Bisection convergence tolerance.
		# @param combination_method  How to combine the matched-pair and reservoir p-values.
		compute_exact_confidence_interval = function(alpha = 0.05, pval_epsilon = 0.005, combination_method = "Fisher"){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			assertChoice(combination_method, c("Fisher", "Stouffer", "min_p"))
			private$ci_exact_zhang_combined(alpha, pval_epsilon, combination_method)
		},

		# @description
		# Computes the exact randomisation-based p-value via Zhang's combined method.
		# @param delta            Null treatment effect (log-odds ratio) to test.
		# @param combination_method  How to combine the matched-pair and reservoir p-values.
		compute_exact_two_sided_pval_for_treatment_effect = function(delta = 0, combination_method = "Fisher"){
			assertNumeric(delta)
			assertChoice(combination_method, c("Fisher", "Stouffer", "min_p"))
			private$compute_combined_exact_pval(delta, combination_method)
		}
	),

	private = list(

		# -----------------------------------------------------------------------
		# Internal hooks for the Zhang bisection algorithm
		# -----------------------------------------------------------------------

		compute_treatment_estimate_internal = function(){
			private$incid_exact_zhang_treatment_estimate()
		},

		compute_mle_confidence_interval_internal = function(alpha){
			private$incid_exact_zhang_mle_ci(alpha)
		},

		# -----------------------------------------------------------------------
		# Implementation logic
		# -----------------------------------------------------------------------

		incid_exact_zhang_treatment_estimate = function(){
			# For bisection starting points, we use the standard LogRegr estimate
			# if available, or a simple odds ratio.
			if (!is.null(private$cached_values$KKstats)){
				y_r = private$cached_values$KKstats$y_reservoir
				w_r = private$cached_values$KKstats$w_reservoir
			} else {
				y_r = private$y
				w_r = private$w
			}

			n11 = sum(y_r[w_r == 1L])
			n01 = sum(y_r[w_r == 0L])
			n10 = sum(w_r == 1L) - n11
			n00 = sum(w_r == 0L) - n01

			# Sample Odds Ratio (with continuity correction)
			log((n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)))
		},

		incid_exact_zhang_mle_ci = function(alpha){
			# Standard error for log-odds ratio (Woolf's formula)
			if (!is.null(private$cached_values$KKstats)){
				y_r = private$cached_values$KKstats$y_reservoir
				w_r = private$cached_values$KKstats$w_reservoir
			} else {
				y_r = private$y
				w_r = private$w
			}

			n11 = sum(y_r[w_r == 1L])
			n01 = sum(y_r[w_r == 0L])
			n10 = sum(w_r == 1L) - n11
			n00 = sum(w_r == 0L) - n01

			est = log((n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)))
			se  = sqrt(1/(n11+0.5) + 1/(n10+0.5) + 1/(n01+0.5) + 1/(n00+0.5))

			z = qnorm(1 - alpha/2)
			c(est - z * se, est + z * se)
		},

		compute_combined_exact_pval = function(delta_0, combination_method = "Fisher"){
			# Design validation: only CRD and KK supported by this method
			is_crd = is(private$seq_des_obj, "SeqDesignCRD")
			is_kk  = is(private$seq_des_obj, "SeqDesignKK14")
			if (!is_crd && !is_kk){
				stop("Zhang incidence inference is only supported for CRD (SeqDesignCRD) and KK (SeqDesignKK14 or subclass) designs.")
			}

			if (!is.null(private$cached_values$KKstats)){
				m   = private$cached_values$KKstats$m
				nRT = private$cached_values$KKstats$nRT
				nRC = private$cached_values$KKstats$nRC
			} else {
				m   = 0L
				nRT = sum(private$w == 1L, na.rm = TRUE)
				nRC = sum(private$w == 0L, na.rm = TRUE)
			}

			p_M = if (m > 0)              private$compute_rand_pval_matched_pairs(delta_0) else NA_real_
			p_R = if (nRT > 0 && nRC > 0) private$compute_rand_pval_reservoir(delta_0)    else NA_real_

			private$combine_rand_pvals(p_M, p_R, m, nRT, nRC, combination_method)
		},

		# Default: no matched-pair test (CRD designs have m = 0 always;
		# KK-specific exact classes override this method when needed).
		compute_rand_pval_matched_pairs = function(delta_0) NA_real_,

		# Exact Fisher test under H0: OR_reservoir = exp(delta_0).
		compute_rand_pval_reservoir = function(delta_0){
			if (!is.null(private$cached_values$KKstats)){
				y_r = private$cached_values$KKstats$y_reservoir
				w_r = private$cached_values$KKstats$w_reservoir
				nRT = private$cached_values$KKstats$nRT
				nRC = private$cached_values$KKstats$nRC
			} else {
				y_r = private$y
				w_r = private$w
				nRT = sum(w_r == 1L, na.rm = TRUE)
				nRC = sum(w_r == 0L, na.rm = TRUE)
			}
			if (nRT == 0L || nRC == 0L) return(NA_real_)

			n11 = sum(y_r[w_r == 1L])              # successes in treatment
			n01 = sum(y_r[w_r == 0L])              # successes in control
			n10 = as.integer(nRT) - n11            # failures  in treatment
			n00 = as.integer(nRC) - n01            # failures  in control

			if (n11 + n01 == 0L || n10 + n00 == 0L) return(NA_real_)

			# 2x2 table: rows = Y (1/0), cols = arm (T/C)
			# OR = (n11 * n00) / (n10 * n01) — standard treatment-vs-control OR
			stats::fisher.test(
				matrix(c(n11, n10, n01, n00), 2L, 2L),
				or          = exp(delta_0),
				alternative = "two.sided"
			)$p.value
		}
	)
)
