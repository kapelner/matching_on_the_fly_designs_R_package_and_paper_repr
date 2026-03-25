# Abstract mixin: randomisation-based CI for binary (incidence) outcomes
#
# @description
# Sits between \code{DesignInference} and incidence-focused inference
# classes. Implements the Zhang (2026) exact test-inversion CI and p-values.
#
# @keywords internal
DesignInferenceIncidExactZhangAbstract = R6::R6Class("DesignInferenceIncidExactZhangAbstract",
	inherit = DesignInferenceAbstractZhangCombinedBase,
	public = list(

		# @description
		# Returns the estimated treatment effect (log-odds ratio).
		compute_treatment_estimate = function(){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		# @description
		# Returns the asymptotic confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect and compute_exact_confidence_interval only. This method is not implemented.")
		},

		# @description
		# Returns the asymptotic p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
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

		compute_asymp_confidence_interval_internal = function(alpha){
			private$incid_exact_zhang_mle_ci(alpha)
		},

		# -----------------------------------------------------------------------
		# Implementation logic
		# -----------------------------------------------------------------------

		incid_exact_zhang_treatment_estimate = function(){
			# For bisection starting points, we use a simple odds ratio on the
			# reservoir counts, which are cached once and reused across CI calls.
			exact_stats = private$get_exact_zhang_stats()
			n11 = exact_stats$n11
			n10 = exact_stats$n10
			n01 = exact_stats$n01
			n00 = exact_stats$n00

			# Sample Odds Ratio (with continuity correction)
			log((n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)))
		},

		incid_exact_zhang_mle_ci = function(alpha){
			# Standard error for log-odds ratio (Woolf's formula)
			exact_stats = private$get_exact_zhang_stats()
			n11 = exact_stats$n11
			n10 = exact_stats$n10
			n01 = exact_stats$n01
			n00 = exact_stats$n00

			est = log((n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)))
			se  = sqrt(1/(n11+0.5) + 1/(n10+0.5) + 1/(n01+0.5) + 1/(n00+0.5))

			z = qnorm(1 - alpha/2)
			c(est - z * se, est + z * se)
		},

		compute_combined_exact_pval = function(delta_0, combination_method = "Fisher"){
			# Design validation: only Bernoulli and KK supported by this method
			is_bernoulli = is(private$seq_des_obj, "SeqDesignBernoulli")
			is_kk  = is(private$seq_des_obj, "SeqDesignKK14")
			if (!is_bernoulli && !is_kk){
				stop("Zhang incidence inference is only supported for Bernoulli (SeqDesignBernoulli) and KK (SeqDesignKK14 or subclass) designs.")
			}

			exact_stats = private$get_exact_zhang_stats()
			p_M = if (exact_stats$m > 0)              private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
			p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) private$compute_exact_pval_reservoir(delta_0)    else NA_real_

			private$combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
		},

		# Default: no matched-pair test (Bernoulli designs have m = 0 always;
		# KK-specific exact classes override this method when needed).
		compute_exact_pval_matched_pairs = function(delta_0) NA_real_,

		# Exact Fisher test under H0: OR_reservoir = exp(delta_0).
		# For Bernoulli designs (KKstats = NULL) all subjects are treated as the reservoir.
		compute_exact_pval_reservoir = function(delta_0){
			exact_stats = private$get_exact_zhang_stats()
			if (exact_stats$nRT == 0L || exact_stats$nRC == 0L) return(NA_real_)
			if (exact_stats$n11 + exact_stats$n01 == 0L || exact_stats$n10 + exact_stats$n00 == 0L) return(NA_real_)

			zhang_exact_fisher_pval_cpp(
				exact_stats$n11,
				exact_stats$n10,
				exact_stats$n01,
				exact_stats$n00,
				delta_0
			)
		},

		get_exact_zhang_stats = function(){
			if (!is.null(private$cached_values$incid_exact_zhang_stats)){
				return(private$cached_values$incid_exact_zhang_stats)
			}

			if (!is.null(private$cached_values$KKstats)){
				KKstats = private$cached_values$KKstats
				exact_stats = list(
					m       = as.integer(KKstats$m),
					nRT     = as.integer(KKstats$nRT),
					nRC     = as.integer(KKstats$nRC),
					d_plus  = as.integer(KKstats$d_plus),
					d_minus = as.integer(KKstats$d_minus),
					n11     = as.integer(KKstats$n11),
					n10     = as.integer(KKstats$n10),
					n01     = as.integer(KKstats$n01),
					n00     = as.integer(KKstats$n00)
				)
			} else {
				y = private$y
				w = private$w
				nRT = sum(w == 1L, na.rm = TRUE)
				nRC = sum(w == 0L, na.rm = TRUE)
				n11 = sum(y[w == 1L])
				n01 = sum(y[w == 0L])
				exact_stats = list(
					m       = 0L,
					nRT     = as.integer(nRT),
					nRC     = as.integer(nRC),
					d_plus  = 0L,
					d_minus = 0L,
					n11     = as.integer(n11),
					n10     = as.integer(nRT - n11),
					n01     = as.integer(n01),
					n00     = as.integer(nRC - n01)
				)
			}

			private$cached_values$incid_exact_zhang_stats = exact_stats
			exact_stats
		}
	)
)
