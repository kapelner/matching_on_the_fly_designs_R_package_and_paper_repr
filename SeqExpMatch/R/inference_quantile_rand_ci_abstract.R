# Abstract mixin: Zhang combined randomisation CI for quantile regression
#
# @description
# Analogous to \code{SeqDesignInferenceAbstractIncidRandCI} but for quantile
# regression estimators.  Provides \code{compute_confidence_interval_rand()}
# via Zhang's combined test-inversion method for both CRD (\eqn{m = 0},
# all subjects in the reservoir) and KK matching-on-the-fly designs
# (\eqn{m > 0}).
#
# @keywords internal
SeqDesignInferenceAbstractQuantileRandCI = R6::R6Class("SeqDesignInferenceAbstractQuantileRandCI",
	inherit = SeqDesignInferenceAbstractZhangCombinedBase,
	public = list(

		# @description
		# Computes a randomization-based confidence interval via Zhang's combined test.
		#
		# @param alpha					The confidence level is 1 - \code{alpha}.
		# @param nsim_exact_test		Number of random sign-flips / permutations.
		# @param pval_epsilon			Bisection convergence tolerance.
		# @param show_progress			Ignored.
		# @return 	A length-2 numeric vector giving the lower and upper CI boundary.
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 499, pval_epsilon = 0.005, show_progress = TRUE){
			if (!is.null(private[["custom_randomization_statistic_function"]])){
				stop("Custom randomization statistic functions are not supported for the Zhang combined CI method used by ", class(self)[1], ". The method uses its own fixed QR-based test statistics.")
			}
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			assertCount(nsim_exact_test, positive = TRUE)
			private$nsim_rand = as.integer(nsim_exact_test)
			private$ci_exact_zhang_combined(alpha, pval_epsilon)
		}
	),

	private = list(
		nsim_rand = 499L
	)
)
