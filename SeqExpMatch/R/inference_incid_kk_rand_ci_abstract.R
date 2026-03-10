#' Extension of SeqDesignInferenceAbstractIncidRandCI for KK matched-pair designs
#'
#' @description
#' Adds the exact McNemar-style matched-pair p-value needed by the
#' \code{"zhang_combined"} randomisation CI when \eqn{m > 0}.
#'
#' @keywords internal
SeqDesignInferenceAbstractIncidKKRandCI = R6::R6Class("SeqDesignInferenceAbstractIncidKKRandCI",
	inherit = SeqDesignInferenceAbstractIncidRandCI,
	private = list(

		# Exact McNemar test under H0: OR_pairs = exp(delta_0).
		# Given k = d_plus + d_minus discordant pairs,
		# d_plus | k ~ Binomial(k, expit(delta_0)) under H0.
		# Concordant pairs contribute no information and are ignored.
		compute_rand_pval_matched_pairs = function(delta_0){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats) || KKstats$m == 0L) return(NA_real_)

			yTs = KKstats$yTs_matched
			yCs = KKstats$yCs_matched
			d_plus  = sum(yTs == 1L & yCs == 0L)   # (T=1, C=0) discordant
			d_minus = sum(yTs == 0L & yCs == 1L)   # (T=0, C=1) discordant
			k = d_plus + d_minus
			if (k == 0L) return(NA_real_)

			p0 = exp(delta_0) / (1 + exp(delta_0))
			stats::binom.test(d_plus, k, p = p0, alternative = "two.sided")$p.value
		}
	)
)
