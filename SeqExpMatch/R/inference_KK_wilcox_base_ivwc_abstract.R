# Abstract base class for KK Wilcoxon-based compound inference
#
# @description
# Shared base for all KK Wilcoxon inference classes. Overrides the per-permutation
# statistic used in randomization tests with standardized Wilcoxon W statistics
# (O(n log n), conf.int = FALSE), avoiding the O(n^2) Walsh-average computation
# required by the full Hodges-Lehmann estimate.
#
# @keywords internal
SeqDesignInferenceAbstractKKWilcoxBaseIVWC = R6::R6Class("SeqDesignInferenceAbstractKKWilcoxBaseIVWC",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		# @description
		# Override to avoid O(n^2) per-resample HL computation during the bootstrap warm-start
		# inside compute_confidence_interval_rand. The asymptotic MLE CI is a perfectly
		# adequate starting bound for the bisection and is computed in O(1).
		# @param alpha					The confidence level. Default is 0.05.
		# @param ... 					Additional arguments passed to super.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, ...){
			self$compute_mle_confidence_interval(alpha)
		}
	),
	private = list(

		# Override the per-permutation statistic to avoid the O(n^2) conf.int = TRUE cost.
		# Uses standardized Wilcoxon W statistics (conf.int = FALSE, O(n log n)) as the
		# rank test statistic; monotone with the HL / rank-regression estimate under the null.
		compute_treatment_estimate_during_randomization_inference = function(){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}

			stat = 0
			n_components = 0

			# Matched pairs: signed-rank W (standardized)
			diffs = KKstats$y_matched_diffs
			m_pairs = length(diffs)
			if (m_pairs > 0 && !all(diffs == 0)){
				W_m = tryCatch(wilcox.test(diffs, conf.int = FALSE)$statistic, error = function(e) NA_real_)
				if (is.finite(W_m)){
					E_W = m_pairs * (m_pairs + 1) / 4
					V_W = m_pairs * (m_pairs + 1) * (2 * m_pairs + 1) / 24
					stat = stat + (W_m - E_W) / sqrt(V_W)
					n_components = n_components + 1
				}
			}

			# Reservoir: rank-sum W (standardized)
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (nRT > 0 && nRC > 0){
				y_r = KKstats$y_reservoir
				w_r = KKstats$w_reservoir
				yT = y_r[w_r == 1]
				yC = y_r[w_r == 0]
				W_r = tryCatch(wilcox.test(yT, yC, conf.int = FALSE)$statistic, error = function(e) NA_real_)
				if (is.finite(W_r)){
					E_W = nRT * nRC / 2
					V_W = nRT * nRC * (nRT + nRC + 1) / 12
					stat = stat + (W_r - E_W) / sqrt(V_W)
					n_components = n_components + 1
				}
			}

			if (n_components == 0L) NA_real_ else stat
		}
	)
)
