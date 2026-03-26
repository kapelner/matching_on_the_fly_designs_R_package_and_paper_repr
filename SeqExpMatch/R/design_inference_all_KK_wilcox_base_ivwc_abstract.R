# Abstract base class for KK Wilcoxon-based compound inference
#
# @description
# Shared base for all KK Wilcoxon inference classes. Overrides the per-permutation
# statistic used in randomization tests with standardized Wilcoxon W statistics
# (O(n log n), conf.int = FALSE), avoiding the O(n^2) Walsh-average computation
# required by the full Hodges-Lehmann estimate.
#
# @keywords internal
DesignInferenceAbstractKKWilcoxBaseIVWC = R6::R6Class("DesignInferenceAbstractKKWilcoxBaseIVWC",
	inherit = DesignInferenceKKPassThrough,
	public = list(

		# @description
		# Override to avoid O(n^2) per-resample HL computation during the bootstrap warm-start
		# inside compute_confidence_interval_rand. The asymptotic MLE CI is a perfectly
		# adequate starting bound for the bisection and is computed in O(1).
		# @param alpha					The confidence level. Default is 0.05.
		# @param ... 					Additional arguments passed to super.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, ...){
			self$compute_asymp_confidence_interval(alpha)
		}
	),
	private = list(

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			# Optimization: w_mat and m_mat are already pre-computed matrices
			w_mat = permutations$w_mat
			m_mat = permutations$m_mat
			nsim = ncol(w_mat)

			# Reconstruct m_mat if NULL (KK permutations return fixed matching as NULL to save memory)
			if (is.null(m_mat)) {
				m_vec = as.integer(private$des_obj_priv_int$m)
				m_mat = matrix(rep(m_vec, nsim), nrow = length(m_vec), ncol = nsim)
				is_fixed_matching = TRUE
			} else {
				# Check if all matchings are identical
				is_fixed_matching = TRUE
				if (nsim > 1) {
					first_match = m_mat[, 1]
					for (j in 2:min(nsim, 5)) {
						if (!all(m_mat[, j] == first_match)) {
							is_fixed_matching = FALSE
							break
						}
					}
				}
			}

			y_sim = as.numeric(y)
			
			# Map transform_responses to transform_code
			t_code = 0L # none
			if (transform_responses == "log") {
				t_code = 1L
			} else if (transform_responses == "logit") {
				t_code = 2L
			} else if (transform_responses == "log1p") {
				t_code = 3L
			}
			
			res = compute_kk_wilcox_distr_parallel_cpp(
				y_sim,
				w_mat,
				m_mat,
				as.numeric(delta),
				t_code,
				is_fixed_matching,
				private$num_cores
			)
			return(res)
		},

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
			if (m_pairs > 0){
				# Optimization: manually compute Wilcoxon rank sum statistic to avoid R overhead
				# This is equivalent to wilcox.test(diffs)$statistic
				abs_diffs = abs(diffs)
				signs = sign(diffs)
				if (!all(signs == 0)){
					# Fast ranking
					rks = rank(abs_diffs, ties.method = "average")
					W_plus = sum(rks[signs > 0]) + 0.5 * sum(rks[signs == 0])
					
					E_W = m_pairs * (m_pairs + 1) / 4
					V_W = m_pairs * (m_pairs + 1) * (2 * m_pairs + 1) / 24
					stat = stat + (W_plus - E_W) / sqrt(V_W)
					n_components = n_components + 1
				}
			}

			# Reservoir: rank-sum W (standardized)
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (nRT > 0 && nRC > 0){
				y_r = KKstats$y_reservoir
				w_r = KKstats$w_reservoir
				# Optimization: use rank() directly
				rks_r = rank(y_r, ties.method = "average")
				W_r = sum(rks_r[w_r == 1])
				
				E_W = nRT * (nRT + nRC + 1) / 2
				V_W = nRT * nRC * (nRT + nRC + 1) / 12
				stat = stat + (W_r - E_W) / sqrt(V_W)
				n_components = n_components + 1
			}

			if (n_components == 0L) NA_real_ else stat
		}
	)
)
