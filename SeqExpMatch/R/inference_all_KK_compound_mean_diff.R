#' Inference based on Maximum Likelihood for KK designs
#'
#' @description
#' Inference for mean difference
#'
#'
#' @export
#' @examples
#' \dontrun{
#' seq_des = DesignSeqOneByOneKK14$new(n = 6, response_type = "continuous")
#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
#'
#' seq_des_inf = InferenceAllKKCompoundMeanDiff$
#'   new(seq_des)
#' seq_des_inf$compute_treatment_estimate()
#' seq_des_inf$compute_asymp_confidence_interval()
#' seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect()
#' }
InferenceAllKKCompoundMeanDiff = R6::R6Class("InferenceAllKKCompoundMeanDiff",
	inherit = InferenceKKPassThroughCompound,
	public = list(


		#'
		#' @description
		#' Computes the appropriate estimate for compound mean difference across pairs and reservoir
		#'
		#' @return	The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		compute_treatment_estimate = function(){
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
				private$compute_reservoir_and_match_statistics()
			}
			private$cached_values$beta_hat_T = 	if (private$cached_values$KKstats$nRT <= 1 || private$cached_values$KKstats$nRC <= 1){
													private$cached_values$KKstats$d_bar
												} else if (private$cached_values$KKstats$m == 0){ #sometimes there's no matches
													private$cached_values$KKstats$r_bar
												} else {
													private$cached_values$KKstats$w_star * private$cached_values$KKstats$d_bar + (1 - private$cached_values$KKstats$w_star) * private$cached_values$KKstats$r_bar #proper weighting
												}
			private$cached_values$beta_hat_T
		},

		#' @description


		#' Computes a 1-alpha level frequentist confidence interval
		#'
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal (except in
		#' the case
		#' of estimat_type "median difference" where a nonparametric bootstrap confidence interval
		#' (see the \code{controlTest::quantileControlTest} method)
		#' is employed. Hence these confidence intervals are asymptotically valid and thus approximate
		#' for any sample size.
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#'
		#' @return	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		#'
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)

			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description


		#' Computes a 2-sided p-value
		#'
		#' @param delta   The null difference to test against. For any treatment effect at all this is
		#'   set to zero (the default).
		#'
		#' @return	The approximate frequentist p-value
		#'
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#' @param	r		The number of randomization vectors. The default is 501.
		#' @param	pval_epsilon			The bisection algorithm tolerance. The default is 0.005.
		#' @param	show_progress			Show a text progress indicator.
		#' @return	A 1 - alpha sized frequentist confidence interval
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE){
			if (private$des_obj_priv_int$response_type %in% c("proportion", "count", "survival")) {
				stop("Randomization confidence intervals are not supported for InferenceAllKKCompoundMeanDiff with proportion, count, or survival response types due to inconsistent estimator units on the transformed scale.")
			}
			super$compute_confidence_interval_rand(alpha = alpha, r = r, pval_epsilon = pval_epsilon, show_progress = show_progress)
		}
	),

	private = list(
		compute_fast_bootstrap_distr = function(B, i_reservoir, n_reservoir, m, y, w, m_vec) {
			# Only safe for simple additive/linear combinations right now.
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)

			n = length(y)
			y_mat = matrix(0.0, nrow = n, ncol = B)
			w_mat = matrix(0L, nrow = n, ncol = B)
			m_mat = matrix(0L, nrow = n, ncol = B)

			for (b in 1:B) {
				# Resample reservoir with replacement
				i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)

				# For matched pairs, sample which pairs to include (with replacement)
				if (m > 0) {
					pairs_to_include = sample(1:m, m, replace = TRUE)
					i_matched_b = integer(0)
					m_vec_b_matched = integer(0)
					for (new_pair_id in 1:m) {
						original_pair_id = pairs_to_include[new_pair_id]
						pair_indices = which(m_vec == original_pair_id)
						i_matched_b = c(i_matched_b, pair_indices)
						m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
					}
				} else {
					i_matched_b = integer(0)
					m_vec_b_matched = integer(0)
				}

				# Combine reservoir and matched indices
				i_b = c(i_reservoir_b, i_matched_b)

				y_mat[, b] = y[i_b]
				w_mat[, b] = w[i_b]
				m_mat[, b] = c(rep(0L, n_reservoir), m_vec_b_matched)
			}

			res = compute_kk_compound_bootstrap_parallel_cpp(
				y_mat,
				w_mat,
				m_mat,
				private$num_cores
			)

			return(res)
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses) {
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			if (delta != 0) return(NULL)

			n = length(y)
			w_mat = permutations$w_mat
			m_mat = permutations$m_mat
			if (is.null(m_mat)) {
				m_mat = matrix(0L, nrow = n, ncol = ncol(w_mat))
			} else {
				m_mat[is.na(m_mat)] = 0L
			}

			res = compute_kk_compound_distr_parallel_cpp(
				as.numeric(y),
				w_mat,
				m_mat,
				private$num_cores
			)
			return(res)
		},

		shared = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				self$compute_treatment_estimate()
			}

			ssqD = private$cached_values$KKstats$ssqD_bar
			ssqR = private$cached_values$KKstats$ssqR

			private$cached_values$s_beta_hat_T =
				if (private$cached_values$KKstats$nRT <= 1 || private$cached_values$KKstats$nRC <= 1){
					# Only matched pairs are usable; fall back to ssqR if ssqD is degenerate
					if (is.finite(ssqD) && ssqD > 0) sqrt(ssqD) else if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
				} else if (private$cached_values$KKstats$m == 0){
					# No matched pairs
					if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
				} else {
					# Combined: require both components to be positive and finite.
					# If one collapses (e.g. tiny reservoir with near-identical responses),
					# fall back to the other rather than pulling the combined SE to zero.
					if (!is.finite(ssqD) || ssqD <= 0) {
						if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
					} else if (!is.finite(ssqR) || ssqR <= 0) {
						sqrt(ssqD)
					} else {
						sqrt(ssqR * ssqD / (ssqR + ssqD))
					}
				}
		}
	)
)
