#' Non-parametric Wilcoxon-based Compound Inference for KK Designs
#'
#' @description
#' Fits a non-parametric compound estimator for KK matching-on-the-fly designs.
#' For matched pairs, it uses the Wilcoxon Signed-Rank Hodges-Lehmann estimate.
#' For reservoir subjects, it uses the Wilcoxon Rank-Sum (Mann-Whitney U) Hodges-Lehmann
#' estimate. The two estimates are combined via a variance-weighted linear combination.
#' This method is robust to outliers and does not assume a specific parametric
#' distribution for the response.
#'
#' @inherit InferenceRand methods
#' @inherit InferenceBoot methods
#' @inherit InferenceAsymp methods
#' @inherit InferenceRandCI methods
#' @export
InferenceAllKKWilcoxIVWC = R6::R6Class("InferenceAllKKWilcoxIVWC",
	inherit = InferenceAbstractKKWilcoxBaseIVWC,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param	des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- DesignSeqOneByOneKK14$new(n = nrow(x_dat), response_type = "continuous", verbose =
		#' FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- InferenceAllKKWilcoxIVWC$new(seq_des, verbose = FALSE)
		#' infer
		#'
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			res_type = des_obj$get_response_type()
			if (res_type == "incidence"){
				stop("Rank-based compound inference is not recommended for incidence data; clogit or compound mean difference estimators are preferred.")
			}
			assertResponseType(res_type, c("continuous", "count", "proportion", "survival", "ordinal"))
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
			if (private$any_censoring){
				stop(class(self)[1], " does not currently support censored survival data. Use restricted mean or Cox-based methods instead.")
			}
		},

		#' @description
		#' Returns the estimated treatment effect (Hodges-Lehmann median shift).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the non-parametric confidence interval.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			# Even though estimates are non-parametric, the combined estimator is asymptotically normal
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the non-parametric p-value.
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("Testing non-zero delta is not yet implemented for the combined rank estimator.")
			}
		},

		#' @description
		#' Computes the bootstrap confidence interval.
		#' @param alpha The confidence level. Default is 0.05.
		#' @param ... Additional arguments passed to the superclass method.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, ...){
			super$compute_bootstrap_confidence_interval(alpha = alpha, ...)
		}
	),

	private = list(

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			# Recompute KKstats if cache was cleared (e.g., after y transformation for rand CI)
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: Wilcoxon Signed-Rank HL Estimate ---
			if (m > 0){
				private$rank_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: Wilcoxon Rank-Sum HL Estimate ---
			if (nRT > 0 && nRC > 0){
				private$rank_for_reservoir()
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0

			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Rank compound estimator: could not compute a finite standard error.")
			}
		},

		rank_for_matched_pairs = function(){
			diffs = private$cached_values$KKstats$y_matched_diffs
			# signed-rank test requires at least some non-zero differences
			if (all(diffs == 0)) return(invisible(NULL))

			mod = tryCatch({
				stats::wilcox.test(diffs, conf.int = TRUE)
			}, error = function(e) NULL)

			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(mod$estimate)
			# Back-calculate SE from 95% CI width; guard against NULL/empty conf.int
			ci = mod$conf.int
			se = if (length(ci) == 2L) (ci[2] - ci[1]) / (2 * 1.96) else NA_real_

			private$cached_values$beta_T_matched     = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (length(se) == 1L && is.finite(se) && se > 0) se^2 else NA_real_
		},

		rank_for_reservoir = function(){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			yT = y_r[w_r == 1]
			yC = y_r[w_r == 0]

			mod = tryCatch({
				stats::wilcox.test(yT, yC, conf.int = TRUE)
			}, error = function(e) NULL)

			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(mod$estimate)
			# Back-calculate SE from 95% CI width; guard against NULL/empty conf.int
			ci = mod$conf.int
			se = if (length(ci) == 2L) (ci[2] - ci[1]) / (2 * 1.96) else NA_real_

			private$cached_values$beta_T_reservoir     = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (length(se) == 1L && is.finite(se) && se > 0) se^2 else NA_real_
		},

		compute_fast_bootstrap_distr = function(B, i_reservoir, n_reservoir, m, y, w, m_vec){
			# Generate bootstrap indices for KK design in R first
			indices_mat = matrix(NA_integer_, nrow = length(y), ncol = B)
			m_mat = matrix(NA_integer_, nrow = length(y), ncol = B)
			w_mat = matrix(NA_integer_, nrow = length(y), ncol = B)

			for (b in 1:B) {
				i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
				w_b_res = w[i_reservoir_b]
				
				i_matched_b = integer(0)
				m_vec_b_matched = integer(0)
				w_b_matched = integer(0)
				
				if (m > 0) {
					pairs_to_include = sample(1:m, m, replace = TRUE)
					for (new_pair_id in 1:m) {
						original_pair_id = pairs_to_include[new_pair_id]
						pair_indices = which(m_vec == original_pair_id)
						i_matched_b = c(i_matched_b, pair_indices)
						m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
						w_b_matched = c(w_b_matched, w[pair_indices])
					}
				}
				
				w_mat[, b] = c(w_b_res, w_b_matched)
				m_mat[, b] = c(rep(0L, n_reservoir), m_vec_b_matched)
				# Note: y is retrieved inside C++ using indices_mat
				indices_mat[, b] = c(i_reservoir_b, i_matched_b)
			}

			compute_wilcox_kk_ivwc_bootstrap_parallel_cpp(
				as.numeric(y),
				as.integer(w),
				as.integer(m_vec),
				indices_mat, # Back to 1-based, let C++ handle it if needed or check impl
				m_mat,
				private$num_cores
			)
		}
	)
)
