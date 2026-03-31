#' Abstract class for Clayton Copula / Standard Weibull Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Clayton copula with Weibull AFT margins for matched pairs.
#' The matched-pair component is fitted by maximizing the pairwise copula likelihood
#' under right censoring. The reservoir component uses standard Weibull AFT regression.
#' The two treatment-effect estimates are combined by inverse-variance weighting.
#'
#' @details
#' This is the package's first copula-based bivariate survival implementation. The matched-pair
#' contribution uses a Clayton survival copula with Weibull AFT marginal models, so the
#' treatment effect is reported on the AFT log-time-ratio scale, matching the existing Weibull
#' survival classes. The copula dependence parameter is estimated jointly with the matched-pair
#' regression coefficients by direct likelihood maximization.
#'
#' @keywords internal
InferenceAbstractKKClaytonCopulaIVWC = R6::R6Class("InferenceAbstractKKClaytonCopulaIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "survival")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
		},

		#' @description
		#' Returns the estimated treatment effect (log-time ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				# For non-zero delta, we can still use the same logic if we have beta_hat_T and s_beta_hat_T
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			}
		},

		#' @description
		#' Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		duplicate = function(verbose = FALSE){
			inf_obj = super$duplicate(verbose = verbose)
			inf_obj
		},

		#' @description
		#' Overridden to use the best design matrix and starting parameters from the initial fit
		#' to speed up randomization-based inference.
		compute_treatment_estimate_during_randomization_inference = function(){
			# Ensure we have the best design and parameters from the original data
			if (is.null(private$best_Xmm_colnames_matched) && is.null(private$best_Xmm_colnames_reservoir)){
				private$shared()
			}

			# If we still don't have enough (e.g., initial fit failed), fall back to standard
			if (is.null(private$best_Xmm_colnames_matched) && is.null(private$best_Xmm_colnames_reservoir)){
				return(self$compute_treatment_estimate())
			}

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			X_data = private$get_X()

			# Matched pairs component
			beta_m = NA_real_
			ssq_m = NA_real_
			if (m > 0 && !is.null(private$best_Xmm_colnames_matched)){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_matched = which(m_vec > 0L)

				X_cov = X_data[i_matched, intersect(private$best_Xmm_colnames_matched, colnames(X_data)), drop = FALSE]
				Xmm = cbind(w = private$w[i_matched], X_cov)

				fit_m = .fit_clayton_weibull_aft(
					y = private$y[i_matched],
					dead = private$dead[i_matched],
					Xmm = Xmm,
					pair_id = m_vec[i_matched],
					include_singletons = FALSE,
					starts = list(private$best_par_matched)
				)
				if (!is.null(fit_m) && is.finite(fit_m$beta) && is.finite(fit_m$ssq) && fit_m$ssq > 0){
					beta_m = fit_m$beta
					ssq_m = fit_m$ssq
				}
			}

			# Reservoir component
			beta_r = NA_real_
			ssq_r = NA_real_
			if (nRT > 0 && nRC > 0 && !is.null(private$best_Xmm_colnames_reservoir)){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_reservoir = which(m_vec == 0L)

				X_cov = X_data[i_reservoir, intersect(private$best_Xmm_colnames_reservoir, colnames(X_data)), drop = FALSE]
				Xmm = cbind(w = private$w[i_reservoir], X_cov)

				fit_r = .fit_standard_weibull_aft_from_matrix(
					y = private$y[i_reservoir],
					dead = private$dead[i_reservoir],
					Xmm = Xmm
				)
				if (!is.null(fit_r) && is.finite(fit_r$beta) && is.finite(fit_r$ssq) && fit_r$ssq > 0){
					beta_r = fit_r$beta
					ssq_r = fit_r$ssq
				}
			}

			# Inverse-variance weighted pooling
			m_ok = is.finite(beta_m) && is.finite(ssq_m)
			r_ok = is.finite(beta_r) && is.finite(ssq_r)

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		}
	),

	private = list(

		filtered_cov_cache = NULL,
		best_Xmm_colnames_matched = NULL,
		best_par_matched = NULL,
		best_Xmm_colnames_reservoir = NULL,

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		# Pre-calculate filtered covariate matrices to avoid redundant correlation checks
		# during repeated calls (e.g., in randomization inference).
		filtered_covariate_candidates = function(){
			if (!is.null(private$filtered_cov_cache)){
				return(private$filtered_cov_cache)
			}

			X_cov_orig = as.matrix(private$X)
			if (is.null(colnames(X_cov_orig))){
				colnames(X_cov_orig) = paste0("x", seq_len(ncol(X_cov_orig)))
			}

			thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			candidates = list()
			keys = character()

			for (thresh in thresholds){
				X_cov = if (is.finite(thresh)) drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M else X_cov_orig
				key = if (ncol(X_cov) > 0) paste(colnames(X_cov), collapse = "|") else ""
				if (!(key %in% keys)){
					candidates[[length(candidates) + 1L]] = X_cov
					keys = c(keys, key)
				}
			}
			private$filtered_cov_cache = candidates
			candidates
		},

		design_matrix_candidates = function(){
			if (!is.null(private$cached_values$clayton_design_candidates)){
				return(private$cached_values$clayton_design_candidates)
			}

			if (!private$include_covariates() || is.null(private$X) || ncol(private$X) == 0L){
				M = matrix(private$w, ncol = 1)
				colnames(M) = "w"
				private$cached_values$clayton_design_candidates = list(M)
				return(private$cached_values$clayton_design_candidates)
			}

			cov_candidates = private$filtered_covariate_candidates()
			candidates = list()

			for (X_cov in cov_candidates){
				if (ncol(X_cov) == 0L){
					M = matrix(private$w, ncol = 1)
					colnames(M) = "w"
				} else {
					M = cbind(w = private$w, X_cov)
					qr_M = qr(M)
					if (qr_M$rank < ncol(M)){
						keep = qr_M$pivot[seq_len(qr_M$rank)]
						if (!(1L %in% keep)) keep = c(1L, keep)
						keep = sort(unique(keep))
						M = M[, keep, drop = FALSE]
					}
					colnames(M)[1] = "w"
				}
				candidates[[length(candidates) + 1L]] = M
			}

			private$cached_values$clayton_design_candidates = candidates
			candidates
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			if (m > 0){
				private$clayton_copula_for_matched_pairs()
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
			       !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			if (nRT > 0 && nRC > 0){
				private$weibull_for_reservoir()
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) &&
			       !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		clayton_copula_for_matched_pairs = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0L)
			if (length(i_matched) == 0L) return(invisible(NULL))

			for (Xcand in private$design_matrix_candidates()){
				fit = .fit_clayton_weibull_aft(
					y = private$y[i_matched],
					dead = private$dead[i_matched],
					Xmm = Xcand[i_matched, , drop = FALSE],
					pair_id = m_vec[i_matched],
					include_singletons = FALSE
				)
				if (!is.null(fit) && is.finite(fit$beta) && is.finite(fit$ssq) && fit$ssq > 0){
					private$cached_values$beta_T_matched = fit$beta
					private$cached_values$ssq_beta_T_matched = fit$ssq
					private$cached_values$theta_matched = fit$theta
					private$best_par_matched = fit$best_par
					private$best_Xmm_colnames_matched = colnames(Xcand)
					return(invisible(NULL))
				}
			}
		},

		weibull_for_reservoir = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_reservoir = which(m_vec == 0L)
			if (length(i_reservoir) == 0L) return(invisible(NULL))

			for (Xcand in private$design_matrix_candidates()){
				fit = .fit_standard_weibull_aft_from_matrix(
					y = private$y[i_reservoir],
					dead = private$dead[i_reservoir],
					Xmm = Xcand[i_reservoir, , drop = FALSE]
				)
				if (!is.null(fit) && is.finite(fit$beta) && is.finite(fit$ssq) && fit$ssq > 0){
					private$cached_values$beta_T_reservoir = fit$beta
					private$cached_values$ssq_beta_T_reservoir = fit$ssq
					private$best_Xmm_colnames_reservoir = colnames(Xcand)
					return(invisible(NULL))
				}
			}
		}
	)
)
