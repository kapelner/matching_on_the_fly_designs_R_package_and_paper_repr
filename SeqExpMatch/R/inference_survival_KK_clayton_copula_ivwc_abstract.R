# Abstract class for Clayton Copula / Standard Weibull Compound Inference
#
# @description
# This class implements a compound estimator for KK matching-on-the-fly designs with
# survival responses using a Clayton copula with Weibull AFT margins for matched pairs.
# The matched-pair component is fitted by maximizing the pairwise copula likelihood
# under right censoring. The reservoir component uses standard Weibull AFT regression.
# The two treatment-effect estimates are combined by inverse-variance weighting.
#
# @details
# This is the package's first copula-based bivariate survival implementation. The matched-pair
# contribution uses a Clayton survival copula with Weibull AFT marginal models, so the
# treatment effect is reported on the AFT log-time-ratio scale, matching the existing Weibull
# survival classes. The copula dependence parameter is estimated jointly with the matched-pair
# regression coefficients by direct likelihood maximization.
#
# @keywords internal
SeqDesignInferenceAbstractKKClaytonCopulaIVWC = R6::R6Class("SeqDesignInferenceAbstractKKClaytonCopulaIVWC",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "survival")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		# @description
		# Returns the estimated treatment effect (log-time ratio).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes the MLE-based confidence interval.
		# @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes the MLE-based p-value.
		# @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("Testing non-zero delta is not yet implemented for this class.")
			}
		}
	),

	private = list(

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

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

			thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			candidates = list()
			keys = character()
			X_cov_orig = as.matrix(private$X)
			if (is.null(colnames(X_cov_orig))){
				colnames(X_cov_orig) = paste0("x", seq_len(ncol(X_cov_orig)))
			}

			for (thresh in thresholds){
				X_cov = if (is.finite(thresh)) drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M else X_cov_orig
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
				key = paste(colnames(M), collapse = "|")
				if (!(key %in% keys)){
					candidates[[length(candidates) + 1L]] = M
					keys = c(keys, key)
				}
			}

			private$cached_values$clayton_design_candidates = candidates
			private$cached_values$clayton_design_candidates
		},

		shared = function(){
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
				stop("Clayton copula compound estimator: could not compute a finite standard error.")
			}
		},

		clayton_copula_for_matched_pairs = function(){
			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			i_matched = which(match_indic > 0L)
			if (length(i_matched) == 0L) return(invisible(NULL))

			for (Xcand in private$design_matrix_candidates()){
				fit = .fit_clayton_weibull_aft(
					y = private$y[i_matched],
					dead = private$dead[i_matched],
					Xmm = Xcand[i_matched, , drop = FALSE],
					pair_id = match_indic[i_matched],
					include_singletons = FALSE
				)
				if (!is.null(fit) && is.finite(fit$beta) && is.finite(fit$ssq) && fit$ssq > 0){
					private$cached_values$beta_T_matched = fit$beta
					private$cached_values$ssq_beta_T_matched = fit$ssq
					private$cached_values$theta_matched = fit$theta
					return(invisible(NULL))
				}
			}
		},

		weibull_for_reservoir = function(){
			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			i_reservoir = which(match_indic == 0L)
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
					return(invisible(NULL))
				}
			}
		}
	)
)
