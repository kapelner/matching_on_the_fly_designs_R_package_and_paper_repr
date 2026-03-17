# Abstract class for Clayton Copula Combined-Likelihood Inference
#
# @description
# Fits a single joint likelihood over all KK design data for survival responses using
# a Clayton survival copula with Weibull AFT margins for matched pairs and standard
# Weibull AFT singleton contributions for reservoir subjects. The treatment effect is
# reported on the AFT log-time-ratio scale.
#
# @details
# This keeps the package's established \code{CombinedLikelihood} terminology even though
# the matched and reservoir parts are combined through a joint copula-based likelihood
# rather than a classical Gaussian-style likelihood. The matched-pair contribution uses
# a Clayton copula with censored bivariate survival contributions; reservoir subjects
# contribute ordinary univariate Weibull AFT terms.
#
# @keywords internal
SeqDesignInferenceAbstractKKClaytonCopulaCombinedLikelihood = R6::R6Class("SeqDesignInferenceAbstractKKClaytonCopulaCombinedLikelihood",
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
		# Returns the combined-likelihood estimate of the treatment effect (log-time ratio).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		# @description
		# Returns a 1 - alpha confidence interval for the treatment effect.
		# @param alpha Significance level; default 0.05 gives a 95 percent CI.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Returns a 2-sided p-value for H0: beta_T = delta.
		# @param delta Null value; default 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
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

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Clayton copula combined-likelihood estimator: could not compute a finite standard error.")
			}
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			for (Xcand in private$design_matrix_candidates()){
				fit = .fit_clayton_weibull_aft(
					y = private$y,
					dead = private$dead,
					Xmm = Xcand,
					pair_id = match_indic,
					include_singletons = TRUE
				)
				if (!is.null(fit) && is.finite(fit$beta) && is.finite(fit$ssq) && fit$ssq > 0){
					private$cached_values$beta_hat_T = fit$beta
					private$cached_values$s_beta_hat_T = sqrt(fit$ssq)
					private$cached_values$theta_hat = fit$theta
					private$cached_values$is_z = TRUE
					return(invisible(NULL))
				}
			}

			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
		}
	)
)
