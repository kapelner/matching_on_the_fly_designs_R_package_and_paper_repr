#' Abstract class for Clayton Copula Combined-Likelihood Inference
#'
#' Fits a single joint likelihood over all KK design data for survival responses using
#' a Clayton survival copula with Weibull AFT margins for matched pairs and standard
#' Weibull AFT singleton contributions for reservoir subjects. The treatment effect is
#' reported on the AFT log-time-ratio scale.
#'
#' @details
#' This keeps the package's established \code{CombinedLikelihood} terminology even though
#' the matched and reservoir parts are combined through a joint copula-based likelihood
#' rather than a classical Gaussian-style likelihood. The matched-pair contribution uses
#' a Clayton copula with censored bivariate survival contributions; reservoir subjects
#' contribute ordinary univariate Weibull AFT terms.
#'
#' @keywords internal
InferenceAbstractKKClaytonCopulaCombinedLikelihood = R6::R6Class("InferenceAbstractKKClaytonCopulaCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "survival")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
		},

		#' @description
		#' Returns the combined-likelihood estimate of the treatment effect (log-time ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an asymptotic confidence interval for the treatment effect.
		#' @param alpha Significance level. Default 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Returns a 2-sided p-value for H0: beta_T = delta.
		#' @param delta Null value; default 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		duplicate = function(verbose = FALSE){
			inf_obj = super$duplicate(verbose = verbose)
			inf_obj
		}
	),

	private = list(

		filtered_cov_cache = NULL,
		best_Xmm_colnames = NULL,
		best_par = NULL,

		# Overridden to use the best design matrix and starting parameters from the initial fit
		# to speed up randomization-based inference.
		compute_treatment_estimate_during_randomization_inference = function(){
			# Ensure we have the best design and parameters from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared()
			}

			# If we still don't have it (e.g., initial fit failed), fall back to standard
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate())
			}

			# Use the same design matrix columns as the original fit
			X_data = private$get_X()
			Xmm_cols = private$best_Xmm_colnames
			X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
			Xmm = cbind(w = private$w, X_cov)

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Use the best parameters from the original fit as the ONLY starting point
			# This significantly speeds up each randomization draw.
			fit = .fit_clayton_weibull_aft(
				y = private$y,
				dead = private$dead,
				Xmm = Xmm,
				pair_id = m_vec,
				include_singletons = TRUE,
				starts = list(private$best_par)
			)

			if (!is.null(fit) && is.finite(fit$beta)){
				return(fit$beta)
			}
			NA_real_
			},

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

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			for (Xcand in private$design_matrix_candidates()){
				fit = .fit_clayton_weibull_aft(
					y = private$y,
					dead = private$dead,
					Xmm = Xcand,
					pair_id = m_vec,
					include_singletons = TRUE
				)
				if (!is.null(fit) && is.finite(fit$beta) && is.finite(fit$ssq) && fit$ssq > 0){
					private$cached_values$beta_hat_T = fit$beta
			if (estimate_only) return(invisible(NULL))
					private$cached_values$s_beta_hat_T = sqrt(fit$ssq)
					private$cached_values$theta_hat = fit$theta
					private$best_par = fit$best_par
					private$best_Xmm_colnames = colnames(Xcand)
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
