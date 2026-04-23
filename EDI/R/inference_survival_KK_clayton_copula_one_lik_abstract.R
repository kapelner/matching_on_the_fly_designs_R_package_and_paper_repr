#' Abstract class for Clayton Copula Combined-Likelihood Inference
#'
#' This class implements a combined-likelihood estimator for KK matching-on-the-fly
#' designs with survival responses using a Clayton copula with Weibull AFT margins.
#' The likelihood includes contributions from both matched pairs (bivariate) and
#' reservoir subjects (univariate).
#'
#' @keywords internal
InferenceAbstractKKClaytonCopulaOneLik = R6::R6Class("InferenceAbstractKKClaytonCopulaOneLik",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		},

		#' @description
		#' Returns the estimated treatment effect (log-time ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. Default is 0.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			inf_obj = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			inf_obj
		}
	),

	private = list(

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared()
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			X_data = private$get_X()
			X_cov = X_data[, intersect(private$best_Xmm_colnames, colnames(X_data)), drop = FALSE]
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = "w"
			if (ncol(X_cov) > 0){
				Xmm = cbind(Xmm, X_cov)
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			fit = .fit_clayton_weibull_aft(
				y = private$y,
				dead = private$dead,
				Xmm = Xmm,
				pair_id = m_vec,
				include_singletons = TRUE,
				starts = list(private$best_par),
				estimate_only = estimate_only
			)
			if (is.null(fit) || !is.finite(fit$beta)){
				return(NA_real_)
			}
			as.numeric(fit$beta)
		},

		best_Xmm_colnames = NULL,
		best_par = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Re-read design variables which might have been transformed during randomization
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			
			# Recompute basic match data for the new w/y/dead
			private$compute_basic_match_data()
			
			# Use the same joint-likelihood logic for the point estimate
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		filtered_covariate_candidates = function(){
			if (ncol(as.matrix(private$X)) == 0L) return(list(matrix(nrow = private$n, ncol = 0L)))
			# Use the helper to get candidates (dropping linearly dependent or high-correlation columns)
			get_reduced_design_matrix_candidates(as.matrix(private$X))
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Optimization candidates
			if (ncol(as.matrix(private$X)) == 0L){
				Xcand = matrix(private$w, ncol = 1)
				colnames(Xcand) = "w"
				candidates = list(Xcand)
			} else {
				cov_candidates = private$filtered_covariate_candidates()
				candidates = list()
				for (X_cov in cov_candidates){
					if (ncol(X_cov) == 0L){
						M = matrix(private$w, ncol = 1)
						colnames(M) = "w"
					} else {
						M = matrix(private$w, ncol = 1)
						colnames(M) = "w"
						M = cbind(M, X_cov)
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
			}

			for (Xmm in candidates){
				fit = .fit_clayton_weibull_aft(
					y = private$y,
					dead = private$dead,
					Xmm = Xmm,
					pair_id = m_vec,
					include_singletons = TRUE,
					estimate_only = estimate_only
				)
				if (!is.null(fit) && is.finite(fit$beta) && (isTRUE(estimate_only) || (is.finite(fit$ssq) && fit$ssq > 0))){
					private$cached_values$beta_hat_T = fit$beta
					private$cached_values$s_beta_hat_T = if (is.finite(fit$ssq) && fit$ssq > 0) sqrt(fit$ssq) else NA_real_
					private$cached_values$theta = fit$theta
					private$best_par = fit$best_par
					private$best_Xmm_colnames = colnames(Xmm)
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
