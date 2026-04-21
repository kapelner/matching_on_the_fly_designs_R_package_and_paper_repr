#' Abstract class for Conditional Logistic Compound Inference for KK Designs
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' binary (incidence) responses using conditional logistic regression for matched pairs.
#' The matched-pair component uses the conditional likelihood to account for pair-level
#' fixed effects. The reservoir component uses standard logistic regression.
#' The two treatment-effect estimates are combined by inverse-variance weighting.
#'
#' @keywords internal
InferenceAbstractKKClogitIVWC = R6::R6Class("InferenceAbstractKKClogitIVWC",
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
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		},

		#' @description
		#' Returns the estimated treatment effect (log-odds ratio).
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
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
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
			# Ensure we have the best design and parameters from the original data
			if (is.null(private$best_Xmm_colnames_matched) && is.null(private$best_Xmm_colnames_reservoir)){
				private$shared()
			}

			# If we still don't have enough (e.g., initial fit failed), fall back to standard
			if (is.null(private$best_Xmm_colnames_matched) && is.null(private$best_Xmm_colnames_reservoir)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			X_data = private$get_X()

			# Matched pairs component (Clogit)
			beta_m = NA_real_
			ssq_m = NA_real_
			if (m > 0 && !is.null(private$best_Xmm_colnames_matched)){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_matched = which(m_vec > 0L)

				# Discordant pairs check
				y_matched = private$y[i_matched]
				pair_ids_matched = m_vec[i_matched]
				
				# Only discordant pairs contribute to clogit
				disc_pair_ids = integer(0)
				for (pid in unique(pair_ids_matched)){
					if (length(unique(y_matched[pair_ids_matched == pid])) > 1L){
						disc_pair_ids = c(disc_pair_ids, pid)
					}
				}
				
				if (length(disc_pair_ids) > 0L){
					i_disc = which(pair_ids_matched %in% disc_pair_ids)
					X_cov = X_data[i_matched[i_disc], intersect(private$best_Xmm_colnames_matched, colnames(X_data)), drop = FALSE]
					Xmm = cbind(w = private$w[i_matched[i_disc]], X_cov)

					fit_m = .fit_clogit(
						y = y_matched[i_disc],
						Xmm = Xmm,
						pair_id = pair_ids_matched[i_disc],
						estimate_only = estimate_only
					)
					if (!is.null(fit_m) && is.finite(fit_m$beta)){
						beta_m = fit_m$beta
						if (!estimate_only && !is.null(fit_m$ssq) && is.finite(fit_m$ssq) && fit_m$ssq > 0){
							ssq_m = fit_m$ssq
						}
					}
				}
			}

			# Reservoir component (Logistic)
			beta_r = NA_real_
			ssq_r = NA_real_
			if (nRT > 0 && nRC > 0 && !is.null(private$best_Xmm_colnames_reservoir)){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_reservoir = which(m_vec == 0L)

				X_cov = X_data[i_reservoir, intersect(private$best_Xmm_colnames_reservoir, colnames(X_data)), drop = FALSE]
				Xmm = cbind(`(Intercept)` = 1, w = private$w[i_reservoir], X_cov)

				fit_r = .fit_logistic_from_matrix(
					y = private$y[i_reservoir],
					Xmm = Xmm,
					estimate_only = estimate_only
				)
				if (!is.null(fit_r) && is.finite(fit_r$beta)){
					beta_r = fit_r$beta
					if (!estimate_only && !is.null(fit_r$ssq) && is.finite(fit_r$ssq) && fit_r$ssq > 0){
						ssq_r = fit_r$ssq
					}
				}
			}

			# Pooling
			m_ok = is.finite(beta_m) && (estimate_only || is.finite(ssq_m))
			r_ok = is.finite(beta_r) && (estimate_only || is.finite(ssq_r))

			if (m_ok && r_ok){
				if (estimate_only) {
					ssq_m_orig = private$cached_values$ssq_beta_T_matched
					ssq_r_orig = private$cached_values$ssq_beta_T_reservoir
					if (is.finite(ssq_m_orig) && is.finite(ssq_r_orig)){
						w_star = ssq_r_orig / (ssq_r_orig + ssq_m_orig)
						return(w_star * beta_m + (1 - w_star) * beta_r)
					}
					return(0.5 * beta_m + 0.5 * beta_r)
				}
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		},

		best_Xmm_colnames_matched = NULL,
		best_Xmm_colnames_reservoir = NULL,

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
				private$clogit_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
			       (!estimate_only && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0 || estimate_only)

			if (nRT > 0 && nRC > 0){
				private$logistic_for_reservoir(estimate_only = estimate_only)
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

		clogit_for_matched_pairs = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0L)
			if (length(i_matched) == 0L) return(invisible(NULL))

			if (ncol(as.matrix(private$X)) == 0L){
				Xmm = matrix(private$w[i_matched], ncol = 1L)
				colnames(Xmm) = "w"
				fit = .fit_clogit(
					y = private$y[i_matched],
					Xmm = Xmm,
					pair_id = m_vec[i_matched],
					estimate_only = estimate_only
				)
				if (!is.null(fit) && is.finite(fit$beta) && (isTRUE(estimate_only) || (is.finite(fit$ssq) && fit$ssq > 0))){
					private$cached_values$beta_T_matched = fit$beta
					private$cached_values$ssq_beta_T_matched = fit$ssq
					private$best_Xmm_colnames_matched = "w"
					return(invisible(NULL))
				}
			} else {
				X_data = private$get_X()
				X_matched = X_data[i_matched, , drop = FALSE]
				X_full = cbind(w = private$w[i_matched], X_matched)
				
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = 1L, # treatment only
					fit_fun = function(X_fit){
						.fit_clogit(
							y = private$y[i_matched],
							Xmm = X_fit,
							pair_id = m_vec[i_matched],
							estimate_only = estimate_only
						)
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod) || !is.finite(mod$beta)) return(FALSE)
						if (estimate_only) return(TRUE)
						is.finite(mod$ssq) && mod$ssq > 0
					}
				)
				
				if (!is.null(attempt$fit)){
					private$cached_values$beta_T_matched = attempt$fit$beta
					private$cached_values$ssq_beta_T_matched = attempt$fit$ssq
					private$best_Xmm_colnames_matched = colnames(attempt$X)
					return(invisible(NULL))
				}
			}
		},

		logistic_for_reservoir = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_reservoir = which(m_vec == 0L)
			if (length(i_reservoir) == 0L) return(invisible(NULL))

			X_data = private$get_X()
			X_reservoir = X_data[i_reservoir, , drop = FALSE]
			X_full = cbind(`(Intercept)` = 1, w = private$w[i_reservoir], X_reservoir)

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L, # intercept and treatment
				fit_fun = function(X_fit){
					.fit_logistic_from_matrix(
						y = private$y[i_reservoir],
						Xmm = X_fit,
						estimate_only = estimate_only
					)
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || !is.finite(mod$beta)) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq) && mod$ssq > 0
				}
			)
			
			if (!is.null(attempt$fit)){
				private$cached_values$beta_T_reservoir = attempt$fit$beta
				private$cached_values$ssq_beta_T_reservoir = attempt$fit$ssq
				private$best_Xmm_colnames_reservoir = setdiff(colnames(attempt$X), c("(Intercept)", "w"))
				return(invisible(NULL))
			}
		}
	)
)
