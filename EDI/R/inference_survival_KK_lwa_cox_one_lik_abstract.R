#' Abstract class for LWA-style Marginal Cox Combined-Likelihood Inference
#'
#' Fits a single joint marginal Cox model over all KK design data for survival
#' responses. Matched subjects share their pair ID as a cluster, and reservoir
#' subjects are assigned singleton cluster IDs. The treatment effect is estimated by
#' the Cox partial likelihood on all data, with Lee-Wei-Amato style cluster-robust
#' variance treating within-pair correlation as a nuisance.
#'
#' Under \code{harden = TRUE}, the multivariate fit preserves the treatment column
#' and retries reduced covariate sets after QR-based rank reduction. Extreme finite
#' coefficients / standard errors are rejected and treated as non-estimable.
#'
#' @keywords internal
InferenceAbstractKKLWACoxOneLik = R6::R6Class("InferenceAbstractKKLWACoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize LWA-style combined-likelihood Cox inference for KK survival data.
		#' @param des_obj A completed KK design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new inference object.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose)
		},

		#' @description
		#' Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		#' @return The treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute an asymptotic confidence interval.
		#' @param alpha Significance level.
		#' @return A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute an asymptotic two-sided p-value for the treatment effect.
		#' @param delta Null treatment effect value.
		#' @return A two-sided p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		max_abs_reasonable_coef = 1e4,
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared_combined_likelihood(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			
			X_cov = if (length(Xmm_cols) == 0L){
				matrix(0, nrow = private$n, ncol = 0)
			} else {
				X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
			}
			X_fit = cbind(w = private$w, X_cov)

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			cluster_id = m_vec
			reservoir_idx = which(cluster_id == 0L)
			if (length(reservoir_idx) > 0L){
				cluster_id[reservoir_idx] = max(cluster_id) + seq_along(reservoir_idx)
			}

			dat = data.frame(y = private$y, dead = private$dead, w = X_fit[, "w"], cluster = factor(cluster_id))
			formula_str = "survival::Surv(y, dead) ~ w"
			if (ncol(X_cov) > 0){
				colnames(X_cov) = paste0("x", seq_len(ncol(X_cov)))
				dat = cbind(dat, X_cov)
				formula_str = paste(formula_str, "+", paste(colnames(X_cov), collapse = " + "))
			}
			formula_str = paste(formula_str, "+ cluster(cluster)")
			
			mod = tryCatch(
				suppressWarnings(survival::coxph(as.formula(formula_str), data = dat, robust = TRUE)),
				error = function(e) NULL
			)
			if (is.null(mod)) return(NA_real_)
			
			coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
			if (is.null(coefs) || !("w" %in% names(coefs))) return(NA_real_)
			as.numeric(coefs["w"])
		},

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		build_design_matrix = function(){
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			if (private$include_covariates() && ncol(private$get_X()) > 0L){
				X_full = cbind(X_full, as.matrix(private$get_X()))
			}
			X_full
		},

		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			cluster_id = m_vec
			reservoir_idx = which(cluster_id == 0L)
			if (length(reservoir_idx) > 0L){
				cluster_id[reservoir_idx] = max(cluster_id) + seq_along(reservoir_idx)
			}

			if (sum(private$dead) == 0L){
				private$cache_nonestimable_estimate("kk_lwa_cox_combined_no_events")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					dat = data.frame(y = private$y, dead = private$dead, w = X_fit[, "w"], cluster = factor(cluster_id))
					formula_str = "survival::Surv(y, dead) ~ w"
					X_covs = X_fit[, colnames(X_fit) != "w", drop = FALSE]
					if (ncol(X_covs) > 0){
						colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
						dat = cbind(dat, X_covs)
						formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
					}
					formula_str = paste(formula_str, "+ cluster(cluster)")
					mod = tryCatch(
						suppressWarnings(survival::coxph(as.formula(formula_str), data = dat, robust = TRUE)),
						error = function(e) NULL
					)
					if (is.null(mod) && ncol(X_covs) > 0){
						mod = tryCatch(
							suppressWarnings(survival::coxph(survival::Surv(y, dead) ~ w + cluster(cluster), data = dat, robust = TRUE)),
							error = function(e) NULL
						)
					}
					mod
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod)) return(FALSE)
					coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
					vcv = tryCatch(stats::vcov(mod), error = function(e) NULL)
					if (is.null(coefs) || is.null(vcv) || !("w" %in% names(coefs)) || !("w" %in% rownames(vcv))) return(FALSE)
					beta = suppressWarnings(as.numeric(coefs["w"]))
					se = suppressWarnings(sqrt(as.numeric(vcv["w", "w"])))
					is.finite(beta) && abs(beta) <= private$max_abs_reasonable_coef &&
						is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
				}
			)
			mod = attempt$fit
			if (!is.null(mod)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), "w")
			}
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_lwa_cox_combined_fit_failed")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
			if (is.null(coefs) || !("w" %in% names(coefs))){
				private$cache_nonestimable_estimate("kk_lwa_cox_combined_treatment_missing")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			beta = as.numeric(coefs["w"])
			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			if (!is.finite(private$cached_values$beta_hat_T) || abs(private$cached_values$beta_hat_T) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_lwa_cox_combined_extreme_estimate")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}
			
			if (!estimate_only) {
				vcv = tryCatch(stats::vcov(mod), error = function(e) NULL)
				if (is.null(vcv) || !("w" %in% rownames(vcv))){
					private$cache_nonestimable_se("kk_lwa_cox_combined_standard_error_unavailable")
					private$cached_values$is_z = TRUE
					return(invisible(NULL))
				} else {
					se = sqrt(as.numeric(vcv["w", "w"]))
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef) se else NA_real_
					if (!is.finite(private$cached_values$s_beta_hat_T)){
						private$cache_nonestimable_se("kk_lwa_cox_combined_standard_error_unavailable")
						private$cached_values$is_z = TRUE
						return(invisible(NULL))
					}
				}
			}
			private$clear_nonestimable_state()
			private$cached_values$is_z = TRUE
			invisible(NULL)
		}
	)
)
