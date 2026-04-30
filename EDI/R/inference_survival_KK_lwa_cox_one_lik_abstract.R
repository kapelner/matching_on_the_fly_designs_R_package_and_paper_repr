#' Abstract class for LWA-style Marginal Cox Combined-Likelihood Inference
#'
#' Fits a single joint marginal Cox model over all KK design data for survival
#' responses. Matched subjects share their pair ID as a cluster, and reservoir
#' subjects are treated as independent (unique clusters). Standard errors are
#' obtained via the Huber-White cluster-robust sandwich estimator (LWA style).
#'
#' @keywords internal
InferenceAbstractKKLWACoxOneLik = R6::R6Class("InferenceAbstractKKLWACoxOneLik",
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
		#' Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute an asymptotic confidence interval.
		#' @param alpha Significance level.
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
		#' Compute an asymptotic two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
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
		optimization_alg = "lbfgs",

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Re-read w, y, dead because they might have been transformed for randomization
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			
			# Recompute basic match data for the new w/y/dead
			private$compute_basic_match_data()
			
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared_combined_likelihood(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(NA_real_)
			}

			X_data = private$get_X()
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			X_covs_filtered = X_data[, intersect(private$best_Xmm_colnames, colnames(X_data)), drop = FALSE]
			if (ncol(X_covs_filtered) > 0){
				X_full = cbind(X_full, X_covs_filtered)
			}
			
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Clusters: pairs + reservoir
			cluster_ids = m_vec
			res_idx = which(cluster_ids == 0L)
			if (length(res_idx) > 0L){
				max_m = if (any(cluster_ids > 0)) max(cluster_ids) else 0L
				cluster_ids[res_idx] = max_m + seq_along(res_idx)
			}

			fit = tryCatch(
				fast_coxph_regression_cpp(
					y = private$y,
					dead = private$dead,
					X = as.matrix(X_full),
					cluster = as.integer(cluster_ids),
					estimate_only = estimate_only,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || length(fit$coefficients) < 1 || !is.finite(fit$coefficients[1])) return(NA_real_)
			as.numeric(fit$coefficients[1])
		},

		design_matrix_candidates = function(){
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			X_covs = private$get_X()
			if (ncol(as.matrix(X_covs)) > 0){
				X_full = cbind(X_full, as.matrix(X_covs))
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

			# Clusters: pairs + reservoir
			cluster_ids = m_vec
			res_idx = which(cluster_ids == 0L)
			if (length(res_idx) > 0L){
				max_m = max(cluster_ids)
				cluster_ids[res_idx] = max_m + seq_along(res_idx)
			}

			X_full = private$design_matrix_candidates()

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L, # treatment
				fit_fun = function(X_fit){
					fast_coxph_regression_cpp(
						y = private$y,
						dead = private$dead,
						X = as.matrix(X_fit),
						cluster = as.integer(cluster_ids),
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					)
				},
				fit_ok = function(res, X_fit, keep){
					if (is.null(res) || !isTRUE(res$converged)) return(FALSE)
					beta = res$coefficients[1L]
					if (!is.finite(beta) || abs(beta) > private$max_abs_reasonable_coef) return(FALSE)
					if (estimate_only) return(TRUE)
					se   = tryCatch(sqrt(res$vcov[1L, 1L]), error = function(e) NA_real_)
					is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
				}
			)
			res = attempt$fit
			if (!is.null(res)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), "w")
				private$cached_values$beta_hat_T = as.numeric(res$coefficients[1])
				if (!estimate_only) {
					se = sqrt(res$vcov[1L, 1L])
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
				return(invisible(NULL))
			}

			private$cache_nonestimable_estimate("kk_lwa_cox_combined_fit_failed")
			invisible(NULL)
		}
	)
)
