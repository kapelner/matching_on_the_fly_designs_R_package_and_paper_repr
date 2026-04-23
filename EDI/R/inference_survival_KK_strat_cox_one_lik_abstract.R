#' Abstract class for Stratified Cox Combined-Likelihood Inference
#'
#' Fits a single joint stratified Cox model over all KK design data for survival
#' responses. Matched subjects share their pair ID as a stratum; reservoir subjects
#' are assigned to a single shared "reservoir" stratum. The treatment effect beta_T
#' and covariate slopes beta_xs are shared across all strata.
#'
#' Under \code{harden = TRUE}, the multivariate fit preserves the treatment column
#' and retries reduced covariate sets after QR-based rank reduction. Extreme finite
#' coefficients / standard errors are rejected and treated as non-estimable.
#'
#' @keywords internal
InferenceAbstractKKStratCoxOneLik = R6::R6Class("InferenceAbstractKKStratCoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		},

		#' @description
		#' Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only Whether to skip standard-error calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an asymptotic confidence interval for the treatment effect.
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
		#' Returns a 2-sided p-value for H0: beta_T = delta.
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

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared_combined_likelihood(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()

			X_cov = if (length(Xmm_cols) == 0L){
				matrix(0, nrow = private$n, ncol = 0)
			} else {
				X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
			}
			X_fit = if (ncol(X_cov) > 0L) cbind(w = private$w, X_cov) else matrix(private$w, ncol = 1L, dimnames = list(NULL, "w"))

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			strata_id = m_vec
			strata_id[strata_id == 0L] = max(strata_id) + 1L

			res = tryCatch(
				fast_stratified_coxph_regression_cpp(private$y, private$dead, X_fit,
				                                     as.integer(strata_id), estimate_only = TRUE),
				error = function(e) NULL
			)
			if (is.null(res)) return(NA_real_)
			res$coefficients[1L]
		},

		build_design_matrix = function(){
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			if (ncol(as.matrix(private$X)) > 0 && ncol(private$get_X()) > 0L){
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

			strata_id = m_vec
			strata_id[strata_id == 0L] = max(strata_id) + 1L
			strata_int = as.integer(strata_id)

			if (sum(private$dead) == 0L){
				private$cache_nonestimable_estimate("kk_strat_cox_combined_no_events")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					tryCatch(
						fast_stratified_coxph_regression_cpp(private$y, private$dead, X_fit, strata_int),
						error = function(e) NULL
					)
				},
				fit_ok = function(res, X_fit, keep){
					if (is.null(res)) return(FALSE)
					beta = tryCatch(res$coefficients[1L], error = function(e) NA_real_)
					se   = tryCatch(sqrt(res$vcov[1L, 1L]), error = function(e) NA_real_)
					is.finite(beta) && abs(beta) <= private$max_abs_reasonable_coef &&
						is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
				}
			)
			res = attempt$fit
			if (!is.null(res)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), "w")
			}
			if (is.null(res)){
				private$cache_nonestimable_estimate("kk_strat_cox_combined_fit_failed")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			beta = tryCatch(res$coefficients[1L], error = function(e) NA_real_)
			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			if (!is.finite(private$cached_values$beta_hat_T) || abs(private$cached_values$beta_hat_T) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_strat_cox_combined_extreme_estimate")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			if (!estimate_only) {
				se = tryCatch(sqrt(res$vcov[1L, 1L]), error = function(e) NA_real_)
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef) se else NA_real_
				if (!is.finite(private$cached_values$s_beta_hat_T)){
					private$cache_nonestimable_se("kk_strat_cox_combined_standard_error_unavailable")
					private$cached_values$is_z = TRUE
					return(invisible(NULL))
				}
			}
			private$clear_nonestimable_state()
			private$cached_values$is_z = TRUE
			invisible(NULL)
		}
	)
)
