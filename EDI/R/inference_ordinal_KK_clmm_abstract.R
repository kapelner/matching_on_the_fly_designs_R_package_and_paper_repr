#' Abstract class for ordinal CLMM-based Inference in KK designs
#'
#' @keywords internal
InferenceAbstractKKOrdinalCLMM = R6::R6Class("InferenceAbstractKKOrdinalCLMM",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj,  verbose = FALSE, harden = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose, harden)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts()) {
				if (!check_package_installed("ordinal")){
					stop("Package 'ordinal' is required for ", class(self)[1], ". Please install it.")
				}
			}
		},

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute asymp two sided pval for treatment effect
		#' @param delta Description for delta
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			if (should_run_asserts()) {
				if (delta == 0){
					private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
				} else {
					stop("TO-DO")
				}
			}
		}
	),

	private = list(
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			
			X_fit = if (length(Xmm_cols) == 0L){
				# Univariate case
				matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				cbind(treatment = private$w, X_cov)
			}
			# Add intercept column for clmm internal expectations
			X_fit_full = cbind("(Intercept)" = 1, X_fit)

			mod = private$fit_clmm(X_fit_full)
			if (is.null(mod)) mod = private$fit_clm_fallback(X_fit_full)
			if (is.null(mod)) return(NA_real_)
			
			as.numeric(stats::coef(mod)["w"])
		},

		clmm_link = function() stop(class(self)[1], " must implement clmm_link()"),

		clmm_predictors_df = function(){
			full_X = private$create_design_matrix()
			private$clmm_predictors_df_from_design(full_X)
		},

		clmm_predictors_df_from_design = function(full_X){
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			full_X = private$create_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = full_X,
				required_cols = c(1L, 2L),
				fit_fun = function(X_fit){
					mod = private$fit_clmm(X_fit)
					summ = if (!is.null(mod)) tryCatch(summary(mod), error = function(e) NULL) else NULL
					se = if (!is.null(summ)) as.numeric(summ$coefficients["w", "Std. Error"]) else NA_real_
					if (is.null(mod) || (!estimate_only && (!is.finite(se) || se <= 0))){
						mod = private$fit_clm_fallback(X_fit)
						summ = if (!is.null(mod)) tryCatch(summary(mod), error = function(e) NULL) else NULL
						se = if (!is.null(summ)) as.numeric(summ$coefficients["w", "Std. Error"]) else NA_real_
					}
					list(mod = mod, summ = summ, se = se)
				},
				fit_ok = function(fit, X_fit, keep){
					if (is.null(fit) || is.null(fit$mod)) return(FALSE)
					beta = tryCatch(as.numeric(stats::coef(fit$mod)["w"]), error = function(e) NA_real_)
					if (!is.finite(beta)) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(fit$se) && fit$se > 0
				}
			)
			mod = attempt$fit$mod
			summ = attempt$fit$summ
			se = attempt$fit$se
			if (!is.null(mod)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), c("(Intercept)", "treatment"))
			}
			if (is.null(mod) || is.null(summ)){
				private$cached_values$beta_hat_T   = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(stats::coef(mod)["w"])
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				return(invisible(NULL))
		},

		fit_clmm = function(full_X = private$create_design_matrix()){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			dat = data.frame(
				y = factor(private$y, ordered = TRUE),
				private$clmm_predictors_df_from_design(full_X),
				group_id = factor(group_id)
			)
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			clmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					ordinal::clmm(
						clmm_formula,
						data = dat,
						link = private$clmm_link()
					)
				)))
				mod
			}, error = function(e) NULL)
		},

		fit_clm_fallback = function(full_X = private$create_design_matrix()){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			dat = data.frame(
				y = factor(private$y, ordered = TRUE),
				private$clmm_predictors_df_from_design(full_X)
			)
			fixed_terms = setdiff(colnames(dat), "y")
			clm_formula = stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					ordinal::clm(
						clm_formula,
						data = dat,
						link = private$clmm_link()
					)
				)))
				mod
			}, error = function(e) NULL)
		}
	)
)
