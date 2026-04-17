#' Abstract class for GLMM-based Inference
#'
#' @keywords internal
InferenceAbstractKKGLMM = R6::R6Class("InferenceAbstractKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), private$glmm_response_type())
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or FixedDesignBinaryMatch.")
				}
			}
			super$initialize(des_obj, verbose)
			if (is(des_obj, "FixedDesignBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_bms_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			if (identical(private$glmm_response_type(), "proportion")) {
				private$y = .sanitize_proportion_response(private$y, interior = FALSE)
			}
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts()) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
				}
			}
		},

		#' @description
		#' Returns the estimated treatment effect.
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
		#' Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
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
		},

		#' @description
		#' Overridden to provide a warning about slowness.
		#' @param r		Number of randomization iterations.
		#' @param ... 					Additional arguments passed to super.
		compute_confidence_interval_rand = function(r = 501, ...){
			warning("Randomization-based confidence intervals for GLMM models are extremely slow because they require many hundreds or thousands of model fits. Consider using asymptotic or bootstrap intervals instead.")
			super$compute_confidence_interval_rand(r = r, ...)
		},

		#' @description
		#' Overridden to provide a warning about slowness.
		#' @param r		Number of randomization iterations.
		#' @param ... 					Additional arguments passed to super.
		compute_two_sided_pval_for_treatment_effect_rand = function(r = 501, ...){
			warning("Randomization-based p-values for GLMM models are slow because they require many model fits.")
			super$compute_two_sided_pval_for_treatment_effect_rand(r = r, ...)
		}
	),

	private = list(

		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the fixed-effect coefficient for "w" directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(){
			mod = private$fit_glmm(se = FALSE)
			if (is.null(mod)) return(NA_real_)
			
			# glmmTMB fixed effects for the conditional model
			beta = glmmTMB::fixef(mod)$cond
			if ("w" %in% names(beta)){
				return(as.numeric(beta["w"]))
			}
			NA_real_
		},

		# Abstract: subclasses must return the expected response type string.
		glmm_response_type = function() stop(class(self)[1], " must implement glmm_response_type()"),

		# Abstract: subclasses must return the glm family object for glmmTMB.
		glmm_family = function() stop(class(self)[1], " must implement glmm_family()"),

		# Default (multivariate): intercept dropped, treatment column named "w".
		# Univariate subclasses override this to return data.frame(w = private$w).
		glmm_predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},

		glmm_predictors_df_candidates = function(){
			predictors_df = private$glmm_predictors_df()
			if (!private$harden || is.null(predictors_df) || ncol(predictors_df) <= 1L){
				return(list(predictors_df))
			}

			X_cov_orig = as.matrix(predictors_df[, setdiff(colnames(predictors_df), "w"), drop = FALSE])
			if (ncol(X_cov_orig) == 0L){
				return(list(data.frame(w = predictors_df$w)))
			}

			thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			candidates = list()
			keys = character()

			for (thresh in thresholds){
				X_cov = if (is.finite(thresh)) drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M else X_cov_orig
				X_full = cbind(w = predictors_df$w, X_cov)
				reduced = qr_reduce_preserve_cols_cpp(as.matrix(X_full), required_cols = 1L)
				X_red = reduced$X_reduced
				keep_idx = as.integer(reduced$keep)
				if (length(keep_idx) == 0L) next
				colnames(X_red) = colnames(X_full)[keep_idx]
				candidate_df = as.data.frame(X_red, check.names = FALSE)
				key = paste(colnames(candidate_df), collapse = "|")
				if (!(key %in% keys)){
					candidates[[length(candidates) + 1L]] = candidate_df
					keys = c(keys, key)
				}
			}

			if (!("w" %in% keys)){
				candidates[[length(candidates) + 1L]] = data.frame(w = predictors_df$w)
			}
			candidates
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			mod = private$fit_glmm(se = TRUE)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_glmm_fit_unavailable")
				return(invisible(NULL))
			}
			coef_table = summary(mod)$coefficients$cond
			private$cached_values$beta_hat_T   = as.numeric(coef_table["w", "Estimate"])
			se = as.numeric(coef_table["w", "Std. Error"])
			# Store NA when the SE is non-finite; SE-dependent methods detect this via assert_finite_se()
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z         = TRUE
		},

		assert_finite_se = function(){
			se = private$cached_values$s_beta_hat_T
			if (should_run_asserts()) {
				if (!is.finite(se) || se <= 0){
					stop("GLMM inference produced a non-finite standard error.")
				}
			}
			invisible(NULL)
		},

		fit_glmm = function(se = TRUE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Build group ID: matched pairs share their m_vec value;
			# reservoir subjects (m_vec == 0) each get a unique singleton ID.
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)

			for (predictors_df in private$glmm_predictors_df_candidates()){
				mod = private$.fit_one_glmm_candidate(predictors_df, group_id, glmm_control, se)
				if (private$.is_usable_glmm_fit(mod, se)) return(mod)
			}
			NULL
		},

		.fit_one_glmm_candidate = function(predictors_df, group_id, glmm_control, se){
			dat = data.frame(y = private$y, predictors_df, group_id = factor(group_id))
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			glmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))
			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					glmmTMB::glmmTMB(
						glmm_formula,
						family  = private$glmm_family(),
						data    = dat,
						control = glmm_control,
						se      = se
					)
				)))
				mod
			}, error = function(e) NULL)
		},

		.is_usable_glmm_fit = function(mod, se){
			if (is.null(mod)) return(FALSE)
			beta = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
			if (is.null(beta) || !("w" %in% names(beta)) || !is.finite(beta["w"])) return(FALSE)
			if (!se) return(TRUE)
			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table)) || !("Std. Error" %in% colnames(coef_table))) return(FALSE)
			se_w = suppressWarnings(as.numeric(coef_table["w", "Std. Error"]))
			is.finite(se_w) && se_w > 0
		}
	)
)
