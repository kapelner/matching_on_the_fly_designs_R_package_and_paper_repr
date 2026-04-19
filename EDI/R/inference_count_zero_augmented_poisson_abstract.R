#' Zero-Augmented Poisson Inference for Count Responses
#'
#' Internal base class for non-KK zero-inflated and hurdle Poisson regression
#' models fit using the \pkg{glmmTMB} fitter. The reported treatment effect is the
#' treatment coefficient from the conditional count component, on the log-rate
#' scale.
#'
#' @keywords internal
#' @noRd
InferenceCountZeroAugmentedPoissonAbstract = R6::R6Class("InferenceCountZeroAugmentedPoissonAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the optimized Rcpp
		#'   implementation. If \code{FALSE}, use \pkg{glmmTMB}.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj, include_covariates = NULL, use_rcpp = TRUE, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(include_covariates, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			super$initialize(des_obj, verbose = verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
			private$use_rcpp = use_rcpp
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
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
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
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		}
	),

	private = list(
		best_Xmm_colnames = NULL,
		include_covariates = NULL,
		use_rcpp = TRUE,

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
			
			if (length(Xmm_cols) == 0L){
				# Univariate case
				X_fit = cbind(1, treatment = private$w)
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				X_fit = cbind(1, treatment = private$w, X_cov)
			}

			if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				# Xcond and Xzi are the same here
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(y = as.numeric(private$y), Xcond = X_fit, Xzi = X_fit, is_hurdle = is_hurdle, estimate_only = estimate_only),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				return(as.numeric(fit$coefficients$cond[2]))
			} else {
				dat = if (length(Xmm_cols) == 0L) data.frame(y = private$y, w = private$w) else data.frame(y = private$y, w = private$w, X_cov)
				mod = private$fit_zero_augmented_model(dat)
				if (is.null(mod)) return(NA_real_)
				
				cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (is.null(cond_coef) || !("w" %in% names(cond_coef))) return(NA_real_)
				return(as.numeric(cond_coef["w"]))
			}
		},

		za_family = function() stop(class(self)[1], " must implement za_family()."),

		za_description = function() stop(class(self)[1], " must implement za_description()."),

		predictors_df = function(){
			if (private$include_covariates) {
				full_X = private$create_design_matrix()
				X_model = full_X[, -1, drop = FALSE]
				colnames(X_model)[1] = "w"
				as.data.frame(X_model)
			} else {
				data.frame(w = private$w)
			}
		},

		build_formula = function(dat){
			fixed_terms = setdiff(colnames(dat), "y")
			stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))
		},

		build_zi_formula = function(dat){
			fixed_terms = setdiff(colnames(dat), "y")
			stats::as.formula(paste("~", paste(fixed_terms, collapse = " + ")))
		},

		fit_zero_augmented_model = function(dat){
			formula_cond = private$build_formula(dat)
			formula_zi = private$build_zi_formula(dat)

			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)

			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = formula_zi,
						family = private$za_family(),
						data = dat,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)

			if (!is.null(mod)) return(mod)
			if (ncol(dat) <= 2L) return(NULL)

			dat_fallback = dat[, c("y", "w"), drop = FALSE]
			tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						y ~ w,
						ziformula = ~ w,
						family = private$za_family(),
						data = dat_fallback,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (is.null(private$best_Xmm_colnames)) {
				X_full = cbind(1, private$w, as.matrix(private$predictors_df()[, setdiff(colnames(private$predictors_df()), "w"), drop = FALSE]))
				colnames(X_full)[1:2] = c("(Intercept)", "w")
				
				# Keep track of original names
				original_colnames = colnames(X_full)
				
				res_reduced = private$reduce_design_matrix_preserving_treatment(X_full)
				X_reduced = res_reduced$X
				if (is.null(X_reduced)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_design_unusable")
					return(invisible(NULL))
				}
				# Restore names for X_reduced
				colnames(X_reduced) = original_colnames[res_reduced$keep]
				
				X_fit = X_reduced
				private$best_Xmm_colnames = setdiff(colnames(X_reduced), c("(Intercept)", "w"))
			} else {
				X_data = private$get_X()
				Xmm_cols = private$best_Xmm_colnames
				if (length(Xmm_cols) == 0L){
					X_fit = cbind(1, w = private$w)
				} else {
					X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
					X_fit = cbind(1, w = private$w, X_cov)
				}
			}

			if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(y = as.numeric(private$y), Xcond = X_fit, Xzi = X_fit, is_hurdle = is_hurdle, estimate_only = estimate_only),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) {
					private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
					return(invisible(NULL))
				}
				private$cached_mod = fit
				
				private$cached_values$beta_hat_T = as.numeric(fit$coefficients$cond[2])
				if (!estimate_only) {
					se = tryCatch(sqrt(fit$vcov[2, 2]), error = function(e) NA_real_)
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
			} else {
				dat = data.frame(y = private$y, X_fit[, -1, drop = FALSE])
				mod = private$fit_zero_augmented_model(dat)
				if (is.null(mod)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
					return(invisible(NULL))
				}
				private$cached_mod = mod

				cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (is.null(cond_coef) || !("w" %in% names(cond_coef))){
					private$cache_nonestimable_estimate("zero_augmented_poisson_treatment_missing")
					return(invisible(NULL))
				}

				private$cached_values$beta_hat_T = as.numeric(cond_coef["w"])
				if (!estimate_only) {
					coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
					se = if (!is.null(coef_table) && ("w" %in% rownames(coef_table))) as.numeric(coef_table["w", "Std. Error"]) else NA_real_
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)
)
