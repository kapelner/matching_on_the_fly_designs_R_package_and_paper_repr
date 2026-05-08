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
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Whether to use Rcpp speedup.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFormula(model_formula, null.ok = TRUE)
				assertFlag(use_rcpp)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
			
			
			private$use_rcpp = use_rcpp
		},

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha The significance level (default 0.05).
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
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			se = private$compute_standard_error_from_information_matrix()
			if (is.finite(se)) return(se)
			private$cached_values$s_beta_hat_T
		},

		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df
		},

		best_X_colnames = NULL,
		use_rcpp = TRUE,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			X_cols = private$best_X_colnames
			X_data = private$get_X()
			
			if (length(X_cols) == 0L){
				# Univariate case
				X_fit = cbind(1, treatment = private$w)
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X_fit = cbind(1, treatment = private$w, X_cov)
			}

			if (private$use_rcpp && identical(private$za_description(), "Zero-Inflated Negative Binomial")) {
				start_params = private$get_fit_warm_start_for_length("params", 2L * ncol(X_fit) + 1L)
				fit = tryCatch(
					fast_zinb_cpp(X = X_fit, y = as.numeric(private$y), Xzi = X_fit, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				private$set_fit_warm_start(as.numeric(fit$params), "params")
				return(as.numeric(fit$coefficients$cond[2]))
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				# X and Xzi are the same here
				start_params = private$get_fit_warm_start_for_length("params", 2L * ncol(X_fit))
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(X = X_fit, y = as.numeric(private$y), Xzi = X_fit, is_hurdle = is_hurdle, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				private$set_fit_warm_start(as.numeric(c(fit$coefficients$cond, fit$coefficients$zi)), "params")
				return(as.numeric(fit$coefficients$cond[2]))
			} else {
				dat = if (length(X_cols) == 0L) data.frame(y = private$y, w = private$w) else data.frame(y = private$y, w = private$w, X_cov)
				mod = private$fit_zero_augmented_model(dat)
				if (is.null(mod)) return(NA_real_)
				
				cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (is.null(cond_coef) || !("w" %in% names(cond_coef))) return(NA_real_)
				return(as.numeric(cond_coef["w"]))
			}
		},

		za_family = function() stop(class(self)[1], " must implement za_family()."),

		za_description = function() stop(class(self)[1], " must implement za_description()."),

		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp) && identical(private$za_description(), "Zero-Inflated Negative Binomial")
		},

		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					start_params = start %||% private$get_fit_warm_start_for_length("params", 2L * ncol(X_fit) + 1L)
					fast_zinb_cpp(
						X = X_fit,
						y = y,
						Xzi = X_fit,
						start_params = start_params,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg,
						fixed_idx = j_treat,
						fixed_values = delta
					)
				},
				extract_start = function(fit){
					as.numeric(fit$params)
				},
				score = function(fit){
					as.numeric(fit$score %||% get_zinb_score_cpp(X = X_fit, y = y, Xzi = X_fit, as.numeric(fit$params)))
				},
				observed_information = function(fit){
					as.matrix(fit$information)
				},
				information = function(fit){
					as.matrix(fit$information)
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_zinb_neg_loglik_cpp(X = X_fit, y = y, Xzi = X_fit, as.numeric(fit$params)))
				}
			)
		},

		predictors_df = function(){
			if (ncol(as.matrix(private$X)) > 0){
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
			private$cached_values$likelihood_test_context = NULL

			if (is.null(private$best_X_colnames)) {
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
				private$best_X_colnames = setdiff(colnames(X_reduced), c("(Intercept)", "w"))
			} else {
				X_data = private$get_X()
				X_cols = private$best_X_colnames
				if (length(X_cols) == 0L){
					X_fit = cbind(1, w = private$w)
				} else {
					X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
					X_fit = cbind(1, w = private$w, X_cov)
				}
			}

			if (private$use_rcpp && identical(private$za_description(), "Zero-Inflated Negative Binomial")) {
				start_params = private$get_fit_warm_start_for_length("params", 2L * ncol(X_fit) + 1L)
				fit = tryCatch(
					fast_zinb_cpp(X = X_fit, y = as.numeric(private$y), Xzi = X_fit, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) {
					private$cache_nonestimable_estimate("zinb_fit_unavailable")
					return(invisible(NULL))
				}
				private$set_fit_warm_start(as.numeric(fit$params), "params")
				private$cached_mod = fit
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					j_treat = 2L
				)
				private$cached_values$beta_hat_T = as.numeric(fit$coefficients$cond[2])
				if (!estimate_only) {
					se = tryCatch(sqrt(fit$vcov[2, 2]), error = function(e) NA_real_)
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				start_params = private$get_fit_warm_start_for_length("params", 2L * ncol(X_fit))
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(X = X_fit, y = as.numeric(private$y), Xzi = X_fit, is_hurdle = is_hurdle, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) {
					private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
					return(invisible(NULL))
				}
				private$set_fit_warm_start(as.numeric(c(fit$coefficients$cond, fit$coefficients$zi)), "params")
				private$cached_mod = fit
				private$cached_values$likelihood_test_context = NULL

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
					private$cached_values$likelihood_test_context = NULL

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
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)
)
