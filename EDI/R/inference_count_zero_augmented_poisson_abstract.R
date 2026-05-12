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
	inherit = InferenceAsympLikStdModCache,
	public = list(
				
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param model_formula_zero Formula for the zero/hurdle auxiliary component. Defaults
		#'   to \code{~ .}, meaning treatment plus all available covariates are used in that
		#'   auxiliary submodel.
		#' @param use_rcpp Whether to use Rcpp speedup.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg  Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, model_formula_zero = ~ ., use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFormula(model_formula, null.ok = TRUE)
				assertFormula(model_formula_zero, null.ok = FALSE)
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
			private$model_formula_zero = model_formula_zero
		},
		#' @description Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute asymp confidence interval
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
		#' @description Compute asymp two sided pval for treatment effect
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
		cached_mod = NULL,
		get_standard_error = function(){
			private$shared(estimate_only = FALSE)
			se = private$compute_standard_error_from_information_matrix()
			if (is.finite(se)) return(se)
			private$cached_values$s_beta_hat_T
		},
		get_degrees_of_freedom = function(){
			private$shared(estimate_only = FALSE)
			private$cached_values$df %||% NA_real_
		},
		best_X_colnames = NULL,
		best_Xzi_colnames = NULL,
		use_rcpp = TRUE,
		model_formula_zero = NULL,
		build_component_matrix = function(model_formula, selected_colnames = NULL, treatment_name = "treatment"){
			if (is.null(selected_colnames)) {
				if (identical(model_formula, ~ .)) {
					X_cov = private$get_X()
				} else {
					X_imp = private$des_obj$get_X_imp()
					if (is.null(X_imp)) {
						X_cov = matrix(NA_real_, nrow = private$n, ncol = 0)
					} else {
						X_cov = create_model_matrix_from_features(model_formula, X_imp)
					}
				}
			} else {
				if (identical(model_formula, ~ .)) {
					X_cov_all = private$get_X()
				} else {
					X_imp = private$des_obj$get_X_imp()
					X_cov_all = if (is.null(X_imp)) matrix(NA_real_, nrow = private$n, ncol = 0) else create_model_matrix_from_features(model_formula, X_imp)
				}
				if (is.null(X_cov_all) || length(selected_colnames) == 0L) {
					X_cov = matrix(NA_real_, nrow = private$n, ncol = 0)
				} else {
					X_cov = as.matrix(X_cov_all[, intersect(selected_colnames, colnames(X_cov_all)), drop = FALSE])
				}
			}
			if (is.null(X_cov) || ncol(as.matrix(X_cov)) == 0L) {
				X_fit = cbind(1, private$w)
				colnames(X_fit) = c("(Intercept)", treatment_name)
				return(X_fit)
			}
			X_cov = as.matrix(X_cov)
			if (isTRUE(private$harden)) {
				X_cov = drop_highly_correlated_cols(X_cov, threshold = 0.999)$M
			}
			X_fit = cbind(1, private$w, X_cov)
			colnames(X_fit)[1:2] = c("(Intercept)", treatment_name)
			if (isTRUE(private$harden)) {
				res = drop_linearly_dependent_cols(X_fit)
				X_fit = res$M
				colnames(X_fit) = c("(Intercept)", treatment_name, colnames(X_cov))[res$js]
			}
			X_fit
		},
		build_component_frame = function(X_cond, Xzi){
			dat = data.frame(y = private$y, w = private$w)
			if (ncol(X_cond) > 2L) {
				Xc = as.data.frame(X_cond[, -c(1, 2), drop = FALSE])
				names(Xc) = make.names(colnames(X_cond)[-c(1, 2)], unique = TRUE)
				dat = cbind(dat, Xc)
			}
			if (ncol(Xzi) > 2L) {
				Xz = as.data.frame(Xzi[, -c(1, 2), drop = FALSE])
				names(Xz) = make.names(colnames(Xzi)[-c(1, 2)], unique = TRUE)
				for (nm in names(Xz)) {
					if (!nm %in% names(dat)) dat[[nm]] = Xz[[nm]]
				}
			}
			dat
		},
		build_formula_from_matrix = function(X_fit, response = "y"){
			vars = colnames(X_fit)[-1]
			rhs = if (!length(vars)) "1" else paste(vars, collapse = " + ")
			if (is.null(response)) return(stats::as.formula(paste("~", rhs)))
			stats::as.formula(paste(response, "~", rhs))
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_fit = private$build_component_matrix(private$model_formula, private$best_X_colnames, treatment_name = "treatment")
			Xzi_fit = private$build_component_matrix(private$model_formula_zero, private$best_Xzi_colnames, treatment_name = "treatment")
			if (private$use_rcpp && identical(private$za_description(), "Zero-Inflated Negative Binomial")) {
				start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + ncol(Xzi_fit) + 1L)
				fit = tryCatch(
					fast_zinb_cpp(X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				private$set_fit_warm_start(as.numeric(fit$params), "params")
				return(as.numeric(fit$coefficients$cond[2]))
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + ncol(Xzi_fit))
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit, is_hurdle = is_hurdle, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				private$set_fit_warm_start(as.numeric(c(fit$coefficients$cond, fit$coefficients$zi)), "params")
				return(as.numeric(fit$coefficients$cond[2]))
			} else {
				dat = private$build_component_frame(X_fit, Xzi_fit)
				mod = private$fit_zero_augmented_model(dat, X_fit, Xzi_fit)
				if (is.null(mod)) return(NA_real_)
				
				cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (is.null(cond_coef) || !("w" %in% names(cond_coef))) return(NA_real_)
				return(as.numeric(cond_coef["w"]))
			}
		},
		za_family = function() stop(class(self)[1], " must implement za_family()."),
		za_description = function() stop(class(self)[1], " must implement za_description()."),
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp) && private$za_description() %in% c(
				"Zero-Inflated Negative Binomial",
				"Zero-Inflated Poisson",
				"Hurdle Poisson"
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			Xzi_fit = ctx$Xzi
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			is_hurdle = isTRUE(ctx$is_hurdle)
			is_zinb = identical(private$za_description(), "Zero-Inflated Negative Binomial")
			start_len = if (is_zinb) ncol(X_fit) + ncol(Xzi_fit) + 1L else ncol(X_fit) + ncol(Xzi_fit)
			list(
				X = X_fit,
				Xzi = Xzi_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					start_params = start %||% private$get_fit_warm_start_for_length("params", start_len)
					fit = if (is_zinb) {
						fast_zinb_cpp(
							X = X_fit,
							y = y,
							Xzi = Xzi_fit,
							start_params = start_params,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					} else {
						fast_zero_augmented_poisson_cpp(
							X = X_fit,
							y = y,
							Xzi = Xzi_fit,
							is_hurdle = is_hurdle,
							start_params = start_params,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					}
					if (!is.null(fit) && is.null(fit$b) && is.null(fit$params) && !is.null(fit$coefficients)) {
						fit$params = as.numeric(c(fit$coefficients$cond, fit$coefficients$zi))
					}
					fit
				},
				extract_start = function(fit){
					as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
				},
				score = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb) {
						as.numeric(fit$score %||% get_zinb_score_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params))
					} else {
						as.numeric(fit$score %||% get_zero_augmented_poisson_score_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params, is_hurdle = is_hurdle))
					}
				},
				observed_information = function(fit){
					mat = fit$information %||% fit$observed_information
					if (!is.null(mat)) return(as.matrix(mat))
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb) NULL
					else as.matrix(-get_zero_augmented_poisson_hessian_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params = params, is_hurdle = is_hurdle))
				},
				information = function(fit){
					mat = fit$information %||% fit$observed_information
					if (!is.null(mat)) return(as.matrix(mat))
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb) NULL
					else as.matrix(-get_zero_augmented_poisson_hessian_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params = params, is_hurdle = is_hurdle))
				},
				neg_loglik = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb) {
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_zinb_neg_loglik_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params))
					} else {
						as.numeric(fit$neg_loglik %||% fit$neg_ll)
					}
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
		fit_zero_augmented_model = function(dat, X_fit, Xzi_fit){
			formula_cond = private$build_formula_from_matrix(X_fit)
			formula_zi = private$build_formula_from_matrix(Xzi_fit, response = NULL)
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
				X_full = private$build_component_matrix(private$model_formula, treatment_name = "w")
				res_reduced = private$reduce_design_matrix_preserving_treatment(X_full)
				X_fit = res_reduced$X
				if (is.null(X_fit)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_design_unusable")
					return(invisible(NULL))
				}
				colnames(X_fit) = colnames(X_full)[res_reduced$keep]
				private$best_X_colnames = setdiff(colnames(X_fit), c("(Intercept)", "w"))
			} else {
				X_fit = private$build_component_matrix(private$model_formula, private$best_X_colnames, treatment_name = "w")
			}
			if (is.null(private$best_Xzi_colnames)) {
				Xzi_full = private$build_component_matrix(private$model_formula_zero, treatment_name = "w")
				res_zi = private$reduce_design_matrix_preserving_treatment(Xzi_full)
				Xzi_fit = res_zi$X
				if (is.null(Xzi_fit)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_aux_design_unusable")
					return(invisible(NULL))
				}
				colnames(Xzi_fit) = colnames(Xzi_full)[res_zi$keep]
				private$best_Xzi_colnames = setdiff(colnames(Xzi_fit), c("(Intercept)", "w"))
			} else {
				Xzi_fit = private$build_component_matrix(private$model_formula_zero, private$best_Xzi_colnames, treatment_name = "w")
			}
			if (private$use_rcpp && identical(private$za_description(), "Zero-Inflated Negative Binomial")) {
				start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + ncol(Xzi_fit) + 1L)
				fit = tryCatch(
					fast_zinb_cpp(X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
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
					Xzi = Xzi_fit,
					j_treat = 2L,
					is_hurdle = FALSE
				)
				private$cached_values$beta_hat_T = as.numeric(fit$coefficients$cond[2])
				if (!estimate_only) {
					se = tryCatch(sqrt(fit$vcov[2, 2]), error = function(e) NA_real_)
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + ncol(Xzi_fit))
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit, is_hurdle = is_hurdle, start_params = start_params, estimate_only = estimate_only, optimization_alg = private$optimization_alg),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) {
					private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
					return(invisible(NULL))
				}
				private$set_fit_warm_start(as.numeric(c(fit$coefficients$cond, fit$coefficients$zi)), "params")
				private$cached_mod = fit
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					Xzi = Xzi_fit,
					j_treat = 2L,
					is_hurdle = is_hurdle
				)
				private$cached_values$beta_hat_T = as.numeric(fit$coefficients$cond[2])
				if (!estimate_only) {
					se = tryCatch(sqrt(fit$vcov[2, 2]), error = function(e) NA_real_)
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
			} else {
				dat = private$build_component_frame(X_fit, Xzi_fit)
				mod = private$fit_zero_augmented_model(dat, X_fit, Xzi_fit)
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
