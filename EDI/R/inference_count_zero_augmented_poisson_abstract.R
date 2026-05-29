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
	inherit = InferenceCountLikelihood,
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
		initialize = function(des_obj, model_formula = NULL, model_formula_zero = ~ ., use_rcpp = TRUE, verbose = FALSE, smart_cold_start_default = NULL, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFormula(model_formula, null.ok = TRUE)
				assertFormula(model_formula_zero, null.ok = FALSE)
				assertFlag(use_rcpp)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
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
		#' @description Compute asymp confidence interval
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			se = private$get_standard_error()
			if (is.finite(se) && se > 0) {
				private$cached_values$s_beta_hat_T = se
			}
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
			se = private$get_standard_error()
			if (is.finite(se) && se > 0) {
				private$cached_values$s_beta_hat_T = se
			}
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Computes the treatment effect estimate for a weighted bootstrap sample.
		#' @param subject_or_block_weights Bootstrap weights at the subject or block level.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			if (!check_package_installed("glmmTMB")) {
				stop(class(self)[1], " weighted bootstrap estimation requires package 'glmmTMB'.")
			}
			row_weights = as.numeric(private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights))
			X_fit = private$build_component_matrix(private$model_formula, private$best_X_colnames, treatment_name = "w")
			if (is.null(X_fit)) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_design_unusable")
				return(NA_real_)
			}
			Xzi_fit = private$build_component_matrix(private$model_formula_zero, private$best_Xzi_colnames, treatment_name = "w")
			if (is.null(Xzi_fit)) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_aux_design_unusable")
				return(NA_real_)
			}
			dat = private$build_component_frame(X_fit, Xzi_fit)
			mod = private$fit_zero_augmented_model(dat, X_fit, Xzi_fit, weights = row_weights)
			if (is.null(mod)) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_weighted_fit_unavailable")
				return(NA_real_)
			}
			cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
			if (is.null(cond_coef) || !("w" %in% names(cond_coef)) || !is.finite(cond_coef["w"])) {
				private$cache_nonestimable_estimate("zero_augmented_poisson_weighted_treatment_missing")
				return(NA_real_)
			}
			private$cached_mod = mod
			private$cached_values$likelihood_test_context = NULL
			private$cached_values$beta_hat_T = as.numeric(cond_coef["w"])
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$full_coefficients = cond_coef
			private$cached_values$beta_hat_T
		}
	),
		private = list(
		supports_reusable_bootstrap_worker = function(){
			FALSE
		},
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
		get_complexity_tier = function() "heavy",
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
				n_params = ncol(X_fit) + ncol(Xzi_fit) + 1L
				ws_args = private$get_backend_warm_start_args(n_params)
				fit = tryCatch(
					fast_zinb_cpp(
						X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				private$set_fit_warm_start(as.numeric(fit$params), "params", fisher = fit$fisher_information)
				return(as.numeric(fit$params[2]))
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				n_params = ncol(X_fit) + ncol(Xzi_fit)
				ws_args = private$get_backend_warm_start_args(n_params)
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(
						X_fit, as.numeric(private$y), Xzi_fit,
						is_hurdle = is_hurdle,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
				private$set_fit_warm_start(as.numeric(c(fit$coefficients$cond, fit$coefficients$zi)), "params", fisher = fit$fisher_information)
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
					warm_start_params = start %||% private$get_fit_warm_start_for_length("params", start_len)
					warm_fisher = private$get_fit_warm_start_fisher(start_len)
					fit = if (is_zinb) {
						fast_zinb_cpp(
							X = X_fit,
							y = y,
							Xzi = Xzi_fit,
							warm_start_params = warm_start_params,
							warm_start_fisher_info = warm_fisher,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					} else {
						fast_zero_augmented_poisson_cpp(
							X_fit,
							y,
							Xzi_fit,
							is_hurdle = is_hurdle,
							warm_start_params = warm_start_params,
							warm_start_fisher_info = warm_fisher,
							smart_cold_start = private$smart_cold_start_default,
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
						as.numeric(fit$score %||% get_zero_augmented_poisson_score_cpp(X_fit, y, Xzi_fit, params, is_hurdle = is_hurdle))
					}
				},
				observed_information = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					as.matrix(fit$observed_information %||% (if (is_zinb) -get_zinb_hessian_cpp(X_fit, y, Xzi_fit, params) else -get_zero_augmented_poisson_hessian_cpp(X_fit, y, Xzi_fit, params, is_hurdle)))
				},
				information = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					as.matrix(fit$information %||% fit$fisher_information %||% fit$observed_information %||% (if (is_zinb) -get_zinb_hessian_cpp(X_fit, y, Xzi_fit, params) else -get_zero_augmented_poisson_hessian_cpp(X_fit, y, Xzi_fit, params, is_hurdle)))
				},
				neg_loglik = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb) {
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_zinb_neg_loglik_cpp(X = X_fit, y = y, Xzi = Xzi_fit, params))
					} else {
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% fit$mod$neg_ll %||% fit$mod$neg_loglik)
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
		fit_zero_augmented_model = function(dat, X_fit, Xzi_fit, weights = NULL){
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
						weights = weights,
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
						weights = weights,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
		},
		generate_mod = function(estimate_only = FALSE){
			private$cached_values$likelihood_test_context = NULL
			if (is.null(private$best_X_colnames)) {
				X_full = private$build_component_matrix(private$model_formula, treatment_name = "w")
				res_reduced = private$reduce_design_matrix_preserving_treatment(X_full)
				X_fit = res_reduced$X
				if (is.null(X_fit)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_design_unusable")
					return(NULL)
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
					return(NULL)
				}
				colnames(Xzi_fit) = colnames(Xzi_full)[res_zi$keep]
				private$best_Xzi_colnames = setdiff(colnames(Xzi_fit), c("(Intercept)", "w"))
			} else {
				Xzi_fit = private$build_component_matrix(private$model_formula_zero, private$best_Xzi_colnames, treatment_name = "w")
			}
			
			out = list()
			if (private$use_rcpp && identical(private$za_description(), "Zero-Inflated Negative Binomial")) {
				n_params = ncol(X_fit) + ncol(Xzi_fit) + 1L
				ws_args = private$get_backend_warm_start_args(n_params)
				fit = tryCatch(
					fast_zinb_cpp(
						X = X_fit, y = as.numeric(private$y), Xzi = Xzi_fit,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) {
					private$cache_nonestimable_estimate("zinb_fit_unavailable")
					return(NULL)
				}
				
				private$cached_mod = fit
				private$set_fit_warm_start(as.numeric(fit$params), "params", fisher = fit$fisher_information)
				
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					Xzi = Xzi_fit,
					j_treat = 2L,
					is_hurdle = FALSE
				)
				out$beta_hat_T = as.numeric(fit$params[2])
				if (!estimate_only) {
					se = tryCatch(sqrt(fit$vcov[2, 2]), error = function(e) NA_real_)
					out$ssq_b_j = if (is.finite(se) && se > 0) se^2 else NA_real_
				}
				out$params = as.numeric(fit$params)
				out$fisher_information = fit$fisher_information
				out$mod = fit
			} else if (private$use_rcpp && !grepl("Negative Binomial", private$za_description())) {
				is_hurdle = identical(private$za_description(), "Hurdle Poisson")
				n_params = ncol(X_fit) + ncol(Xzi_fit)
				ws_args = private$get_backend_warm_start_args(n_params)
				fit = tryCatch(
					fast_zero_augmented_poisson_cpp(
						X_fit, as.numeric(private$y), Xzi_fit,
						is_hurdle = is_hurdle,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = estimate_only, optimization_alg = private$optimization_alg
					),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged)) {
					private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
					return(NULL)
				}
				
				private$cached_mod = fit
				full_params = as.numeric(c(fit$coefficients$cond, fit$coefficients$zi))
				private$set_fit_warm_start(full_params, "params", fisher = fit$fisher_information)
				
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					Xzi = Xzi_fit,
					j_treat = 2L,
					is_hurdle = is_hurdle
				)
				out$beta_hat_T = as.numeric(fit$coefficients$cond[2])
				if (!estimate_only) {
					se = tryCatch(sqrt(fit$vcov[2, 2]), error = function(e) NA_real_)
					out$ssq_b_j = if (is.finite(se) && se > 0) se^2 else NA_real_
				}
				out$params = full_params
				out$fisher_information = fit$fisher_information
				out$mod = fit
			} else {
				dat = private$build_component_frame(X_fit, Xzi_fit)
				mod = private$fit_zero_augmented_model(dat, X_fit, Xzi_fit)
				if (is.null(mod)){
					private$cache_nonestimable_estimate("zero_augmented_poisson_fit_unavailable")
					return(NULL)
				}
				
				private$cached_values$likelihood_test_context = NULL
				cond_coef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
				if (is.null(cond_coef) || !("w" %in% names(cond_coef))){
					private$cache_nonestimable_estimate("zero_augmented_poisson_treatment_missing")
					return(NULL)
				}
				out$beta_hat_T = as.numeric(cond_coef["w"])
				if (!estimate_only) {
					coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
					se = if (!is.null(coef_table) && ("w" %in% rownames(coef_table))) as.numeric(coef_table["w", "Std. Error"]) else NA_real_
					out$ssq_b_j = if (is.finite(se) && se > 0) se^2 else NA_real_
				}
				out$mod = mod
			}
			out
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		supports_lik_ratio_param_bootstrap = function(){
			isTRUE(private$use_rcpp) && private$za_description() %in% c(
				"Zero-Inflated Negative Binomial",
				"Zero-Inflated Poisson",
				"Hurdle Poisson"
			)
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			is_zinb   = identical(private$za_description(), "Zero-Inflated Negative Binomial")
			is_hurdle = identical(private$za_description(), "Hurdle Poisson")
			X   = spec$X
			Xzi = spec$Xzi
			j   = spec$j
			n   = nrow(X)
			if (is_zinb) {
				p_null = as.numeric(null_fit$params)
				b_cond = p_null[seq_len(ncol(X))]
				b_zi   = p_null[ncol(X) + seq_len(ncol(Xzi))]
			} else {
				b_cond = as.numeric(null_fit$coefficients$cond)
				b_zi   = as.numeric(null_fit$coefficients$zi)
			}
			lambda = exp(pmin(as.numeric(X %*% b_cond), 20))
			pi     = plogis(as.numeric(Xzi %*% b_zi))
			if (is_zinb){
				full_params = as.numeric(null_fit$params)
				theta = exp(min(full_params[length(full_params)], 15))
				if (!is.finite(theta) || theta <= 0) return(NULL)
				u = rbinom(n, 1L, pi)
				counts = rnbinom(n, size = theta, mu = lambda)
				y_sim = as.integer(ifelse(u == 1L, 0L, counts))
			} else if (is_hurdle){
				p0 = ppois(0, lambda)
				u  = rbinom(n, 1L, pi)
				y_pos = as.integer(qpois(p0 + (1 - p0) * runif(n), lambda))
				y_pos = pmax(y_pos, 1L)
				y_sim = as.integer(ifelse(u == 0L, 0L, y_pos))
			} else {
				u = rbinom(n, 1L, pi)
				y_sim = as.integer(ifelse(u == 1L, 0L, rpois(n, lambda)))
			}
			n_params = if (is_zinb) ncol(X) + ncol(Xzi) + 1L else ncol(X) + ncol(Xzi)
			full_res = tryCatch(
				if (is_zinb) {
					fast_zinb_cpp(
						X = X, y = as.numeric(y_sim), Xzi = Xzi,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg
					)
				} else {
					fast_zero_augmented_poisson_cpp(
						X, as.numeric(y_sim), Xzi,
						is_hurdle = is_hurdle,
						smart_cold_start = private$smart_cold_start_default,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg
					)
				},
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged)) return(NULL)
			beta_j = if (is_zinb) as.numeric(full_res$params[j]) else as.numeric(full_res$coefficients$cond[j])
			if (!is.finite(beta_j)) return(NULL)
			if (is.null(full_res$params) && !is.null(full_res$coefficients)){
				full_res$params = as.numeric(c(full_res$coefficients$cond, full_res$coefficients$zi))
			}
			list(
				worker_data = list(y = as.numeric(y_sim)),
				full_fit = full_res,
				fit_null = function(d, start = NULL){
					ws = start %||% private$get_fit_warm_start_for_length("params", n_params)
					fit = tryCatch(
						if (is_zinb) {
							fast_zinb_cpp(
								X = X, y = as.numeric(y_sim), Xzi = Xzi,
								warm_start_params = ws,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = FALSE,
								optimization_alg = private$optimization_alg,
								fixed_idx = j, fixed_values = d
							)
						} else {
							fast_zero_augmented_poisson_cpp(
								X, as.numeric(y_sim), Xzi,
								is_hurdle = is_hurdle,
								warm_start_params = ws,
								smart_cold_start = private$smart_cold_start_default,
								estimate_only = FALSE,
								optimization_alg = private$optimization_alg,
								fixed_idx = j, fixed_values = d
							)
						},
						error = function(e) NULL
					)
					if (!is.null(fit) && is.null(fit$params) && !is.null(fit$coefficients)){
						fit$params = as.numeric(c(fit$coefficients$cond, fit$coefficients$zi))
					}
					fit
				},
				neg_loglik = function(fit){
					params = as.numeric(fit$params %||% c(as.numeric(fit$coefficients$cond), as.numeric(fit$coefficients$zi)))
					if (is_zinb){
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_zinb_neg_loglik_cpp(X = X, y = as.numeric(y_sim), Xzi = Xzi, params))
					} else {
						as.numeric(fit$neg_loglik %||% fit$neg_ll)
					}
				}
			)
		}
	)
)
