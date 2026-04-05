#' Marginal Standardization / G-Computation for Proportion Responses
#'
#' Internal base class for proportion-outcome g-computation estimators. A
#' fractional logit working model is fit, then potential-outcome mean
#' proportions under all-treated and all-control assignments are standardized
#' over the empirical covariate distribution. Inference uses a sandwich-robust
#' covariance for the regression coefficients and the delta method for the
#' marginal mean-difference estimand.
#'
#' @details
#' The implementation is optimized for resampling-based inference. It utilizes a
#' fast C++ IRLS solver for the underlying fractional logit regression. During
#' resampling draws, it bypasses the calculation of the sandwich covariance
#' matrix and delta-method standard errors, providing a significant speedup when
#' computing bootstrap or randomization distributions.
#'
#' @keywords internal
#' @noRd
InferencePropGCompAbstract = R6::R6Class("InferencePropGCompAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the g-computation inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with a proportion response.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "proportion")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the g-computation treatment-effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = TRUE)
			private$cached_values$md
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			tryCatch({
				private$shared(estimate_only = FALSE)
				private$compute_effect_confidence_interval(alpha)
			}, error = function(e){
				warning("G-computation mean difference: ", conditionMessage(e))
				c(NA_real_, NA_real_)
			})
		},

		#' @description
		#' Computes a two-sided p-value for the treatment effect.
		#' @param delta The null mean difference. Defaults to 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta, len = 1)
			tryCatch({
				private$shared(estimate_only = FALSE)
				private$compute_effect_pvalue(delta)
			}, error = function(e){
				warning("G-computation mean difference: ", conditionMessage(e))
				NA_real_
			})
		},

		#' @description
		#' Computes a bootstrap two-sided p-value for the treatment effect.
		#' @param delta The null mean difference. Defaults to 0.
		#' @param B Number of bootstrap samples.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param boundary_tol Resample screening threshold for boundary mass near 0/1.
		#' @param max_boundary_mass Reject a resample when at least this fraction is near the boundary.
		#' @param sep_tol Separation tolerance used to reject nearly perfectly separated resamples.
		#' @param min_group_n Minimum number of observations required in each treatment arm.
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE,
			boundary_tol = 0.02, max_boundary_mass = 0.95, sep_tol = 0.02, min_group_n = 5L){
			assertNumeric(delta, len = 1)
			old_bootstrap_screening = private$bootstrap_screening_control
			private$bootstrap_screening_control = list(
				boundary_tol = boundary_tol,
				max_boundary_mass = max_boundary_mass,
				sep_tol = sep_tol,
				min_group_n = as.integer(min_group_n)
			)
			on.exit({private$bootstrap_screening_control = old_bootstrap_screening}, add = TRUE)
			super$compute_bootstrap_two_sided_pval(delta = delta, B = B, na.rm = na.rm)
		},

		#' @description
		#' Computes a bootstrap confidence interval.
		#' @param alpha The confidence level 1 - \code{alpha}.
		#' @param B Number of bootstrap samples.
		#' @param type Bootstrap CI type.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		#' @param show_progress Whether to show bootstrap progress.
		#' @param boundary_tol Resample screening threshold for boundary mass near 0/1.
		#' @param max_boundary_mass Reject a resample when at least this fraction is near the boundary.
		#' @param sep_tol Separation tolerance used to reject nearly perfectly separated resamples.
		#' @param min_group_n Minimum number of observations required in each treatment arm.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = "percentile",
			na.rm = TRUE, show_progress = TRUE, boundary_tol = 0.02, max_boundary_mass = 0.95,
			sep_tol = 0.02, min_group_n = 5L){
			old_bootstrap_screening = private$bootstrap_screening_control
			private$bootstrap_screening_control = list(
				boundary_tol = boundary_tol,
				max_boundary_mass = max_boundary_mass,
				sep_tol = sep_tol,
				min_group_n = as.integer(min_group_n)
			)
			on.exit({private$bootstrap_screening_control = old_bootstrap_screening}, add = TRUE)
			super$compute_bootstrap_confidence_interval(
				alpha = alpha, B = B, type = type, na.rm = na.rm, show_progress = show_progress
			)
		},

		#' @description
		#' Abbreviated bootstrap sampler that reuses the reduced column layout.
		#' @param B Description for B
		#' @param show_progress Description for show_progress
		#' @param max_resample_attempts Description for max_resample_attempts
		#' @param boundary_tol Resample screening threshold for boundary mass near 0/1.
		#' @param max_boundary_mass Reject a resample when at least this fraction is near the boundary.
		#' @param sep_tol Separation tolerance used to reject nearly perfectly separated resamples.
		#' @param min_group_n Minimum number of observations required in each treatment arm.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, max_resample_attempts = 50,
			boundary_tol = 0.02, max_boundary_mass = 0.95, sep_tol = 0.02, min_group_n = 5L, debug = FALSE){
			assertCount(B, positive = TRUE)
			assertCount(max_resample_attempts, positive = TRUE)
			assertNumber(boundary_tol, lower = 0, upper = 0.5)
			assertNumber(max_boundary_mass, lower = 0, upper = 1)
			assertNumber(sep_tol, lower = 0)
			assertCount(min_group_n, positive = TRUE)
			private$shared(estimate_only = TRUE)
			cols = private$gcomp_design_colnames
			if (is.null(cols)){
				return(super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug = debug))
			}

			X_full = private$build_named_design_matrix()
			n = nrow(X_full)
			y = private$y
			w = private$w

			#' @description
			#' Draw sample
			#' @param ... Other arguments passed to the method.
				draw_sample = function(...){
					attempt = 1
					repeat {
						i_b = sample_int_replace_cpp(n, n)
						w_b = w[i_b]
						y_b = y[i_b]
						if (private$bootstrap_sample_is_usable(
							w_b, y_b,
							boundary_tol = boundary_tol,
							max_boundary_mass = max_boundary_mass,
							sep_tol = sep_tol,
							min_group_n = min_group_n
						)) break
						attempt = attempt + 1
						if (attempt > max_resample_attempts) return(NA_real_)
					}
					X_b = X_full[i_b, , drop = FALSE]
					private$bootstrap_effect_from_sample(X_b, y_b)
				}

			if (isTRUE(debug)) {
				debug_results = vector("list", B)
				for (b in seq_len(B)) {
					iter_warns = character(0)
					iter_val = withCallingHandlers(
						tryCatch(draw_sample(), error = function(e) list(val = NA_real_, error = conditionMessage(e))),
						warning = function(wrn) { iter_warns <<- c(iter_warns, conditionMessage(wrn)); invokeRestart("muffleWarning") }
					)
					debug_results[[b]] = list(
						val = if (is.list(iter_val) && !is.null(iter_val$val)) as.numeric(iter_val$val)[1L] else as.numeric(iter_val)[1L],
						errors = if (is.list(iter_val) && !is.null(iter_val$error)) iter_val$error else character(0),
						warnings = iter_warns
					)
				}
				values = sapply(debug_results, `[[`, "val")
				errors_list = lapply(debug_results, `[[`, "errors")
				warnings_list = lapply(debug_results, `[[`, "warnings")
				num_errors_vec = lengths(errors_list)
				num_warnings_vec = lengths(warnings_list)
				return(list(
					values = values,
					errors = errors_list,
					warnings = warnings_list,
					num_errors = num_errors_vec,
					num_warnings = num_warnings_vec,
					prop_iterations_with_errors = mean(num_errors_vec > 0),
					prop_iterations_with_warnings = mean(num_warnings_vec > 0),
					prop_illegal_values = mean(!is.finite(values))
				))
			} else if (self$num_cores == 1){
				vapply(seq_len(B), function(b) draw_sample(), numeric(1))
			} else {
				unlist(private$par_lapply(seq_len(B), function(b) draw_sample(), 
				n_cores = min(2L, self$num_cores),
					export_list = list(draw_sample = draw_sample, n = n, y = y, w = w, X_full = X_full, sample_int_replace_cpp = EDI:::sample_int_replace_cpp)))
			}
		}
	),

	private = list(
		gcomp_design_colnames = NULL,
		gcomp_design_j_treat = NULL,
			bootstrap_screening_control = NULL,
		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),
		build_named_design_matrix = function(){
			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			colnames(X_full) = c(
				"(Intercept)",
				"treatment",
				if (ncol(X_full) > 2L) private$get_covariate_names() else NULL
			)
			X_full
		},

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

		set_failed_fit_cache = function(){
			private$cached_values$summary_table = NULL
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$mean1 = NA_real_
			private$cached_values$mean0 = NA_real_
			private$cached_values$md = NA_real_
			private$cached_values$se_md = NA_real_
		},

		effects_are_usable = function(effects, estimate_only = FALSE){
			if (estimate_only) return(is.finite(effects$md))
			is.finite(effects$md) && is.finite(effects$se_md) && effects$se_md > 0
		},

		select_covariate_to_drop = function(X_curr, coef_hat){
			covariate_cols = seq.int(3L, ncol(X_curr))
			if (length(covariate_cols) == 0L) return(NA_integer_)

			coef_mags = abs(coef_hat[covariate_cols])
			if (length(coef_mags) == 0L || all(!is.finite(coef_mags))){
				return(tail(covariate_cols, 1L))
			}
			covariate_cols[which.max(replace(coef_mags, !is.finite(coef_mags), -Inf))]
		},

		fit_fractional_logit_with_sandwich = function(X_full, estimate_only = FALSE){
			X_curr = X_full

			repeat {
				reduced = private$reduce_design_matrix_preserving_treatment(X_curr)
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
					return(NULL)
				}

				mod = tryCatch(
					fast_logistic_regression_cpp(X = X_fit, y = as.numeric(private$y)),
					error = function(e) NULL
				)
				if (is.null(mod)){
					return(NULL)
				}

				coef_hat = as.numeric(mod$b)
				converged = all(is.finite(coef_hat))
				if (!converged){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				if (estimate_only){
					private$gcomp_design_colnames = colnames(X_fit)
					private$gcomp_design_j_treat = j_treat
					return(list(
						X = X_fit,
						j_treat = j_treat,
						coefficients = coef_hat,
						estimate_only = TRUE
					))
				}

				mu_hat = inv_logit(X_fit %*% coef_hat)
				mu_hat = pmin(pmax(as.numeric(mu_hat), .Machine$double.eps), 1 - .Machine$double.eps)
				W = mu_hat * (1 - mu_hat)
				if (any(!is.finite(W)) || any(W <= 0)){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				post_fit = tryCatch(
					gcomp_fractional_logit_post_fit_cpp(
						X_fit = X_fit,
						y = as.numeric(private$y),
						coef_hat = coef_hat,
						mu_hat = mu_hat,
						j_treat = j_treat
					),
					error = function(e) NULL
				)
				if (is.null(post_fit)){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				coef_names = colnames(X_fit)
				names(coef_hat) = coef_names
				vcov_robust = post_fit$vcov
				colnames(vcov_robust) = rownames(vcov_robust) = coef_names

				private$gcomp_design_colnames = coef_names
				private$gcomp_design_j_treat = j_treat

				return(list(
					X = X_fit,
					j_treat = j_treat,
					coefficients = coef_hat,
					vcov = vcov_robust,
					post_fit = post_fit,
					estimate_only = FALSE
				))
			}
		},

			compute_standardized_effect_components = function(X_fit, coef_hat, j_treat, clip_lower = NULL, clip_upper = NULL){
				X1 = X_fit
				X0 = X_fit
				X1[, j_treat] = 1
				X0[, j_treat] = 0

				eta1 = as.numeric(X1 %*% coef_hat)
				eta0 = as.numeric(X0 %*% coef_hat)
				mean1_i = stats::plogis(eta1)
				mean0_i = stats::plogis(eta0)
				if (!is.null(clip_lower) || !is.null(clip_upper)) {
					if (is.null(clip_lower)) clip_lower = 0
					if (is.null(clip_upper)) clip_upper = 1
					mean1_i = pmin(clip_upper, pmax(clip_lower, mean1_i))
					mean0_i = pmin(clip_upper, pmax(clip_lower, mean0_i))
				}

				list(
					X1 = X1,
					X0 = X0,
					mean1_i = mean1_i,
					mean0_i = mean0_i,
					mean1 = mean(mean1_i),
					mean0 = mean(mean0_i),
					md = mean(mean1_i) - mean(mean0_i)
				)
			},

		stabilize_covariance_matrix = function(vcov_mat, ridge_factor = 1e-8){
			if (is.null(vcov_mat) || !is.matrix(vcov_mat) || any(!is.finite(vcov_mat))) return(NULL)
			vcov_sym = (vcov_mat + t(vcov_mat)) / 2
			eig = tryCatch(eigen(vcov_sym, symmetric = TRUE), error = function(e) NULL)
			if (is.null(eig) || any(!is.finite(eig$values))) return(NULL)
			scale = max(1, max(abs(eig$values), na.rm = TRUE))
			floor_val = ridge_factor * scale
			vals = pmax(eig$values, floor_val)
			vcov_psd = eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
			vcov_psd = (vcov_psd + t(vcov_psd)) / 2
			dimnames(vcov_psd) = dimnames(vcov_mat)
			vcov_psd
		},

		invert_information_matrix = function(info_mat, ridge_factor = 1e-8){
			if (is.null(info_mat) || !is.matrix(info_mat) || any(!is.finite(info_mat))) return(NULL)
			info_sym = (info_mat + t(info_mat)) / 2
			eig = tryCatch(eigen(info_sym, symmetric = TRUE), error = function(e) NULL)
			if (is.null(eig) || any(!is.finite(eig$values))) return(NULL)
			scale = max(1, max(abs(eig$values), na.rm = TRUE))
			floor_val = ridge_factor * scale
			vals = pmax(eig$values, floor_val)
			inv_info = eig$vectors %*% diag(1 / vals, nrow = length(vals)) %*% t(eig$vectors)
			inv_info = (inv_info + t(inv_info)) / 2
			dimnames(inv_info) = dimnames(info_mat)
			inv_info
		},

			model_based_covariance_matrix = function(X_fit, coef_hat, clip_lower = 1e-6, clip_upper = 1 - 1e-6){
				eta = as.numeric(X_fit %*% coef_hat)
				mu_hat = stats::plogis(eta)
				mu_hat = pmin(clip_upper, pmax(clip_lower, mu_hat))
				W = pmin(pmax(as.numeric(mu_hat * (1 - mu_hat)), clip_lower * (1 - clip_upper)), 0.25)
				info_mat = crossprod(X_fit, X_fit * W)
				colnames(info_mat) = rownames(info_mat) = names(coef_hat)
				private$invert_information_matrix(info_mat)
			},

			finite_difference_md_gradient = function(X_fit, coef_hat, j_treat, clip_lower = NULL, clip_upper = NULL){
				p = length(coef_hat)
				grad = rep(NA_real_, p)
				base_step = .Machine$double.eps^(1/3)
				for (j in seq_len(p)) {
					h = base_step * (abs(coef_hat[j]) + 1)
					if (!is.finite(h) || h <= 0) h = 1e-6
					beta_plus = coef_hat
					beta_minus = coef_hat
					beta_plus[j] = beta_plus[j] + h
					beta_minus[j] = beta_minus[j] - h
					md_plus = private$compute_standardized_effect_components(X_fit, beta_plus, j_treat, clip_lower = clip_lower, clip_upper = clip_upper)$md
					md_minus = private$compute_standardized_effect_components(X_fit, beta_minus, j_treat, clip_lower = clip_lower, clip_upper = clip_upper)$md
					if (is.finite(md_plus) && is.finite(md_minus)) {
						grad[j] = (md_plus - md_minus) / (2 * h)
					}
				}
				names(grad) = names(coef_hat)
				grad
			},

		variance_from_gradient = function(grad, vcov_mat){
			if (is.null(vcov_mat) || any(!is.finite(grad)) || any(!is.finite(vcov_mat))) return(NA_real_)
			var_md = suppressWarnings(as.numeric(t(grad) %*% vcov_mat %*% grad))
			if (is.finite(var_md) && var_md > 0) var_md else NA_real_
		},

			resolve_effect_variance = function(X_fit, coef_hat, j_treat, vcov_robust, components){
				grad1 = as.numeric(crossprod(components$X1, components$mean1_i * (1 - components$mean1_i))) / nrow(components$X1)
				grad0 = as.numeric(crossprod(components$X0, components$mean0_i * (1 - components$mean0_i))) / nrow(components$X0)
				grad_md = grad1 - grad0

				var_md = private$variance_from_gradient(grad_md, vcov_robust)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_robust, method = "robust"))

				vcov_stable = private$stabilize_covariance_matrix(vcov_robust)
				var_md = private$variance_from_gradient(grad_md, vcov_stable)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_stable, method = "stabilized_robust"))

				vcov_model = private$model_based_covariance_matrix(X_fit, coef_hat)
				var_md = private$variance_from_gradient(grad_md, vcov_model)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_model, method = "model_based"))

				grad_fd = private$finite_difference_md_gradient(X_fit, coef_hat, j_treat)
				var_md = private$variance_from_gradient(grad_fd, vcov_stable)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_stable, method = "stabilized_robust_fd"))

				var_md = private$variance_from_gradient(grad_fd, vcov_model)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_model, method = "model_based_fd"))

				strong_clip = 1e-4
				components_clipped = private$compute_standardized_effect_components(X_fit, coef_hat, j_treat, clip_lower = strong_clip, clip_upper = 1 - strong_clip)
				grad1_clip = as.numeric(crossprod(components_clipped$X1, components_clipped$mean1_i * (1 - components_clipped$mean1_i))) / nrow(components_clipped$X1)
				grad0_clip = as.numeric(crossprod(components_clipped$X0, components_clipped$mean0_i * (1 - components_clipped$mean0_i))) / nrow(components_clipped$X0)
				grad_md_clip = grad1_clip - grad0_clip

				var_md = private$variance_from_gradient(grad_md_clip, vcov_stable)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_stable, method = "stabilized_robust_strong_clip"))

				vcov_model_clip = private$model_based_covariance_matrix(X_fit, coef_hat, clip_lower = strong_clip, clip_upper = 1 - strong_clip)
				var_md = private$variance_from_gradient(grad_md_clip, vcov_model_clip)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_model_clip, method = "model_based_strong_clip"))

				grad_fd_clip = private$finite_difference_md_gradient(X_fit, coef_hat, j_treat, clip_lower = strong_clip, clip_upper = 1 - strong_clip)
				var_md = private$variance_from_gradient(grad_fd_clip, vcov_model_clip)
				if (is.finite(var_md)) return(list(var_md = var_md, vcov = vcov_model_clip, method = "model_based_fd_strong_clip"))

				list(
					var_md = NA_real_,
					vcov = if (!is.null(vcov_model_clip)) vcov_model_clip else if (!is.null(vcov_stable)) vcov_stable else vcov_model,
					method = "failed"
				)
			},

		compute_standardized_effects_r = function(fit){
			X_fit = fit$X
			coef_hat = fit$coefficients
			vcov_robust = fit$vcov
			j_treat = fit$j_treat
			estimate_only = isTRUE(fit$estimate_only)

			components = private$compute_standardized_effect_components(X_fit, coef_hat, j_treat)

			if (estimate_only){
				return(list(
					mean1 = components$mean1,
					mean0 = components$mean0,
					md = components$md,
					se_md = NA_real_,
					full_coefficients = coef_hat,
					full_vcov = NULL,
					summary_table = NULL
				))
			}

			var_res = private$resolve_effect_variance(X_fit, coef_hat, j_treat, vcov_robust, components)
			vcov_used = var_res$vcov
			std_err = if (!is.null(vcov_used)) sqrt(pmax(diag(vcov_used), 0)) else rep(NA_real_, length(coef_hat))
			std_err = as.numeric(std_err)
			names(std_err) = names(coef_hat)
			z_vals = rep(NA_real_, length(coef_hat))
			ok_se = is.finite(std_err) & std_err > 0
			z_vals[ok_se] = coef_hat[ok_se] / std_err[ok_se]
			names(z_vals) = names(coef_hat)
			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)

			list(
				mean1 = components$mean1,
				mean0 = components$mean0,
				md = components$md,
				se_md = if (is.finite(var_res$var_md)) sqrt(var_res$var_md) else NA_real_,
				full_coefficients = coef_hat,
				full_vcov = vcov_used,
				summary_table = summary_table
			)
		},

		compute_standardized_effects = function(fit){
			if (isTRUE(fit$estimate_only)){
				return(private$compute_standardized_effects_r(fit))
			}
			coef_hat = fit$coefficients
			fast = fit$post_fit
			if (is.null(fast)){
				return(private$compute_standardized_effects_r(fit))
			}

			fit$vcov = fast$vcov
			colnames(fit$vcov) = rownames(fit$vcov) = names(coef_hat)
			private$compute_standardized_effects_r(fit)
		},

			bootstrap_sample_is_usable = function(w_b, y_b, boundary_tol = 0.02, max_boundary_mass = 0.95, sep_tol = 0.02, min_group_n = 5L){
				if (length(w_b) != length(y_b) || length(y_b) == 0L) return(FALSE)
				if (any(!is.finite(y_b))) return(FALSE)
				if (!any(w_b == 1, na.rm = TRUE) || !any(w_b == 0, na.rm = TRUE)) return(FALSE)

				boundary_mass = mean(y_b <= boundary_tol | y_b >= 1 - boundary_tol)
				if (!is.finite(boundary_mass) || boundary_mass >= max_boundary_mass) return(FALSE)

				y1 = y_b[w_b == 1]
				y0 = y_b[w_b == 0]
				if (length(y1) < min_group_n || length(y0) < min_group_n) return(FALSE)

				sep_hi = min(y1, na.rm = TRUE) >= max(y0, na.rm = TRUE) + sep_tol
				sep_lo = min(y0, na.rm = TRUE) >= max(y1, na.rm = TRUE) + sep_tol
				if (isTRUE(sep_hi) || isTRUE(sep_lo)) return(FALSE)

				TRUE
			},

			bootstrap_fit_from_sample = function(X_b_full, y_b){
				X_curr = X_b_full

				repeat {
					reduced = private$reduce_design_matrix_preserving_treatment(X_curr)
					X_fit = reduced$X
					j_treat = reduced$j_treat
					if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
						return(NULL)
					}

					mod = tryCatch(
						fast_logistic_regression_cpp(X = X_fit, y = as.numeric(y_b)),
						error = function(e) NULL
					)
					if (is.null(mod)){
						return(NULL)
					}

					coef_hat = as.numeric(mod$b)
					if (all(is.finite(coef_hat))){
						coef_names = colnames(X_fit)
						names(coef_hat) = coef_names
						return(list(X = X_fit, j_treat = j_treat, coefficients = coef_hat))
					}

					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
				}
			},

			bootstrap_effect_from_sample = function(X_b_full, y_b){
				if (nrow(X_b_full) == 0) return(NA_real_)
				fit = private$bootstrap_fit_from_sample(X_b_full, y_b)
				if (is.null(fit)) return(NA_real_)

				md_val = private$compute_standardized_effect_components(fit$X, fit$coefficients, fit$j_treat)$md
				if (!is.finite(md_val)){
					return(NA_real_)
				}
				md_val
			},

		compute_effect_confidence_interval = function(alpha){
			z = stats::qnorm(1 - alpha / 2)
			est = private$cached_values$md
			se = private$cached_values$se_md
			if (!is.finite(est) || !is.finite(se) || se <= 0){
				return(c(NA_real_, NA_real_))
			}
			ci = est + c(-1, 1) * z * se
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		compute_effect_pvalue = function(delta){
			est = private$cached_values$md
			se = private$cached_values$se_md
			if (!is.finite(est) || !is.finite(se) || se <= 0){
				return(NA_real_)
			}
			z_stat = (est - delta) / se
			2 * stats::pnorm(-abs(z_stat))
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$md) && (estimate_only || !is.null(private$cached_values$summary_table))) return(invisible(NULL))

			X_full = private$build_named_design_matrix()

			fit = private$fit_fractional_logit_with_sandwich(X_full, estimate_only = estimate_only)
			effects = if (!is.null(fit)) private$compute_standardized_effects(fit) else NULL
			if ((is.null(fit) || is.null(effects) || !private$effects_are_usable(effects, estimate_only)) && ncol(X_full) > 2L){
				fit = private$fit_fractional_logit_with_sandwich(X_full[, 1:2, drop = FALSE], estimate_only = estimate_only)
				effects = if (!is.null(fit)) private$compute_standardized_effects(fit) else NULL
			}
			if (is.null(fit) || is.null(effects) || !private$effects_are_usable(effects, estimate_only)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$summary_table = effects$summary_table
			private$cached_values$full_coefficients = effects$full_coefficients
			private$cached_values$full_vcov = effects$full_vcov
			private$cached_values$mean1 = effects$mean1
			private$cached_values$mean0 = effects$mean0
			private$cached_values$md = effects$md
			private$cached_values$se_md = effects$se_md
		}
	)
)
