#' Marginal Standardization / G-Computation for Proportion Responses
#'
#' @description
#' Internal base class for proportion-outcome g-computation estimators. A
#' fractional logit working model is fit, then potential-outcome mean
#' proportions under all-treated and all-control assignments are standardized
#' over the empirical covariate distribution. Inference uses a sandwich-robust
#' covariance for the regression coefficients and the delta method for the
#' marginal mean-difference estimand.
#'
#' @keywords internal
#' @noRd
DesignInferencePropGCompAbstract = R6::R6Class("DesignInferencePropGCompAbstract",
	inherit = DesignInference,
	public = list(

		# @description
		# Initialize the g-computation inference object.
		# @param seq_des_obj A completed \code{SeqDesign} object with a proportion response.
		# @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		# @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "proportion")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		# @description
		# Computes the g-computation treatment-effect estimate.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$md
		},

		# @description
		# Computes a 1 - \code{alpha} confidence interval.
		# @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_effect_confidence_interval(alpha)
		},

		# @description
		# Computes a two-sided p-value for the treatment effect.
		# @param delta The null mean difference. Defaults to 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta, len = 1)
			private$shared()
			private$compute_effect_pvalue(delta)
		},

		# @description
		# Computes a bootstrap two-sided p-value for the treatment effect.
		# @param delta The null mean difference. Defaults to 0.
		# @param B Number of bootstrap samples.
		# @param na.rm Whether to remove non-finite bootstrap replicates.
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE){
			assertNumeric(delta, len = 1)
			super$compute_bootstrap_two_sided_pval(delta = delta, B = B, na.rm = na.rm)
		},

		#' @description
		#' Abbreviated bootstrap sampler that reuses the reduced column layout.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, max_resample_attempts = 50){
			assertCount(B, positive = TRUE)
			assertCount(max_resample_attempts, positive = TRUE)
			private$shared()
			cols = private$gcomp_design_colnames
			if (is.null(cols)){
				return(super$approximate_bootstrap_distribution_beta_hat_T(B, max_resample_attempts))
			}

			X_full = private$build_named_design_matrix()
			X_fit = X_full[, cols, drop = FALSE]
			n = nrow(X_fit)
			y = private$y
			w = private$w

			draw_sample = function(...){
				set_package_threads(1L)
				attempt = 1
				repeat {
					i_b = sample_int_replace_cpp(n, n)
					w_b = w[i_b]
					if (any(w_b == 1, na.rm = TRUE) && any(w_b == 0, na.rm = TRUE)) break
					attempt = attempt + 1
					if (attempt > max_resample_attempts) return(NA_real_)
				}
				X_b = X_fit[i_b, , drop = FALSE]
				y_b = y[i_b]
				private$bootstrap_effect_from_sample(X_b, y_b)
			}

			if (private$num_cores == 1){
				vapply(seq_len(B), function(b) draw_sample(), numeric(1))
			} else {
				unlist(parallel::mclapply(seq_len(B), function(b) draw_sample(), mc.cores = min(2L, private$num_cores)))
			}
		}
	),

	private = list(
		gcomp_design_colnames = NULL,
		gcomp_design_j_treat = NULL,
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

		effects_are_usable = function(effects){
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

		fit_fractional_logit_with_sandwich = function(X_full){
			X_curr = X_full

			repeat {
				reduced = private$reduce_design_matrix_preserving_treatment(X_curr)
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
					return(NULL)
				}

				mod = tryCatch(
					suppressWarnings(stats::glm.fit(x = X_fit, y = as.numeric(private$y), family = stats::binomial(link = "logit"))),
					error = function(e) NULL
				)
				if (is.null(mod)){
					return(NULL)
				}

				coef_hat = as.numeric(mod$coefficients)
				converged = isTRUE(mod$converged) && length(coef_hat) == ncol(X_fit) && all(is.finite(coef_hat))
				if (!converged){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				mu_hat = pmin(pmax(as.numeric(mod$fitted.values), .Machine$double.eps), 1 - .Machine$double.eps)
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
					post_fit = post_fit
				))
			}
		},

		compute_standardized_effects_r = function(fit){
			X_fit = fit$X
			coef_hat = fit$coefficients
			vcov_robust = fit$vcov
			j_treat = fit$j_treat

			X1 = X_fit
			X0 = X_fit
			X1[, j_treat] = 1
			X0[, j_treat] = 0

			eta1 = as.numeric(X1 %*% coef_hat)
			eta0 = as.numeric(X0 %*% coef_hat)
			mean1_i = stats::plogis(eta1)
			mean0_i = stats::plogis(eta0)
			mean1 = mean(mean1_i)
			mean0 = mean(mean0_i)

			grad1 = as.numeric(crossprod(X1, mean1_i * (1 - mean1_i))) / nrow(X1)
			grad0 = as.numeric(crossprod(X0, mean0_i * (1 - mean0_i))) / nrow(X0)

			md = mean1 - mean0
			grad_md = grad1 - grad0
			var_md = as.numeric(t(grad_md) %*% vcov_robust %*% grad_md)

			std_err = sqrt(pmax(diag(vcov_robust), 0))
			z_vals = coef_hat / std_err
			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)

			list(
				mean1 = mean1,
				mean0 = mean0,
				md = md,
				se_md = if (is.finite(var_md) && var_md >= 0) sqrt(var_md) else NA_real_,
				full_coefficients = coef_hat,
				full_vcov = vcov_robust,
				summary_table = summary_table
			)
		},

		compute_standardized_effects = function(fit){
			coef_hat = fit$coefficients
			fast = fit$post_fit
			if (is.null(fast)){
				return(private$compute_standardized_effects_r(fit))
			}

			vcov_robust = fast$vcov
			colnames(vcov_robust) = rownames(vcov_robust) = names(coef_hat)
			std_err = fast$std_err
			names(std_err) = names(coef_hat)
			z_vals = fast$z_vals
			names(z_vals) = names(coef_hat)
			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = std_err,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)

			list(
				mean1 = fast$mean1,
				mean0 = fast$mean0,
				md = fast$md,
				se_md = fast$se_md,
				full_coefficients = coef_hat,
				full_vcov = vcov_robust,
				summary_table = summary_table
			)
		},

		bootstrap_effect_from_sample = function(X_b, y_b){
			if (nrow(X_b) == 0) return(NA_real_)
			mod = tryCatch(
				suppressWarnings(stats::glm.fit(
					x = X_b,
					y = as.numeric(y_b),
					family = stats::binomial(link = "logit")
				)),
				error = function(e) NULL
			)
			if (is.null(mod)){
				return(NA_real_)
			}
			coef_hat = as.numeric(mod$coefficients)
			if (!isTRUE(mod$converged) || length(coef_hat) != ncol(X_b) || any(!is.finite(coef_hat))){
				return(NA_real_)
			}

			mu_hat = pmin(pmax(as.numeric(mod$fitted.values), .Machine$double.eps), 1 - .Machine$double.eps)
			post_fit = tryCatch(
				gcomp_fractional_logit_post_fit_cpp(
					X_fit = X_b,
					y = as.numeric(y_b),
					coef_hat = coef_hat,
					mu_hat = mu_hat,
					j_treat = private$gcomp_design_j_treat
				),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				return(NA_real_)
			}

			md_val = post_fit$md
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
				stop("G-computation mean difference: could not compute a finite delta-method standard error.")
			}
			ci = est + c(-1, 1) * z * se
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		compute_effect_pvalue = function(delta){
			est = private$cached_values$md
			se = private$cached_values$se_md
			if (!is.finite(est) || !is.finite(se) || se <= 0){
				stop("G-computation mean difference: could not compute a finite delta-method standard error.")
			}
			z_stat = (est - delta) / se
			2 * stats::pnorm(-abs(z_stat))
		},

		shared = function(){
			if (!is.null(private$cached_values$summary_table)) return(invisible(NULL))

			X_full = private$build_named_design_matrix()

			fit = private$fit_fractional_logit_with_sandwich(X_full)
			effects = if (!is.null(fit)) private$compute_standardized_effects(fit) else NULL
			if ((is.null(fit) || is.null(effects) || !private$effects_are_usable(effects)) && ncol(X_full) > 2L){
				fit = private$fit_fractional_logit_with_sandwich(X_full[, 1:2, drop = FALSE])
				effects = if (!is.null(fit)) private$compute_standardized_effects(fit) else NULL
			}
			if (is.null(fit) || is.null(effects) || !private$effects_are_usable(effects)){
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
