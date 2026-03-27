#' Constrained Binomial-Link Inference for Binary Responses
#'
#' @description
#' Internal base class for binomial GLMs with constrained mean links on
#' incidence outcomes. Shared logic handles design-matrix reduction, covariate
#' dropping on failed fits, and cached Wald inference.
#'
#' @keywords internal
#' @noRd
InferenceIncidConstrainedBinomialAbstract = R6::R6Class("InferenceIncidConstrainedBinomialAbstract",
	inherit = InferenceAsymp,
	public = list(

		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),
		fit_constrained_binomial = function(X_fit, j_treat) stop(class(self)[1], " must implement fit_constrained_binomial()."),
		method_label = function() stop(class(self)[1], " must implement method_label()."),

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
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

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				stop(private$method_label(), ": could not compute a finite standard error.")
			}
		},

		set_failed_fit_cache = function(){
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},

		fit_model_with_dropping = function(X_full){
			X_curr = X_full

			repeat {
				reduced = private$reduce_design_matrix_preserving_treatment(X_curr)
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
					return(NULL)
				}

				mod = tryCatch(
					private$fit_constrained_binomial(X_fit, j_treat),
					error = function(e) NULL
				)
				if (is.null(mod) || !isTRUE(mod$converged)){
					if (ncol(X_curr) <= 2L) return(NULL)
					coef_hat = if (!is.null(mod) && !is.null(mod$b)) as.numeric(mod$b) else rep(NA_real_, ncol(X_fit))
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				coef_hat = as.numeric(mod$b)
				ssq_treat = as.numeric(mod$ssq_b_j)
				if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat)) || !is.finite(ssq_treat) || ssq_treat <= 0){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				coef_names = colnames(X_fit)
				names(coef_hat) = coef_names
				vcov_model = mod$vcov
				colnames(vcov_model) = rownames(vcov_model) = coef_names
				std_err = as.numeric(mod$std_err)
				names(std_err) = coef_names
				z_vals = as.numeric(mod$z_vals)
				names(z_vals) = coef_names

				return(list(
					beta_hat = coef_hat[j_treat],
					se = sqrt(ssq_treat),
					coefficients = coef_hat,
					vcov = vcov_model,
					summary_table = cbind(
						Value = coef_hat,
						`Std. Error` = std_err,
						`z value` = z_vals,
						`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
					)
				))
			}
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)

			fit = private$fit_model_with_dropping(X_full)
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = fit$beta_hat
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_full) - ncol(fit$vcov)
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$full_vcov = fit$vcov
			private$cached_values$summary_table = fit$summary_table
		}
	)
)
