#' Marginal Standardization / G-Computation for Ordinal Responses
#'
#' Internal base class for ordinal-outcome g-computation (G-Comp) estimators. A
#' proportional odds model is fit, then potential-outcome expected values under
#' all-treated and all-control assignments are standardized over the empirical
#' covariate distribution. Inference currently supports bootstrap.
#'
#' @keywords internal
#' @noRd
InferenceOrdinalGCompAbstract = R6::R6Class("InferenceOrdinalGCompAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the g-computation (G-Comp) inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an ordinal response.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the g-computation (G-Comp) treatment-effect estimate (mean difference).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$md
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval for the G-Comp mean difference.
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$ensure_md_se()
			z_val = stats::qnorm(1 - alpha / 2)
			ci = private$cached_values$md + c(-1, 1) * z_val * private$cached_values$se_md
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		#' @description
		#' Computes a two-sided Wald p-value for the G-Comp mean difference.
		#' @param delta Description for delta
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$ensure_md_se()
			z_val = (private$cached_values$md - delta) / private$cached_values$se_md
			2 * stats::pnorm(-abs(z_val))
		}
	),

	private = list(
		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$md)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			# X_full expected to have treatment in column 2, no intercept
			X_fit = X_full[, -1, drop = FALSE] 
			j_treat = 1 # Treatment is now first column of X_fit

			# Fit proportional odds model
			fit = tryCatch(
				fast_ordinal_regression_with_var_cpp(X = X_fit, y = as.numeric(private$y)),
				error = function(e) NULL
			)

			if (is.null(fit) || length(fit$b) == 0 || is.null(fit$alpha)){
				private$cached_values$mean1 = NA_real_
				private$cached_values$mean0 = NA_real_
				private$cached_values$md = NA_real_
				private$cached_values$se_md = NA_real_
				return(invisible(NULL))
			}

			coef_hat = as.numeric(fit$b)
			alpha_hat = as.numeric(fit$alpha)

			# Post-fit standardization
			res = gcomp_ordinal_proportional_odds_post_fit_cpp(
				X_fit = X_fit,
				coef_hat = coef_hat,
				alpha_hat = alpha_hat,
				j_treat = j_treat
			)

			private$cached_values$mean1 = res$mean1
			private$cached_values$mean0 = res$mean0
			private$cached_values$md = res$md

			theta = c(alpha_hat, coef_hat)
			vcov_mat = if (!is.null(fit$vcov)) as.matrix(fit$vcov) else NULL
			if (!is.null(vcov_mat) && length(theta) > 0L &&
				 nrow(vcov_mat) == length(theta) && ncol(vcov_mat) == length(theta)){
				grad = private$compute_md_gradient(X_fit, theta, length(alpha_hat), j_treat)
				var_md = as.numeric(crossprod(grad, vcov_mat %*% grad))
				private$cached_values$se_md = if (is.finite(var_md) && var_md > 0) sqrt(var_md) else NA_real_
			} else {
				private$cached_values$se_md = NA_real_
			}
		},

		compute_md_gradient = function(X_fit, theta, n_alpha, j_treat, base_step = 1e-6){
			n_params = length(theta)
			grad = numeric(n_params)
			for (j in seq_len(n_params)){
				step = max(base_step, base_step * (1 + abs(theta[j])))
				theta_plus = theta
				theta_minus = theta
				theta_plus[j] = theta[j] + step
				theta_minus[j] = theta[j] - step
				md_plus = private$compute_md_from_theta(X_fit, theta_plus, n_alpha, j_treat)
				md_minus = private$compute_md_from_theta(X_fit, theta_minus, n_alpha, j_treat)
				grad[j] = (md_plus - md_minus) / (2 * step)
			}
			grad
		},

		compute_md_from_theta = function(X_fit, theta, n_alpha, j_treat){
			alpha_vec = theta[seq_len(n_alpha)]
			coef_vec = theta[(n_alpha + 1):length(theta)]
			as.numeric(
				gcomp_ordinal_proportional_odds_post_fit_cpp(
					X_fit = X_fit,
					coef_hat = coef_vec,
					alpha_hat = alpha_vec,
					j_treat = j_treat
				)$md
			)
		},

		ensure_md_se = function(){
			se = private$cached_values$se_md
			if (!(is.finite(se) && se > 0)){
				stop("Ordinal G-computation: could not compute a finite delta-method standard error.")
			}
			invisible(NULL)
		}
	)
)
