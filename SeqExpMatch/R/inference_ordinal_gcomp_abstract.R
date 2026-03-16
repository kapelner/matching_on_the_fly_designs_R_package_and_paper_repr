#' Marginal Standardization / G-Computation for Ordinal Responses
#'
#' @description
#' Internal base class for ordinal-outcome g-computation estimators. A
#' proportional odds model is fit, then potential-outcome expected values under
#' all-treated and all-control assignments are standardized over the empirical
#' covariate distribution. Inference currently supports bootstrap.
#'
#' @keywords internal
#' @noRd
SeqDesignInferenceOrdinalGCompAbstract = R6::R6Class("SeqDesignInferenceOrdinalGCompAbstract",
	inherit = SeqDesignInference,
	public = list(

		# @description
		# Initialize the g-computation inference object.
		# @param seq_des_obj A completed \code{SeqDesign} object with an ordinal response.
		# @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		# @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		# @description
		# Computes the g-computation treatment-effect estimate (mean difference).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$md
		},

		# @description
		# Computes a 1 - \code{alpha} confidence interval.
		# @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_effect_confidence_interval(alpha)
		},

		# @description
		# Computes a two-sided p-value for the treatment effect.
		# @param delta The null mean difference. Defaults to 0.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta, len = 1)
			private$shared()
			private$compute_effect_pvalue(delta)
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

		set_failed_fit_cache = function(){
			private$cached_values$mean1 = NA_real_
			private$cached_values$mean0 = NA_real_
			private$cached_values$md = NA_real_
			private$cached_values$se_md = NA_real_
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},

		effects_are_usable = function(effects){
			is.finite(effects$md) && is.finite(effects$se_md) && effects$se_md > 0
		},

		compute_standardized_effects = function(X_fit, fit){
			coef_hat = as.numeric(fit$b)
			alpha_hat = as.numeric(fit$alpha)
			fast = tryCatch(
				ordinal_gcomp_post_fit_cpp(
					X_fit = X_fit,
					y = as.numeric(private$y),
					coef_hat = coef_hat,
					alpha_hat = alpha_hat,
					j_treat = 1L
				),
				error = function(e) NULL
			)
			if (is.null(fast)){
				return(NULL)
			}

			vcov_beta = fast$vcov
			colnames(vcov_beta) = rownames(vcov_beta) = colnames(X_fit)
			se_beta = fast$std_err
			names(se_beta) = colnames(X_fit)
			z_vals = fast$z_vals
			names(z_vals) = colnames(X_fit)
			summary_table = cbind(
				Value = coef_hat,
				`Std. Error` = se_beta,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)

			list(
				mean1 = fast$mean1,
				mean0 = fast$mean0,
				md = fast$md,
				se_md = fast$se_md,
				full_coefficients = coef_hat,
				full_vcov = vcov_beta,
				summary_table = summary_table
			)
		},

		compute_effect_confidence_interval = function(alpha){
			z = stats::qnorm(1 - alpha / 2)
			est = private$cached_values$md
			se = private$cached_values$se_md
			if (!is.finite(est) || !is.finite(se) || se <= 0){
				stop("Ordinal g-computation mean difference: could not compute a finite delta-method standard error.")
			}
			ci = est + c(-1, 1) * z * se
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		compute_effect_pvalue = function(delta){
			est = private$cached_values$md
			se = private$cached_values$se_md
			if (!is.finite(est) || !is.finite(se) || se <= 0){
				stop("Ordinal g-computation mean difference: could not compute a finite delta-method standard error.")
			}
			z_stat = (est - delta) / se
			2 * stats::pnorm(-abs(z_stat))
		},

		shared = function(){
			if (!is.null(private$cached_values$summary_table)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			# X_full expected to have treatment in column 2, no intercept (fast_ordinal handles alphas)
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)
			X_fit = X_full[, -1, drop = FALSE]

			fit = tryCatch(
				fast_ordinal_regression_cpp(X = X_fit, y = as.numeric(private$y)),
				error = function(e) NULL
			)
			effects = if (!is.null(fit) && length(fit) > 0) private$compute_standardized_effects(X_fit, fit) else NULL

			if ((is.null(fit) || is.null(effects) || !private$effects_are_usable(effects)) && ncol(X_full) > 2L){
				X_fit = X_full[, 2, drop = FALSE]
				fit = tryCatch(
					fast_ordinal_regression_cpp(X = X_fit, y = as.numeric(private$y)),
					error = function(e) NULL
				)
				effects = if (!is.null(fit) && length(fit) > 0) private$compute_standardized_effects(X_fit, fit) else NULL
			}

			if (is.null(fit) || is.null(effects) || !private$effects_are_usable(effects)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$mean1 = effects$mean1
			private$cached_values$mean0 = effects$mean0
			private$cached_values$md = effects$md
			private$cached_values$se_md = effects$se_md
			private$cached_values$full_coefficients = effects$full_coefficients
			private$cached_values$full_vcov = effects$full_vcov
			private$cached_values$summary_table = effects$summary_table
		}
	)
)
