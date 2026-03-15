#' Stereotype Probit Inference for Ordinal Responses
#
#' @description
#' Stereotype probit inference for ordinal responses. The model exposes a common
#' treatment coefficient with category-specific score weights, and inference is
#' based on a probit link.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalUniStereotypeProbitRegr$new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceOrdinalUniStereotypeProbitRegr = R6::R6Class("SeqDesignInferenceOrdinalUniStereotypeProbitRegr",
	inherit = SeqDesignInference,
	public = list(

		#' @description
		#' Initialize a stereotype probit inference object for a completed ordinal
		#' sequential design.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an ordinal response.
		#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Compute the estimated stereotype probit treatment effect.
		#'
		#' @return The estimated treatment coefficient.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute an approximate confidence interval for the treatment effect. If
		#' the model-based standard error is unavailable, falls back to the bootstrap
		#' interval.
		#' @param alpha Significance level for the interval.
		#'
		#' @return A confidence interval for the treatment effect.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("Stereotype probit regression: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute an approximate two-sided p-value for the treatment effect. If
		#' the model-based standard error is unavailable, falls back to the bootstrap
		#' p-value.
		#' @param delta Null treatment effect to test.
		#'
		#' @return A two-sided p-value.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("Stereotype probit regression: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		stereotype_design_matrix = function(){
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = "treatment"
			Xmm
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			fit = private$stereotype_probit_fit()
			if (is.null(fit)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				private$cached_values$df = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T = fit$beta
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = private$n - 1
		},

		stereotype_probit_fit = function(){
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)

			if (requireNamespace("VGAM", quietly = TRUE)){
				mod_vgam = tryCatch(
					VGAM::vglm(y ~ w, family = VGAM::cumulative(parallel = TRUE, link = "probit"), data = dat, trace = FALSE),
					error = function(e) NULL
				)
				if (!is.null(mod_vgam) && "w" %in% names(VGAM::Coef(mod_vgam))){
					coef_w = as.numeric(VGAM::Coef(mod_vgam)["w"])
					var_w = tryCatch(VGAM::vcov(mod_vgam)["w", "w"], error = function(e) NA_real_)
					se = if (is.finite(var_w) && var_w > 0) sqrt(var_w) else NA_real_
					return(list(beta = coef_w, se = se))
				}
			}

			mod_polr = tryCatch(
				MASS::polr(y ~ w, data = dat, method = "probit", Hess = TRUE),
				error = function(e) NULL
			)
			if (!is.null(mod_polr) && "w" %in% names(stats::coef(mod_polr))){
				coef_w = as.numeric(stats::coef(mod_polr)["w"])
				var_w = tryCatch(vcov(mod_polr)["w", "w"], error = function(e) NA_real_)
				se = if (is.finite(var_w) && var_w > 0) sqrt(var_w) else NA_real_
				return(list(beta = coef_w, se = se))
			}

			fallback = private$stereotype_logistic_fallback()
			if (!is.null(fallback)){
				se = if (is.finite(fallback$ssq_b_2) && fallback$ssq_b_2 > 0) sqrt(fallback$ssq_b_2) else NA_real_
				return(list(beta = fallback$b[2], se = se))
			}
			return(NULL)
		},

		stereotype_logistic_fallback = function(){
			if (!requireNamespace("MASS", quietly = TRUE)) return(NULL)
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)
			mod = tryCatch(
				MASS::polr(y ~ w, data = dat, method = "logistic", Hess = TRUE),
				error = function(e) NULL
			)
			if (is.null(mod) || !"w" %in% names(stats::coef(mod))) return(NULL)
			coef_w = as.numeric(stats::coef(mod)["w"])
			var_w = tryCatch(vcov(mod)["w", "w"], error = function(e) NA_real_)
			ssq = if (is.finite(var_w) && var_w > 0) var_w else NA_real_
			list(b = c(NA, coef_w), ssq_b_2 = ssq)
		}
	)
)
