#' Zero/One-Inflated Beta Inference for Proportion Responses
#'
#' @description
#' Internal base class for non-KK zero/one-inflated beta regression models. The
#' response is modeled as a three-component mixture with point masses at 0 and 1
#' plus a beta-distributed interior component on \eqn{(0, 1)}. The reported
#' treatment effect is the treatment coefficient from the beta mean submodel, on
#' the logit scale.
#'
#' @details
#' The inflation probabilities for exact 0 and exact 1 are intercept-only in this
#' first implementation. The beta mean submodel uses treatment alone in the
#' univariate class and treatment plus covariates in the multivariate class.
#'
#' @keywords internal
#' @noRd
InferencePropZeroOneInflatedBetaAbstract = R6::R6Class("InferencePropZeroOneInflatedBetaAbstract",
	inherit = InferenceAsymp,
	public = list(

		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "proportion")
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
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("Zero/one-inflated beta estimator: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("Zero/one-inflated beta estimator: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		build_design_matrix_candidates = function(){
			X_cov_orig = as.matrix(private$predictors_df())
			if (ncol(X_cov_orig) == 0L){
				M = matrix(private$w, ncol = 1)
				colnames(M) = "treatment"
				return(list(M))
			}

			thresholds = c(Inf, 0.95, 0.90, 0.80, 0.70)
			candidates = list()
			keys = character()
			for (thresh in thresholds){
				X_cov = if (is.finite(thresh)) drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M else X_cov_orig
				M = cbind(treatment = private$w, X_cov)
				qr_M = qr(M)
				if (qr_M$rank < ncol(M)){
					keep = qr_M$pivot[seq_len(qr_M$rank)]
					if (!(1L %in% keep)) keep = c(1L, keep)
					keep = sort(unique(keep))
					M = M[, keep, drop = FALSE]
				}
				colnames(M)[1] = "treatment"
				key = paste(colnames(M), collapse = "|")
				if (!(key %in% keys)){
					candidates[[length(candidates) + 1L]] = M
					keys = c(keys, key)
				}
			}
			candidates
		},

		predictors_df = function(){
			matrix(numeric(0), nrow = private$n, ncol = 0)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				stop("Zero/one-inflated beta estimator: could not compute a finite standard error.")
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

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			fit = NULL
			for (Xmm in private$build_design_matrix_candidates()){
				fit = .fit_zero_one_inflated_beta(private$y, Xmm)
				if (!is.null(fit)) break
			}
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			coef_full = fit$coefficients
			vcov_full = fit$vcov
			if (!("treatment" %in% names(coef_full)) || !("treatment" %in% rownames(vcov_full))){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			se_full = sqrt(pmax(diag(vcov_full), 0))
			z_vals = coef_full / se_full
			summary_table = cbind(
				Value = coef_full,
				`Std. Error` = se_full,
				`z value` = z_vals,
				`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
			)

			private$cached_values$beta_hat_T = as.numeric(coef_full["treatment"])
			private$cached_values$s_beta_hat_T = as.numeric(se_full["treatment"])
			private$cached_values$is_z = TRUE
			private$cached_values$df = private$n - length(coef_full)
			private$cached_values$full_coefficients = coef_full
			private$cached_values$full_vcov = vcov_full
			private$cached_values$summary_table = summary_table
			private$cached_values$zoib_best_design_matrix = Xmm
			private$cached_values$zoib_best_fit = list(
				coefficients = coef_full,
				vcov = vcov_full
			)
		}
	)
)
