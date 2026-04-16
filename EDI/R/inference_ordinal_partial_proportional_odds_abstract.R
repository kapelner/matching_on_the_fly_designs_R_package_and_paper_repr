#' Internal Base Class for Ordinal Partial Proportional-Odds Inference
#'
#' Shared implementation for ordinal partial proportional-odds estimators. When
#' the requested model has no nonparallel covariates, this class uses the
#' package's fast Rcpp proportional-odds solver before falling back to the
#' general R fitters.
#'
#' @keywords internal
InferenceOrdinalPartialProportionalOddsAbstract = R6::R6Class(
	"InferenceOrdinalPartialProportionalOddsAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description
		#' Initialize the internal PPO base object.
		#' @param des_obj A completed \code{DesignSeqOneByOne} object with an ordinal
		#'   response.
		#' @param nonparallel Covariate names that may vary across thresholds.
		#' @param verbose Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj,
				nonparallel = character(0),
				
				verbose = FALSE,
				harden = TRUE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, verbose, harden)
			assertNoCensoring(private$any_censoring)
			assertCharacter(nonparallel, null.ok = TRUE)
			private$nonparallel = unique(nonparallel)
		},

		#' @description
		#' Retrieve the estimated treatment log-odds shift.
		#'
		#' @return The estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute a Wald-style confidence interval for the treatment effect. If the
		#' model-based standard error is unavailable, falls back to the bootstrap
		#' interval.
		#' @param alpha Significance level for the interval.
		#'
		#' @return A confidence interval for the treatment effect.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(
				alpha,
				lower = .Machine$double.xmin,
				upper = 1 - .Machine$double.xmin
			)
			private$shared()
			if (!private$has_finite_se()){
				warning(
					"Partial proportional-odds regression: falling back to ",
					"bootstrap because standard error is unavailable."
				)
				return(self$compute_bootstrap_confidence_interval(
					alpha = alpha,
					na.rm = TRUE
				))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute a Wald-style two-sided p-value for the treatment effect. If the
		#' model-based standard error is unavailable, falls back to the bootstrap
		#' p-value.
		#' @param delta Null treatment effect to test.
		#'
		#' @return A two-sided p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			if (!private$has_finite_se()){
				warning(
					"Partial proportional-odds regression: falling back to ",
					"bootstrap because standard error is unavailable."
				)
				return(self$compute_bootstrap_two_sided_pval(
					delta = delta,
					na.rm = TRUE
				))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Benchmark the asymptotic p-value path with a timing breakdown.
		#'
		#' This is a diagnostic helper for performance investigation. It separates
		#' the model fit, cache/SE materialization, and final p-value arithmetic.
		#'
		#' @param delta Null treatment effect to test.
		#'
		#' @return A named list with timing and result details.
		benchmark_asymp_two_sided_pval_breakdown = function(delta = 0){
			assertNumeric(delta)

			t0 = proc.time()[["elapsed"]]
			fit = private$fit_partial_proportional_odds()
			fit_time = round(proc.time()[["elapsed"]] - t0, 6)

			if (is.null(fit) || !is.finite(fit$beta)){
				return(list(
					fit_time = fit_time,
					cache_time = NA_real_,
					pval_math_time = NA_real_,
					total_time = fit_time,
					pval = NA_real_,
					beta_hat_T = NA_real_,
					s_beta_hat_T = NA_real_
				))
			}

			t1 = proc.time()[["elapsed"]]
			private$cached_values$beta_hat_T = fit$beta
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = private$n - 1
			cache_time = round(proc.time()[["elapsed"]] - t1, 6)

			t2 = proc.time()[["elapsed"]]
			pval = private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			pval_math_time = round(proc.time()[["elapsed"]] - t2, 6)

			list(
				fit_time = fit_time,
				cache_time = cache_time,
				pval_math_time = pval_math_time,
				total_time = round(fit_time + cache_time + pval_math_time, 6),
				pval = pval,
				beta_hat_T = private$cached_values$beta_hat_T,
				s_beta_hat_T = private$cached_values$s_beta_hat_T
			)
		}
	),
	private = list(
		nonparallel = character(0),

		ppo_covariate_matrix = function(){
			stop(class(self)[1], " must implement ppo_covariate_matrix()")
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			fit = private$fit_partial_proportional_odds()
			if (is.null(fit) || !is.finite(fit$beta)){
				private$cached_values$beta_hat_T = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				private$cached_values$df = private$n - 1
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = fit$beta
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = private$n - 1
		},

		has_finite_se = function(){
			is.finite(private$cached_values$s_beta_hat_T) &&
				private$cached_values$s_beta_hat_T > 0
		},

		fit_partial_proportional_odds = function(){
			X_cov = private$ppo_covariate_matrix()
			X_full = cbind(treatment = private$w, X_cov)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					X_cov_fit = X_fit[, -1, drop = FALSE]
					private$fit_partial_proportional_odds_from_covariates(X_cov_fit)
				},
				fit_ok = function(fit, X_fit, keep){
					private$ppo_fit_is_usable(fit)
				}
			)
			attempt$fit
		},

		ppo_fit_is_usable = function(fit){
			!is.null(fit) && is.finite(fit$beta)
		},

		fit_partial_proportional_odds_from_covariates = function(X_cov){
			covar_names = colnames(X_cov)
			if (is.null(covar_names)) covar_names = character(0)
			nonparallel_covars = intersect(private$nonparallel, covar_names)
			parallel_covars = setdiff(covar_names, nonparallel_covars)

			if (length(nonparallel_covars) == 0){
				fit = private$fit_fast_proportional_odds(X_cov)
				if (!is.null(fit)) return(fit)
			}

			dat = data.frame(
				y = ordered(private$y, levels = sort(unique(private$y))),
				treatment = private$w,
				as.data.frame(X_cov, check.names = FALSE),
				check.names = FALSE
			)
			if (nlevels(dat$y) < 2) return(NULL)

			fit = private$fit_vgam(dat, parallel_covars, nonparallel_covars)
			if (!is.null(fit)) return(fit)

			fit = private$fit_clm(dat, parallel_covars, nonparallel_covars)
			if (!is.null(fit)) return(fit)

			fit = private$fit_polr(dat, parallel_covars, nonparallel_covars)
			if (!is.null(fit)) return(fit)

			NULL
		},

		fit_fast_proportional_odds = function(X_cov){
			X_fit = cbind(treatment = private$w, X_cov)
			if (is.null(dim(X_fit))){
				X_fit = matrix(X_fit, ncol = 1)
				colnames(X_fit) = "treatment"
			}

			res = tryCatch(
				fast_ordinal_regression_with_var_cpp(
					X = X_fit,
					y = as.numeric(private$y)
				),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 1 || !is.finite(res$b[1]) || (isTRUE(private$harden) && !is.null(res$converged) && !res$converged)){
				return(NULL)
			}

			se_beta = if (is.finite(res$ssq_b_2) && res$ssq_b_2 > 0) {
				sqrt(res$ssq_b_2)
			} else {
				NA_real_
			}

			list(beta = as.numeric(res$b[1]), se = se_beta)
		},

		main_formula = function(term_names){
			stats::reformulate(termlabels = term_names, response = "y")
		},

		parallel_formula = function(term_names){
			stats::reformulate(termlabels = term_names)
		},

		extract_common_treatment_fit = function(mod, coef_getter, vcov_getter){
			coefs = tryCatch(coef_getter(mod), error = function(e) NULL)
			if (is.null(coefs) || !"treatment" %in% names(coefs)) return(NULL)

			beta_hat = as.numeric(coefs[["treatment"]])
			var_beta = tryCatch(
				vcov_getter(mod)["treatment", "treatment"],
				error = function(e) NA_real_
			)
			se_beta = if (is.finite(var_beta) && var_beta > 0) sqrt(var_beta) else NA_real_

			list(beta = beta_hat, se = se_beta)
		},

		fit_vgam = function(dat, parallel_covars, nonparallel_covars){
			if (!check_package_installed("VGAM")) return(NULL)

			all_terms = unique(c("treatment", parallel_covars, nonparallel_covars))
			par_terms = unique(c("treatment", parallel_covars))

			mod = tryCatch(
				suppressWarnings(
					VGAM::vglm(
						formula = private$main_formula(all_terms),
						family = VGAM::cumulative(
							link = "logitlink",
							parallel = private$parallel_formula(par_terms)
						),
						data = dat,
						trace = FALSE,
						model = FALSE
					)
				),
				error = function(e) NULL
			)
			if (is.null(mod)) return(NULL)

			private$extract_common_treatment_fit(
				mod,
				coef_getter = VGAM::Coef,
				vcov_getter = VGAM::vcov
			)
		},

		fit_clm = function(dat, parallel_covars, nonparallel_covars){
			if (!check_package_installed("ordinal")) return(NULL)

			main_terms = unique(c("treatment", parallel_covars))
			nominal_form = if (length(nonparallel_covars) == 0) {
				NULL
			} else {
				stats::reformulate(termlabels = nonparallel_covars)
			}

			mod = tryCatch(
				suppressWarnings(
					ordinal::clm(
						formula = private$main_formula(main_terms),
						nominal = nominal_form,
						data = dat,
						link = "logit",
						Hess = TRUE
					)
				),
				error = function(e) NULL
			)
			if (is.null(mod)) return(NULL)

			private$extract_common_treatment_fit(
				mod,
				coef_getter = stats::coef,
				vcov_getter = stats::vcov
			)
		},

		fit_polr = function(dat, parallel_covars, nonparallel_covars){
			if (length(nonparallel_covars) > 0) return(NULL)
			main_terms = unique(c("treatment", parallel_covars))
			mod = tryCatch(
				suppressWarnings(
					MASS::polr(
						formula = private$main_formula(main_terms),
						data = dat,
						method = "logistic",
						Hess = TRUE
					)
				),
				error = function(e) NULL
			)
			if (is.null(mod)) return(NULL)

			private$extract_common_treatment_fit(
				mod,
				coef_getter = stats::coef,
				vcov_getter = stats::vcov
			)
		}
	)
)
