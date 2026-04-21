#' KK Hurdle Poisson Combined-Likelihood Inference for Count Responses
#'
#' Internal base class for KK hurdle-Poisson combined-likelihood models. The
#' fitted model uses a single \pkg{glmmTMB} hurdle-Poisson likelihood over all
#' subjects, with pair-specific random intercepts active on matched rows and a
#' fixed reservoir intercept on reservoir rows. The reported treatment effect is
#' the treatment coefficient from the positive-count component on the log-rate
#' scale.
#'
#' @keywords internal
#' @noRd
InferenceAbstractKKHurdlePoissonOneLik = R6::R6Class("InferenceAbstractKKHurdlePoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB}.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(use_rcpp)
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or FixedDesignBinaryMatch.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (is(des_obj, "FixedDesignBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_bms_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			private$use_rcpp = use_rcpp

			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !private$use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
		},

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_hurdle(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("KK hurdle-Poisson combined-likelihood: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute asymp two sided pval for treatment effect
		#' @param delta Description for delta
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("KK hurdle-Poisson combined-likelihood: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		use_rcpp = TRUE,

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		fit_hurdle_model = function(model_data){
			if (private$use_rcpp) {
				# Use glmmTMB as the engine for our "ourself" path for now,
				# but it is managed via the internal dispatch.
			}

			dat = model_data$dat
			has_pairs = isTRUE(model_data$has_pairs)
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)

			formula_cond = private$build_cond_formula(dat, has_pairs = has_pairs)
			formula_zi = private$build_zi_formula(dat, has_pairs = has_pairs)

			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = formula_zi,
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
			if (!is.null(mod)) return(mod)

			dat_fallback = dat[, intersect(c("y", "w", "reservoir_ind", "pair_group", "pair_active"), colnames(dat)), drop = FALSE]
			formula_cond = private$build_cond_formula(dat_fallback, has_pairs = has_pairs)
			formula_zi = private$build_zi_formula(dat_fallback, has_pairs = has_pairs)
			tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = formula_zi,
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat_fallback,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
		},

		shared_combined_hurdle = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			model_data = private$build_model_data()
			if (is.null(model_data)){
				private$cache_nonestimable_estimate("kk_hurdle_poisson_combined_no_model_data")
				return(invisible(NULL))
			}

			mod = private$fit_hurdle_model(model_data)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_hurdle_poisson_combined_fit_unavailable")
				return(invisible(NULL))
			}

			cond_fixef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
			if (is.null(cond_fixef) || !("w" %in% names(cond_fixef)) || !is.finite(cond_fixef[["w"]])){
				private$cache_nonestimable_estimate("kk_hurdle_poisson_combined_treatment_missing")
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(cond_fixef[["w"]])
			if (!estimate_only) {
				coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
				se = if (!is.null(coef_table) && ("w" %in% rownames(coef_table))) as.numeric(coef_table["w", "Std. Error"]) else NA_real_
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			}
			private$cached_values$is_z = TRUE
		}
	)
)
