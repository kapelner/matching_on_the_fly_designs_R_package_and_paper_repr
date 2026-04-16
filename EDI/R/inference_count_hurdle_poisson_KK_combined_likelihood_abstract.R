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
InferenceAbstractKKHurdlePoissonCombinedLikelihood = R6::R6Class("InferenceAbstractKKHurdlePoissonCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "count")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
			if (!check_package_installed("glmmTMB")){
				stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
			}
		},

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_combined_hurdle(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
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
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("KK hurdle-Poisson combined-likelihood: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()."),

		predictors_df = function(){
			data.frame(w = private$w)
		},

		build_fixed_term_names = function(dat){
			terms = setdiff(colnames(dat), c("y", "pair_group", "pair_active"))
			terms
		},

		build_cond_formula = function(dat, has_pairs){
			terms = private$build_fixed_term_names(dat)
			rhs = paste(terms, collapse = " + ")
			if (!nzchar(rhs)) rhs = "1"
			if (has_pairs){
				rhs = paste(rhs, "(0 + pair_active | pair_group)", sep = " + ")
			}
			stats::as.formula(paste("y ~", rhs))
		},

		build_zi_formula = function(dat, has_pairs){
			terms = private$build_fixed_term_names(dat)
			rhs = paste(terms, collapse = " + ")
			if (!nzchar(rhs)) rhs = "1"
			if (has_pairs){
				rhs = paste(rhs, "(0 + pair_active | pair_group)", sep = " + ")
			}
			stats::as.formula(paste("~", rhs))
		},

		build_model_data = function(){
			if (is.null(private$m)){
				private$m = rep(0L, private$n)
			}
			m_vec = as.integer(private$m)
			m_vec[is.na(m_vec)] = 0L

			pred_df = private$predictors_df()
			if (!("w" %in% colnames(pred_df))){
				stop(class(self)[1], " predictors_df() must include a treatment column named 'w'.")
			}

			X_full = cbind(1, as.matrix(pred_df))
			colnames(X_full)[1:2] = c("(Intercept)", "w")
			X_reduced = private$reduce_design_matrix_preserving_treatment_matrix(X_full)
			if (is.null(X_reduced)){
				return(NULL)
			}

			pred_df = as.data.frame(X_reduced[, -1, drop = FALSE])
			colnames(pred_df)[1] = "w"

			grouping = compute_kk_grouping_cpp(m_vec)
			dat = data.frame(
				y = private$y,
				pred_df,
				reservoir_ind = grouping$reservoir_ind,
				pair_active = grouping$pair_active,
				pair_group = factor(grouping$cluster_id)
			)

			has_pairs = any(grouping$pair_active > 0.5)
			has_reservoir = any(grouping$reservoir_ind > 0.5)
			if (!has_reservoir){
				dat$reservoir_ind = NULL
			}

			list(
				dat = dat,
				has_pairs = has_pairs,
				has_reservoir = has_reservoir
			)
		},

		fit_hurdle_model = function(model_data){
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
