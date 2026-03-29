#' Stratified Cox PH Inference for Survival Responses
#'
#' Internal base class for all-subject non-KK stratified Cox proportional hazards
#' regression. Stratification variables are chosen automatically from the recorded
#' low-cardinality covariates. If no suitable stratification covariates are found,
#' the fit falls back to the corresponding standard Cox PH model.
#'
#' @keywords internal
#' @noRd
InferenceSurvivalStratCoxPHAbstract = R6::R6Class("InferenceSurvivalStratCoxPHAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(


		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute asymp two sided pval for treatment effect
		#' @param delta Description for delta
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Compute confidence interval rand
		#' @param alpha Description for alpha
		#' @param r Number of vectors to draw.
		#' @param pval_epsilon Description for pval_epsilon
		#' @param show_progress Description for show_progress
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE){
			stop("Randomization confidence intervals are not supported for stratified Cox PH models because the estimator units (Log-Hazard Ratio) are inconsistent with the randomization test's required transformed scale (Log-Time Ratio / AFT effect).")
		}
	),

	private = list(
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			mod = private$generate_mod()
			private$cached_values$beta_hat_T = as.numeric(mod$b[2])
			if (estimate_only) return(invisible(NULL))
			se = if (is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0) sqrt(mod$ssq_b_2) else NA_real_
			private$cached_values$s_beta_hat_T = se
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
		},

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()."),

		compute_strata_info = function(X_full){
			n = length(private$y)
			if (is.null(X_full) || ncol(X_full) == 0){
				return(list(strata_id = rep.int(1L, n), selected_cols = integer(0), num_strata = 1L))
			}

			info = tryCatch(
				compute_survival_strata_ids_cpp(as.matrix(X_full)),
				error = function(e) NULL
			)
			if (is.null(info)){
				return(list(strata_id = rep.int(1L, n), selected_cols = integer(0), num_strata = 1L))
			}

			list(
				strata_id = as.integer(info$strata_id),
				selected_cols = as.integer(info$selected_cols),
				num_strata = as.integer(info$num_strata)
			)
		},

		reduce_covariates_preserving_treatment = function(X_covars){
			if (is.null(X_covars) || ncol(X_covars) == 0){
				return(matrix(numeric(0), nrow = length(private$y), ncol = 0))
			}

			X_covars = as.matrix(X_covars)
			if (ncol(X_covars) == 0){
				return(matrix(numeric(0), nrow = nrow(X_covars), ncol = 0))
			}

			full_design = cbind(w = private$w, X_covars)
			reduced = drop_linearly_dependent_cols(full_design)
			X_keep = reduced$M
			if (ncol(X_keep) == 0){
				return(matrix(numeric(0), nrow = nrow(full_design), ncol = 0))
			}
			if (!("w" %in% colnames(X_keep))){
				return(matrix(numeric(0), nrow = nrow(full_design), ncol = 0))
			}

			X_keep[, colnames(X_keep) != "w", drop = FALSE]
		},

		get_informative_rows = function(strata_id){
			if (length(strata_id) != length(private$y)) return(integer(0))
			good = rep(FALSE, length(strata_id))
			for (s in unique(strata_id)){
				i_s = which(strata_id == s)
				if (length(i_s) < 2) next
				if (length(unique(private$w[i_s])) < 2) next
				if (!any(private$dead[i_s] == 1, na.rm = TRUE)) next
				good[i_s] = TRUE
			}
			which(good)
		},

		fit_cox_with_formula = function(dat, formula_str){
			tryCatch(
				suppressWarnings(survival::coxph(stats::as.formula(formula_str), data = dat)),
				error = function(e) NULL
			)
		},

		format_mod_output = function(mod){
			if (is.null(mod)){
				return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_))
			}
			coef_w = tryCatch(as.numeric(stats::coef(mod)["w"]), error = function(e) NA_real_)
			ssq_w = tryCatch(as.numeric(stats::vcov(mod)["w", "w"]), error = function(e) NA_real_)
			list(
				b = c(0, coef_w),
				ssq_b_2 = if (is.finite(ssq_w) && ssq_w > 0) ssq_w else NA_real_
			)
		},

		generate_mod = function(){
			surv_obj = survival::Surv(private$y, private$dead)
			X_full = private$get_X()
			strata_info = private$compute_strata_info(X_full)

			X_linear = matrix(numeric(0), nrow = length(private$y), ncol = 0)
			if (private$include_covariates() && !is.null(X_full) && ncol(X_full) > 0){
				keep_cols = setdiff(seq_len(ncol(X_full)), strata_info$selected_cols)
				if (length(keep_cols) > 0){
					X_linear = private$reduce_covariates_preserving_treatment(X_full[, keep_cols, drop = FALSE])
				}
			}

			dat_full = data.frame(y = private$y, dead = private$dead, w = private$w)
			if (ncol(X_linear) > 0){
				colnames(X_linear) = paste0("x", seq_len(ncol(X_linear)))
				dat_full = cbind(dat_full, as.data.frame(X_linear))
			}

			base_terms = c("w", colnames(dat_full)[!(colnames(dat_full) %in% c("y", "dead", "w"))])
			base_formula = paste("survival::Surv(y, dead) ~", paste(base_terms, collapse = " + "))

			# Prefer a genuinely stratified fit if the auto-strata carry information.
			informative_rows = integer(0)
			if (!is.null(strata_info$strata_id) && isTRUE(strata_info$num_strata > 1L)){
				informative_rows = private$get_informative_rows(strata_info$strata_id)
			}

			if (length(informative_rows) >= 4){
				dat_strat = dat_full[informative_rows, , drop = FALSE]
				dat_strat$strata_id = factor(strata_info$strata_id[informative_rows])
				strat_formula = paste(base_formula, "+ strata(strata_id)")

				mod = private$fit_cox_with_formula(dat_strat, strat_formula)
				if (!is.null(mod)){
					return(private$format_mod_output(mod))
				}

				if (private$include_covariates() && ncol(X_linear) > 0){
					mod = private$fit_cox_with_formula(
						dat_strat[, c("y", "dead", "w", "strata_id"), drop = FALSE],
						"survival::Surv(y, dead) ~ w + strata(strata_id)"
					)
					if (!is.null(mod)){
						return(private$format_mod_output(mod))
					}
				}
			}

			# Fallback to standard Cox using all subjects.
			mod = private$fit_cox_with_formula(dat_full, base_formula)
			if (is.null(mod) && private$include_covariates() && ncol(X_linear) > 0){
				mod = private$fit_cox_with_formula(dat_full[, c("y", "dead", "w"), drop = FALSE], "survival::Surv(y, dead) ~ w")
			}
			private$format_mod_output(mod)
		}
	)
)
