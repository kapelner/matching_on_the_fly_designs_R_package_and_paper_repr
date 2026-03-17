#' Zero-Augmented Poisson Inference for Count Responses
#'
#' @description
#' Internal base class for non-KK zero-inflated and hurdle Poisson regression
#' models fit using the \pkg{glmmTMB} fitter. The reported treatment effect is the
#' treatment coefficient from the conditional count component, on the log-rate
#' scale.
#'
#' @keywords internal
#' @noRd
SeqDesignInferenceCountZeroAugmentedPoissonAbstract = R6::R6Class("SeqDesignInferenceCountZeroAugmentedPoissonAbstract",
	inherit = SeqDesignInference,
	public = list(

		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "count")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			if (!requireNamespace("glmmTMB", quietly = TRUE)){
				stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
			}
		},

		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning(private$za_description(), ": falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		za_family = function() stop(class(self)[1], " must implement za_family()."),

		za_description = function() stop(class(self)[1], " must implement za_description()."),

		predictors_df = function(){
			data.frame(w = private$w)
		},

		reduce_design_matrix_preserving_treatment = function(X_full){
			qr_X = qr(X_full)
			target_rank = qr_X$rank
			required = c(1L, 2L)
			candidate_order = c(required, setdiff(qr_X$pivot, required))
			keep = integer(0)

			for (j in candidate_order){
				trial_keep = c(keep, j)
				trial_rank = qr(X_full[, trial_keep, drop = FALSE])$rank
				if (trial_rank > length(keep)){
					keep = trial_keep
				}
				if (length(keep) >= target_rank){
					break
				}
			}

			keep = sort(unique(keep))
			if (!(2L %in% keep)){
				return(NULL)
			}
			X_full[, keep, drop = FALSE]
		},

		build_formula = function(dat){
			fixed_terms = setdiff(colnames(dat), "y")
			stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))
		},

		build_zi_formula = function(dat){
			fixed_terms = setdiff(colnames(dat), "y")
			stats::as.formula(paste("~", paste(fixed_terms, collapse = " + ")))
		},

		fit_zero_augmented_model = function(dat){
			formula_cond = private$build_formula(dat)
			formula_zi = private$build_zi_formula(dat)

			glmm_control = if (private$num_cores > 1) {
				glmmTMB::glmmTMBControl(parallel = 1L)
			} else {
				glmmTMB::glmmTMBControl()
			}

			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = formula_zi,
						family = private$za_family(),
						data = dat,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)

			if (!is.null(mod)) return(mod)
			if (ncol(dat) <= 2L) return(NULL)

			dat_fallback = dat[, c("y", "w"), drop = FALSE]
			tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						y ~ w,
						ziformula = ~ w,
						family = private$za_family(),
						data = dat_fallback,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = cbind(1, private$w, as.matrix(private$predictors_df()[, setdiff(colnames(private$predictors_df()), "w"), drop = FALSE]))
			colnames(X_full)[1:2] = c("(Intercept)", "w")
			X_reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			if (is.null(X_reduced)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			pred_df = as.data.frame(X_reduced[, -1, drop = FALSE])
			colnames(pred_df)[1] = "w"
			dat = data.frame(y = private$y, pred_df)

			mod = private$fit_zero_augmented_model(dat)
			if (is.null(mod)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table))){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(coef_table["w", "Estimate"])
			se = as.numeric(coef_table["w", "Std. Error"])
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop(private$za_description(), ": non-finite standard error (possible convergence failure or separation in zero model).")
			}
		}
	)
)
