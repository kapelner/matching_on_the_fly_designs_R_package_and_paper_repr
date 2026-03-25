#' Hurdle Negative Binomial Inference for Count Responses
#'
#' @description
#' Internal base class for non-KK hurdle negative binomial regression models. The
#' hurdle indicator model is fit on all subjects, and the count component is fit on
#' the positive-count subjects using a zero-truncated negative binomial model
#' optimized in C++ via L-BFGS. The reported treatment effect is the treatment
#' coefficient from the conditional count component, on the log-rate scale.
#'
#' @keywords internal
#' @noRd
DesignInferenceCountHurdleNegBinAbstract = R6::R6Class("DesignInferenceCountHurdleNegBinAbstract",
	inherit = DesignInference,
	public = list(

		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "count")
			super$initialize(seq_des_obj, num_cores, verbose)
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
		hurdle_description = function() stop(class(self)[1], " must implement hurdle_description()."),

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
				return(list(X = NULL, keep = keep, j_treat = NA_integer_))
			}

			list(
				X = X_full[, keep, drop = FALSE],
				keep = keep,
				j_treat = match(2L, keep)
			)
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			X_full = cbind(1, private$w, as.matrix(private$predictors_df()[, setdiff(colnames(private$predictors_df()), "w"), drop = FALSE]))
			colnames(X_full)[1:2] = c("(Intercept)", "treatment")
			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			mod = tryCatch(
				fast_hurdle_negbin_with_var_cpp(X_fit, private$y, j = j_treat),
				error = function(e) NULL
			)
			if (is.null(mod)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			coef_hat = as.numeric(mod$b)
			ssq_b_2 = as.numeric(mod$ssq_b_j)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat)) || !is.finite(ssq_b_2) || ssq_b_2 < 0){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			b_full = rep(NA_real_, ncol(X_full))
			b_full[reduced$keep] = coef_hat
			names(b_full) = colnames(X_full)

			private$cached_values$beta_hat_T = coef_hat[j_treat]
			private$cached_values$s_beta_hat_T = sqrt(ssq_b_2)
			private$cached_values$is_z = TRUE
			private$cached_values$full_coefficients = b_full
			private$cached_values$theta_hat = as.numeric(mod$theta_hat)
			private$cached_values$hurdle_coefficients = mod$hurdle_b
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop(private$hurdle_description(), ": non-finite standard error (possible convergence failure in the truncated count component).")
			}
		}
	)
)
