#' Univariate Conditional Proportional-Odds Combined Inference for KK Designs
#'
#' Fits a conditional proportional-odds model for ordinal responses under a KK
#' matching-on-the-fly design, including both matched pairs and reservoir subjects.
#' Each matched pair forms a stratum across thresholds. Each reservoir subject
#' is treated as its own stratum. This is analogous to a stratified proportional odds
#' model where the strata are the pairs and the individual reservoir subjects.
#'
#' @details
#' This estimator uses the conditional likelihood approach by expanding the ordinal
#' response into binary threshold indicators. Matched pairs contribute to the
#' conditional likelihood. Note that singleton strata (reservoir subjects) do not
#' contribute to the conditional logit likelihood directly unless combined with other
#' information or if the model is viewed as a stratified model. In this implementation,
#' we follow the standard conditional logit approach where only strata with variation
#' in both treatment and response contribute.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$new(n = nrow(x_dat), response_type = "ordinal",
#' verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUnivKKCondPropOddsCombinedRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalUnivKKCondPropOddsCombinedRegr = R6::R6Class(
	"InferenceOrdinalUnivKKCondPropOddsCombinedRegr",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize a univariate conditional proportional-odds combined inference object.
		#' @param	des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param	num_cores			Number of CPU cores.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design.")
			}
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha Significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta Null value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("Testing non-zero delta is not yet implemented.")
			}
		}
	),

	private = list(
		max_abs_reasonable_coef = 1e4,

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL)) 
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: conditional logit ---
			private$clogit_for_matched_pairs()
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: proportional odds ---
			private$ord_for_reservoir()
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0

			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		clogit_for_matched_pairs = function(){
			m_vec = private$m
			m = max(m_vec, na.rm = TRUE)
			if (m == 0) return(invisible(NULL))

			y_ord = as.integer(factor(private$y, ordered = TRUE))
			K = max(y_ord)
			if (K < 2L) return(invisible(NULL))

			n_thresholds = K - 1L
			matched_idx = which(m_vec > 0)
			y_m = y_ord[matched_idx]
			w_m = private$w[matched_idx]
			strata_m = m_vec[matched_idx]
			n_m = length(matched_idx)

			y_stack = integer(n_m * n_thresholds)
			w_stack = integer(n_m * n_thresholds)
			strata_stack = integer(n_m * n_thresholds)

			for (k in seq_len(n_thresholds)){
				idx = ((k - 1L) * n_m + 1L):(k * n_m)
				y_stack[idx] = as.integer(y_m > k)
				w_stack[idx] = w_m
				strata_stack[idx] = strata_m + (k - 1L) * m
			}

			mod = clogit_helper(y_stack, data.frame(), w_stack, strata_stack)
			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(mod$b[1])
			ssq  = as.numeric(mod$ssq_b_j)
			private$cached_values$beta_T_matched     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(ssq) && ssq > 0) ssq else NA_real_
		},

		ord_for_reservoir = function(){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			if (length(y_r) < 2 || length(unique(y_r)) < 2 || length(unique(w_r)) < 2) return(invisible(NULL))

			Xmm = matrix(w_r, ncol = 1)
			mod = tryCatch(
				fast_ordinal_regression_with_var_cpp(X = Xmm, y = as.numeric(y_r)),
				error = function(e) NULL
			)
			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(mod$b[1])
			ssq  = as.numeric(mod$ssq_b_2)
			private$cached_values$beta_T_reservoir     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(ssq) && ssq > 0) ssq else NA_real_
		}
	)
)
