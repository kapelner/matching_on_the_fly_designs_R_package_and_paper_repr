#' Abstract Quantile Regression Combined-Likelihood Compound Estimator for KK Designs
#'
#' Fits a single joint quantile regression over all KK design data by stacking
#' matched-pair differences and reservoir observations into one design matrix.
#'
#' Column layout of X_stack: [beta_0 | beta_T | beta_xs (p cols)]
#' Pair rows:      [0 | 1   | Xd_k] → Q_tau(yd_k) = beta_T + Xd_k' beta_xs
#' Reservoir rows: [1 | w_i | X_i]  → Q_tau(y_i)  = beta_0 + w_i*beta_T + X_i'*beta_xs
#' Fitting a single rq() on the stacked dataset minimises the combined check-function loss.
#'
#' Special cases:
#' Pairs only:     beta_0 column is all-zero and dropped; layout [beta_T | beta_xs].
#' Reservoir only: standard quantile regression layout [beta_0 | beta_T | beta_xs].
#'
#' Standard errors use Powell's "nid" sandwich estimator, falling back to "iid".
#'
#' @keywords internal
InferenceAbstractKKQuantileRegrOneLik = R6::R6Class("InferenceAbstractKKQuantileRegrOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractQuantileRandCI,
	public = list(

		#' @description
		#' Initialize KK quantile-regression combined-likelihood inference.
		#' @param des_obj A completed KK design object.
		#' @param tau The quantile level for regression, strictly between 0 and 1. The default
		#'   \code{tau = 0.5} estimates the median treatment effect. Values of exactly 0 or 1
		#'   are excluded because quantile regression is undefined at the boundary (the check
		#'   function \eqn{\rho_\tau(u) = u(\tau - \mathbf{1}_{u < 0})} degenerates there);
		#'   the bound is enforced as \code{(.Machine$double.eps, 1 - .Machine$double.eps)},
		#'   i.e. the smallest representable positive number away from 0 and 1.
		#' @param transform_y_fn Optional response transformation.
		#' @param verbose Whether to print progress messages.
		#' @return A new inference object.
		initialize = function(des_obj, tau = 0.5, transform_y_fn = identity,  verbose = FALSE){
			if (should_run_asserts()) {
				assertNumeric(tau, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			}
			if (should_run_asserts()) {
				if (!check_package_installed("quantreg")) {
					stop("Package 'quantreg' is required. Please install it with install.packages(\"quantreg\").")
				}
			}
			private$tau = tau
			private$transform_y_fn_list = list(fn = transform_y_fn)
			super$initialize(des_obj, verbose)
			if (private$is_KK){
				private$m = des_obj$.__enclos_env__$private$m
				private$compute_basic_match_data()
			}
		},

		#' @description
		#' Compute the quantile-regression treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		#' @return The treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute an asymptotic confidence interval.
		#' @param alpha Significance level.
		#' @return A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute an asymptotic two-sided p-value for the treatment effect.
		#' @param delta Null treatment effect value.
		#' @return A two-sided p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		tau = NULL,
		transform_y_fn_list = NULL,  # list(fn = ...) wrapping avoids R6 treating function as a locked method
		m = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		compute_basic_match_data = function(){
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv     = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		# Fit the combined check-function loss over matched-pair differences and
		# reservoir observations with SHARED covariate effects beta_xs.
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			tau = private$tau
			fn  = private$transform_y_fn_list$fn

			has_reservoir = nRT > 0 && nRC > 0

			y_stack  = NULL
			X_stack  = NULL
			j_beta_T = 2L

			if (m > 0){
				yd = fn(KKstats$yTs_matched) - fn(KKstats$yCs_matched)

				if (has_reservoir){
					# Combined case: pair rows and reservoir rows must share the same
					# covariate columns. Use the full-width pair differences from the
					# shared C++ preprocessing rather than the reduced pair-only matrix.
					Xd = as.matrix(KKstats$X_matched_diffs_full)
					p   = ncol(Xd)
					y_r = fn(KKstats$y_reservoir)
					w_r = KKstats$w_reservoir
					X_r = as.matrix(KKstats$X_reservoir)
					X_pairs = if (p > 0) cbind(0, 1, Xd) else matrix(c(0, 1), nrow = m, ncol = 2, byrow = TRUE)
					X_res   = if (p > 0) cbind(1, w_r, X_r) else cbind(1, w_r)
					X_stack  = rbind(X_pairs, X_res)
					y_stack  = c(yd, y_r)
					j_beta_T = 2L
				} else {
					# Pairs only: drop all-zero beta_0 column; intercept = beta_T.
					Xd = as.matrix(KKstats$X_matched_diffs)
					p  = ncol(Xd)
					X_stack  = if (p > 0) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1)
					y_stack  = yd
					j_beta_T = 1L
				}
			} else if (has_reservoir){
				y_r = fn(KKstats$y_reservoir)
				w_r = KKstats$w_reservoir
				X_r = as.matrix(KKstats$X_reservoir)
				p   = ncol(X_r)
				X_stack  = if (p > 0) cbind(1, w_r, X_r) else cbind(1, w_r)
				y_stack  = y_r
				j_beta_T = 2L
			}

			if (is.null(X_stack)){
				private$cached_values$beta_hat_T   = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# QR-reduce to full rank, preserving beta_T column
			reduced = qr_reduce_preserve_cols_cpp(X_stack, j_beta_T)
			X_stack  = reduced$X_reduced
			j_beta_T = match(j_beta_T, reduced$keep)

			n_total  = nrow(X_stack)
			n_params = ncol(X_stack)
			if (n_total <= n_params){
				private$cached_values$beta_hat_T   = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			cn = paste0("x", seq_len(n_params))
			cn[j_beta_T] = "trt__"
			colnames(X_stack) = cn
			dat = as.data.frame(X_stack)
			dat$y_stack__ = y_stack

			fit = tryCatch(
				suppressWarnings(quantreg::rq(y_stack__ ~ . - 1, tau = tau, data = dat)),
				error = function(e) NULL
			)
			if (is.null(fit)){
				private$cached_values$beta_hat_T   = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			beta = tryCatch(coef(fit)[["trt__"]], error = function(e) NA_real_)
			private$cached_values$beta_hat_T   = if (is.finite(beta)) beta else NA_real_
			if (!estimate_only) {
				se = private$extract_se_from_rq(fit, "trt__")
				private$cached_values$s_beta_hat_T = if (!is.na(se)) se else NA_real_
			}
			private$cached_values$is_z         = TRUE
			invisible(NULL)
		},

		# Helper: extract SE from rq fit by coefficient name, trying "nid" then "iid".
		# SEs above 1e6 are treated as invalid (the "nid" sparsity estimator can return
		# astronomically large but finite values when the density at the quantile is near
		# zero, which bypasses the usual !is.finite() || <= 0 guard in callers).
		extract_se_from_rq = function(fit, coef_name){
			.extract_se_from_rq_fit(fit, coef_name)
		}
	)
)
