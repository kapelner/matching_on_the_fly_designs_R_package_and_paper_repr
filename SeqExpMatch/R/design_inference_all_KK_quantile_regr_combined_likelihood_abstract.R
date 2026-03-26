# Abstract Quantile Regression Combined-Likelihood Compound Estimator for KK Designs
#
# @description
# Fits a single joint quantile regression over all KK design data by stacking
# matched-pair differences and reservoir observations into one design matrix.
#
# Column layout of X_stack: [beta_0 | beta_T | beta_xs (p cols)]
#   Pair rows:      [0 | 1   | Xd_k] → Q_tau(yd_k) = beta_T + Xd_k' beta_xs
#   Reservoir rows: [1 | w_i | X_i]  → Q_tau(y_i)  = beta_0 + w_i*beta_T + X_i'*beta_xs
# Fitting a single rq() on the stacked dataset minimises the combined check-function loss.
#
# Special cases:
#   Pairs only:     beta_0 column is all-zero and dropped; layout [beta_T | beta_xs].
#   Reservoir only: standard quantile regression layout [beta_0 | beta_T | beta_xs].
#
# Standard errors use Powell's "nid" sandwich estimator, falling back to "iid".
#
# @keywords internal
DesignInferenceAbstractKKQuantileRegrCombinedLikelihood = R6::R6Class("DesignInferenceAbstractKKQuantileRegrCombinedLikelihood",
	inherit = DesignInferenceAbstractQuantileRandCI,
	public = list(

		# @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		# @param tau				The quantile level for regression, strictly between 0 and 1.
		# @param transform_y_fn	A function applied to y values before quantile regression.
		# @param num_cores			The number of CPU cores to use to parallelize sampling.
		# @param verbose			A flag indicating whether messages should be displayed. Default is FALSE.
		initialize = function(seq_des_obj, tau = 0.5, transform_y_fn = identity, num_cores = 1, verbose = FALSE){
			assertNumeric(tau, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			if (!requireNamespace("quantreg", quietly = TRUE)) {
				stop("Package 'quantreg' is required. Please install it with install.packages(\"quantreg\").")
			}
			private$tau = tau
			private$transform_y_fn_list = list(fn = transform_y_fn)
			super$initialize(seq_des_obj, num_cores, verbose)
			if (private$is_KK){
				private$m = seq_des_obj$.__enclos_env__$private$m
				private$compute_basic_match_data()
			}
		},

		# @description
		# Returns the combined-likelihood estimate of the treatment effect.
		compute_treatment_estimate = function(){
			private$shared_combined_likelihood()
			private$cached_values$beta_hat_T
		},

		# @description
		# Returns a 1 - alpha confidence interval for beta_T.
		# @param alpha Significance level; default 0.05 gives a 95% CI.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared_combined_likelihood()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Returns a 2-sided p-value for H0: beta_T = delta.
		# @param delta Null value; default 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared_combined_likelihood()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		tau = NULL,
		transform_y_fn_list = NULL,  # list(fn = ...) wrapping avoids R6 treating function as a locked method
		m = NULL,

		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			private$cached_values$KKstats = .compute_kk_basic_match_data(
				X = private$X,
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Quantile regression combined-likelihood: could not compute a finite standard error.")
			}
		},

		# Fit the combined check-function loss over matched-pair differences and
		# reservoir observations with SHARED covariate effects beta_xs.
		#
		# Column layout of X_stack: [beta_0 | beta_T | beta_xs (p cols)]
		#   col 1       : beta_0   (0 for pair rows; 1 for reservoir rows)
		#   col 2       : beta_T   (1 for pair rows; w_i for reservoir rows)
		#   cols 3..p+2 : beta_xs  (Xd_k for pair rows; X_i for reservoir rows)
		#
		# Pairs-only case: beta_0 column is all-zero and dropped; layout [beta_T | beta_xs].
		# Reservoir-only case: standard layout [beta_0 | beta_T | beta_xs].
		shared_combined_likelihood = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

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
				private$cached_values$s_beta_hat_T = NA_real_
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
				private$cached_values$s_beta_hat_T = NA_real_
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
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			beta = tryCatch(coef(fit)[["trt__"]], error = function(e) NA_real_)
			se   = private$extract_se_from_rq(fit, "trt__")

			private$cached_values$beta_hat_T   = if (is.finite(beta)) beta else NA_real_
			private$cached_values$s_beta_hat_T = if (!is.na(se)) se else NA_real_
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
