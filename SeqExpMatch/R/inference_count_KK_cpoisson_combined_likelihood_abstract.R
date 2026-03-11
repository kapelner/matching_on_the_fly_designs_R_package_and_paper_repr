# Abstract class for Conditional Poisson Combined-Likelihood Compound Inference
#
# @description
# Fits a single joint likelihood over all KK design data for count responses.
# The matched-pair component uses the conditional Poisson likelihood; conditioning
# on each pair's total count n_k = yT_k + yC_k reduces it to a weighted binomial
# logistic form: yT_k | n_k ~ Binomial(n_k, p_k) with logit(p_k) = beta_T + X_diff_k' beta_xs.
# The reservoir component uses the marginal Poisson likelihood.
# Both components share beta_T and beta_xs; beta_0 (intercept) is only identified
# from reservoir rows and is profiled out in the pairs-only case.
#
# Because the two components use different link functions (logit vs log), they
# cannot be stacked into a single GLM.  Instead, the combined log-likelihood
#   L_total = L_cond_Poisson(pairs) + L_Poisson(reservoir)
# is maximized jointly via BFGS with an analytic gradient.
# The standard error for beta_T is extracted from the numerical Hessian at the MLE.
#
# Degenerate cases:
#   Pairs only     : reduces to weighted logistic (glm.fit), no BFGS needed.
#   Reservoir only : reduces to Poisson regression (glm.fit), no BFGS needed.
#
# @keywords internal
SeqDesignInferenceAbstractKKCPoissonCombinedLikelihood = R6::R6Class("SeqDesignInferenceAbstractKKCPoissonCombinedLikelihood",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "count")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
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
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared_combined_likelihood()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Returns a 2-sided p-value for H0: beta_T = delta.
		# @param delta Null value; default 0.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared_combined_likelihood()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("CPoisson combined-likelihood: could not compute a finite standard error (possible convergence failure or insufficient data).")
			}
		},

		# Maximize L_cond_Poisson(pairs) + L_Poisson(reservoir) jointly.
		#
		# Pair contribution (conditional Poisson = weighted binomial logistic):
		#   eta_k = beta_T + X_diff_k' beta_xs
		#   ll_pairs = sum(yT_k * eta_k - n_k * log(1 + exp(eta_k)))
		#
		# Reservoir contribution (marginal Poisson):
		#   eta_i = beta_0 + w_i * beta_T + X_i' beta_xs
		#   ll_res = sum(y_i * eta_i - exp(eta_i))
		#
		# Parameter layout: [beta_0, beta_T, beta_xs (p)] when reservoir present;
		#                   [beta_T, beta_xs (p)]           when pairs-only.
		# The combined case is handled by fast_cpoisson_combined_with_var_cpp
		# (Newton's method with analytic Fisher-information Hessian).
		shared_combined_likelihood = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			p             = if (private$include_covariates()) ncol(private$X) else 0L
			has_reservoir = nRT > 0 && nRC > 0

			# ---- Pair data (conditional Poisson, zero-total pairs discarded) ----
			yT_v     = numeric(0)
			n_k_v    = numeric(0)
			X_diff_v = matrix(nrow = 0L, ncol = p)

			if (m > 0){
				yT    = KKstats$yTs_matched
				yC    = KKstats$yCs_matched
				n_k   = yT + yC
				valid = which(n_k > 0)
				if (length(valid) > 0){
					yT_v  = yT[valid]
					n_k_v = n_k[valid]
					if (p > 0L) X_diff_v = as.matrix(KKstats$X_matched_diffs[valid, , drop = FALSE])
				}
			}
			has_pairs = length(yT_v) > 0

			# ---- Reservoir data (marginal Poisson) ----
			y_r_v = numeric(0)
			w_r_v = numeric(0)
			X_r_v = matrix(nrow = 0L, ncol = p)
			if (has_reservoir){
				y_r_v = KKstats$y_reservoir
				w_r_v = KKstats$w_reservoir
				if (p > 0L) X_r_v = as.matrix(KKstats$X_reservoir)
			}

			if (!has_pairs && !has_reservoir){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- Case 1: pairs + reservoir — joint BFGS optimization --------
			if (has_pairs && has_reservoir){
				# Combined: Newton's method in C++ (fast_cpoisson_combined_with_var_cpp)
				# params layout: [beta_0, beta_T, beta_xs]; beta_T is at 1-based index 2
				mod = tryCatch(
					fast_cpoisson_combined_with_var_cpp(
						as.double(yT_v), as.double(n_k_v), X_diff_v,
						as.double(y_r_v), as.double(w_r_v), X_r_v
					),
					error = function(e) NULL
				)
				if (is.null(mod) || !mod$converged || !is.finite(mod$ssq_b_j) || mod$ssq_b_j <= 0){
					private$cached_values$beta_hat_T   = NA_real_
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}
				private$cached_values$beta_hat_T   = as.numeric(mod$b[2L])
				private$cached_values$s_beta_hat_T = sqrt(mod$ssq_b_j)
				private$cached_values$is_z         = TRUE

			# ---- Case 2: pairs only — weighted logistic (glm.fit) ---------------
			} else if (has_pairs){
				y_prop = yT_v / n_k_v
				Xmm    = if (p > 0L) cbind(1, X_diff_v) else matrix(1, nrow = length(yT_v), ncol = 1)
				mod = tryCatch({
					g_mod = stats::glm.fit(x = Xmm, y = y_prop, weights = n_k_v, family = stats::binomial())
					R_mat = qr.R(g_mod$qr)
					Rinv  = backsolve(R_mat, diag(ncol(R_mat)))
					vcov  = Rinv %*% t(Rinv)
					list(b = g_mod$coefficients, ssq_b_1 = vcov[1, 1])
				}, error = function(e) NULL)
				if (is.null(mod)){
					private$cached_values$beta_hat_T   = NA_real_
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}
				private$cached_values$beta_hat_T   = as.numeric(mod$b[1L])
				private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_1))
				private$cached_values$is_z         = TRUE

			# ---- Case 3: reservoir only — Poisson regression (glm.fit) ----------
			} else {
				j_beta_T = 2L
				X_full   = if (p > 0L) cbind(1, w_r_v, X_r_v) else cbind(1, w_r_v)
				qr_X = qr(X_full)
				if (qr_X$rank < ncol(X_full)){
					keep     = qr_X$pivot[seq_len(qr_X$rank)]
					if (!(j_beta_T %in% keep)) keep[qr_X$rank] = j_beta_T
					keep     = sort(keep)
					X_full   = X_full[, keep, drop = FALSE]
					j_beta_T = which(keep == j_beta_T)
				}
				mod = tryCatch({
					g_mod  = stats::glm.fit(x = X_full, y = y_r_v, family = stats::poisson())
					mu_hat = g_mod$fitted.values
					XtWX   = crossprod(X_full * sqrt(mu_hat))   # Fisher information X'WX, W=diag(mu)
					vcov   = solve(XtWX)
					list(b = g_mod$coefficients, ssq_b_j = vcov[j_beta_T, j_beta_T])
				}, error = function(e) NULL)
				if (is.null(mod)){
					private$cached_values$beta_hat_T   = NA_real_
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}
				private$cached_values$beta_hat_T   = as.numeric(mod$b[j_beta_T])
				private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
				private$cached_values$is_z         = TRUE
			}

			invisible(NULL)
		}
	)
)
