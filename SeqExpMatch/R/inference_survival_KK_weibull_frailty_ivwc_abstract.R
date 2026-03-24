# Abstract class for Weibull Frailty / Standard Weibull Compound Inference
#
# @description
# This class implements a compound estimator for KK matching-on-the-fly designs with
# survival responses using the Weibull AFT model. For matched pairs, it uses a frailty
# model (each pair is a cluster) to estimate the log-time ratio (AFT scale). For
# reservoir subjects, it uses standard Weibull AFT regression. The two estimates
# are combined via a variance-weighted linear combination.
#
# @details
# The matched-pairs frailty estimator differs between the univariate and multivariate
# subclasses due to computational constraints:
#
# \strong{Univariate} (\code{include_covariates() == FALSE}): uses \code{parfm::parfm()}
# to fit a fully parametric Weibull gamma-frailty model with only the treatment
# indicator. The PH coefficient and Weibull shape \eqn{\rho} are both taken from
# \code{parfm}, and the AFT log-time ratio is recovered as
# \eqn{\hat\alpha_\mathrm{AFT} = -\hat\beta_\mathrm{PH} / \hat\rho}.
# This is fast because there is only one fixed-effect parameter.
#
# \strong{Multivariate} (\code{include_covariates() == TRUE}): \code{parfm} becomes
# prohibitively slow when many covariates are present (its EM algorithm is
# \eqn{O(p^2)} per iteration and convergence degrades sharply with dimension).
# Instead, the treatment PH coefficient \eqn{\hat\beta_\mathrm{PH}} and its SE are
# obtained from \code{survival::coxph()} with a \code{frailty(cluster,
# distribution = "gamma")} term, which is a highly optimised C implementation.
# The Weibull shape \eqn{\hat\rho} is estimated separately from a standard
# \code{survival::survreg()} fit (no frailty) on the same matched-pair data.
# The AFT log-time ratio is then
# \eqn{\hat\alpha_\mathrm{AFT} = -\hat\beta_\mathrm{PH} / \hat\rho},
# and the SE is \eqn{\widehat{\mathrm{SE}}_\mathrm{AFT} = \widehat{\mathrm{SE}}_\mathrm{PH} / \hat\rho}
# (delta method treating \eqn{\hat\rho} as fixed). Because \code{coxph} is
# semi-parametric, \eqn{\hat\rho} from the auxiliary \code{survreg} fit is used
# only for scale conversion and does not affect the frailty correction itself.
#
# The \pkg{parfm} package is required only for the univariate subclass.
#
# @keywords internal
SeqDesignInferenceAbstractKKWeibullFrailtyIVWC = R6::R6Class("SeqDesignInferenceAbstractKKWeibullFrailtyIVWC",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			if (!requireNamespace("parfm", quietly = TRUE)) {
				stop("Package 'parfm' is required for ", class(self)[1], ". Please install it.")
			}
			assertResponseType(seq_des_obj$get_response_type(), "survival")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		# @description
		# Returns the estimated treatment effect (log-time ratio).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes the asymptotic confidence interval.
		# @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes the asymptotic p-value.
		# @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("Testing non-zero delta is not yet implemented for this class.")
			}
		}
	),

	private = list(

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		# Overridden to use the faster coxph+survreg approximation during randomization inference
		# for the univariate case, avoiding the heavy parfm overhead.
		compute_treatment_estimate_during_randomization_inference = function(){
			if (private$include_covariates()){
				# Multivariate already uses coxph+survreg which is fast
				return(self$compute_treatment_estimate())
			}

			# For univariate, use the fast approximation for the matched-pairs component
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			beta_m = NA_real_
			ssq_m  = NA_real_
			if (m > 0){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(0L, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_matched = which(m_vec > 0)
				y_m       = private$y[i_matched]
				dead_m    = private$dead[i_matched]
				w_m       = private$w[i_matched]
				strata_m  = m_vec[i_matched]

				strata_with_events = unique(strata_m[dead_m == 1])
				if (length(strata_with_events) > 0){
					i_valid = which(strata_m %in% strata_with_events)
					dat = data.frame(
						time = y_m[i_valid],
						event = dead_m[i_valid],
						w = w_m[i_valid],
						cluster = factor(strata_m[i_valid])
					)
					res = private$fit_fast_approx(dat)
					if (!is.null(res)){
						beta_m = res$beta_aft
						ssq_m  = res$ssq_aft
					}
				}
			}

			beta_r = NA_real_
			ssq_r  = NA_real_
			if (nRT > 0 && nRC > 0){
				# Reservoir uses standard Weibull which is already fast
				private$weibull_for_reservoir()
				beta_r = private$cached_values$beta_T_reservoir
				ssq_r  = private$cached_values$ssq_beta_T_reservoir
			}

			m_ok = is.finite(beta_m) && is.finite(ssq_m) && ssq_m > 0
			r_ok = is.finite(beta_r) && is.finite(ssq_r) && ssq_r > 0

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		},

		# High-speed bootstrap implementation using the coxph+survreg approximation.
		compute_fast_bootstrap_distr = function(B, i_reservoir, n_reservoir, m, y, w, m_vec) {
			dead = private$dead
			X = if (private$include_covariates()) as.matrix(private$X) else NULL
			if (!is.null(X) && is.null(colnames(X))) colnames(X) = paste0("x", seq_len(ncol(X)))
			cov_str = if (!is.null(X)) paste(colnames(X), collapse = " + ") else NULL

			cores_to_use = private$num_cores

			beta_hat_T_bs = unlist(parallel::mclapply(1:B, function(b) {
				# --- Matched pairs component ---
				beta_m = NA_real_
				ssq_m  = NA_real_
				if (m > 0) {
					pairs_b = sample(1:m, m, replace = TRUE)
					i_match_b = integer(0)
					cluster_b = integer(0)
					for (new_id in 1:m) {
						orig_id = pairs_b[new_id]
						idx = which(m_vec == orig_id)
						i_match_b = c(i_match_b, idx)
						cluster_b = c(cluster_b, rep(new_id, length(idx)))
					}
					dat_m = data.frame(
						time = y[i_match_b],
						event = dead[i_match_b],
						w = w[i_match_b],
						cluster = factor(cluster_b)
					)
					if (!is.null(X)) dat_m = cbind(dat_m, X[i_match_b, , drop = FALSE])
					res_m = private$fit_fast_approx(dat_m, cov_str)
					if (!is.null(res_m)) {
						beta_m = res_m$beta_aft
						ssq_m  = res_m$ssq_aft
					}
				}

				# --- Reservoir component ---
				beta_r = NA_real_
				ssq_r  = NA_real_
				if (n_reservoir > 0) {
					i_res_b = sample(i_reservoir, n_reservoir, replace = TRUE)
					X_r = if (!is.null(X)) cbind(w = w[i_res_b], X[i_res_b, , drop = FALSE]) else matrix(w[i_res_b], ncol = 1, dimnames = list(NULL, "w"))
					mod_r = tryCatch(fast_weibull_regression(y[i_res_b], dead[i_res_b], X_r), error = function(e) NULL)
					if (!is.null(mod_r)) {
						beta_r = as.numeric(mod_r$coefficients["w"])
						ssq_r  = if (is.matrix(mod_r$vcov) && "w" %in% colnames(mod_r$vcov)) mod_r$vcov["w", "w"] else NA_real_
					}
				}

				# --- Pooling ---
				m_ok = is.finite(beta_m) && is.finite(ssq_m) && ssq_m > 0
				r_ok = is.finite(beta_r) && is.finite(ssq_r) && ssq_r > 0

				if (m_ok && r_ok) {
					w_star = ssq_r / (ssq_r + ssq_m)
					w_star * beta_m + (1 - w_star) * beta_r
				} else if (m_ok) {
					beta_m
				} else if (r_ok) {
					beta_r
				} else {
					NA_real_
				}
			}, mc.cores = cores_to_use))
			beta_hat_T_bs
		},

		# Fast approximation using semi-parametric coxph with frailty + univariate survreg for rho
		fit_fast_approx = function(dat, cov_str = NULL){
			cox_formula = if (is.null(cov_str)){
				survival::Surv(time, event) ~ w + survival::frailty(cluster, distribution = "gamma")
			} else {
				as.formula(paste("survival::Surv(time, event) ~ w +", cov_str,
				                 "+ survival::frailty(cluster, distribution = 'gamma')"))
			}

			cox_mod = tryCatch(survival::coxph(cox_formula, data = dat), error = function(e) NULL)
			if (is.null(cox_mod)) return(NULL)

			cox_coef = tryCatch(summary(cox_mod)$coefficients, error = function(e) NULL)
			if (is.null(cox_coef) || !("w" %in% rownames(cox_coef))) return(NULL)
			beta_ph = cox_coef["w", "coef"]
			se_ph   = cox_coef["w", "se(coef)"]

			survreg_formula = if (is.null(cov_str)){
				survival::Surv(time, event) ~ w
			} else {
				as.formula(paste("survival::Surv(time, event) ~ w +", cov_str))
			}

			survreg_mod = tryCatch(survival::survreg(survreg_formula, data = dat, dist = "weibull"), error = function(e) NULL)
			if (is.null(survreg_mod)) return(NULL)

			rho = 1 / survreg_mod$scale
			list(beta_aft = -beta_ph / rho, ssq_aft = (se_ph / rho)^2)
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			# Recompute KKstats if cache was cleared (e.g., after y transformation for rand CI)
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: Weibull Frailty ---
			if (m > 0){
				private$weibull_frailty_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: Standard Weibull (AFT) ---
			if (nRT > 0 && nRC > 0){
				private$weibull_for_reservoir()
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0

			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
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

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Weibull Frailty compound estimator: could not compute a finite standard error.")
			}
		},

		weibull_frailty_for_matched_pairs = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0)
			y_m       = private$y[i_matched]
			dead_m    = private$dead[i_matched]
			w_m       = private$w[i_matched]
			strata_m  = m_vec[i_matched]

			# Filter strata that have no events (provide no information)
			strata_with_events = unique(strata_m[dead_m == 1])
			if (length(strata_with_events) == 0) return(invisible(NULL))

			i_valid = which(strata_m %in% strata_with_events)

			dat = data.frame(
				time = y_m[i_valid],
				event = dead_m[i_valid],
				w = w_m[i_valid],
				cluster = factor(strata_m[i_valid])
			)

			if (private$include_covariates()){
				# Multivariate: coxph with gamma frailty (optimised C implementation) for the
				# PH treatment coefficient; survreg for Weibull shape rho; convert via
				# alpha_AFT = -beta_PH / rho,  SE_AFT = SE_PH / rho  (delta method, rho fixed).
				X_m = as.matrix(private$X[i_matched[i_valid], , drop = FALSE])
				colnames(X_m) = paste0("x", 1:ncol(X_m))
				dat = cbind(dat, X_m)
				cov_str = paste(colnames(X_m), collapse = " + ")

				res = private$fit_fast_approx(dat, cov_str)
				if (is.null(res)) return(invisible(NULL))
				alpha_aft     = res$beta_aft
				ssq_alpha_aft = res$ssq_aft
			} else {
				# Univariate: parfm parametric Weibull gamma-frailty (fast when p = 1).
				mod = tryCatch({
					parfm::parfm(
						survival::Surv(time, event) ~ w,
						cluster = "cluster",
						data = dat,
						dist = "weibull",
						frailty = "gamma"
					)
				}, error = function(e) NULL)

				if (is.null(mod)) return(invisible(NULL))

				# parfm returns PH coefficients; convert to AFT: alpha_AFT = -beta_PH / rho.
				rho           = as.numeric(mod["rho", "ESTIMATE"])
				beta_ph       = as.numeric(mod["w",   "ESTIMATE"])
				se_ph         = as.numeric(mod["w",   "SE"])
				alpha_aft     = -beta_ph / rho
				ssq_alpha_aft = (se_ph / rho)^2
			}

			private$cached_values$beta_T_matched     = if (is.finite(alpha_aft)) alpha_aft else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(ssq_alpha_aft) && ssq_alpha_aft > 0) ssq_alpha_aft else NA_real_
		},

		weibull_for_reservoir = function(){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			m_vec_safe = private$m
			if (is.null(m_vec_safe)) m_vec_safe = rep(0L, private$n)
			m_vec_safe[is.na(m_vec_safe)] = 0L
			dead_r = private$dead[m_vec_safe == 0]
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)

			X_full = matrix(w_r, ncol = 1)
			colnames(X_full) = "w"

			if (private$include_covariates()){
				X_full = cbind(X_full, X_r)
			}

			# fast_weibull_regression already uses survreg (AFT units).
			mod = tryCatch(
				fast_weibull_regression(y_r, dead_r, X_full),
				error = function(e) NULL
			)

			# Fallback if multivariate failed
			if (is.null(mod) && private$include_covariates()){
				mod = tryCatch(
					fast_weibull_regression(y_r, dead_r, matrix(w_r, ncol = 1, dimnames = list(NULL, "w"))),
					error = function(e) NULL
				)
			}

			if (is.null(mod)) return(invisible(NULL))

			alpha_aft = as.numeric(mod$coefficients["w"])
			ssq_alpha_aft = if (is.matrix(mod$vcov) && "w" %in% colnames(mod$vcov)) mod$vcov["w", "w"] else NA_real_

			private$cached_values$beta_T_reservoir     = if (is.finite(alpha_aft)) alpha_aft else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(ssq_alpha_aft) && ssq_alpha_aft > 0) ssq_alpha_aft else NA_real_
		}
	)
)
