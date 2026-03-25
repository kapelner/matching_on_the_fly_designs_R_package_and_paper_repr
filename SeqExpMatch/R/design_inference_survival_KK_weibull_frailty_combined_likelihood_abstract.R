# Abstract class for Weibull Frailty Combined-Likelihood Compound Inference
#
# @description
# Fits a single joint parametric Weibull frailty model over all KK design data
# for survival (or count) responses. Matched subjects belong to their pair cluster
# (frailty integrates away the pair-level baseline); reservoir subjects are assigned
# unique singleton cluster IDs (singleton frailty degenerates to standard Weibull),
# so the joint likelihood is:
#   L = L_pairs_frailty(beta_T, beta_xs) x L_reservoir_weibull(beta_T, beta_xs).
#
# @details
# As with the IVWC analogue, the matched-pairs frailty estimator differs between the
# univariate and multivariate subclasses:
#
# \strong{Univariate} (\code{include_covariates() == FALSE}): uses \code{parfm::parfm()}
# to fit a fully parametric Weibull gamma-frailty model with only the treatment
# indicator on ALL data (matched + reservoir singletons). The PH coefficient and
# Weibull shape \eqn{\rho} are both taken from \code{parfm}, and the AFT log-time
# ratio is recovered as \eqn{\hat\alpha_\mathrm{AFT} = -\hat\beta_\mathrm{PH} / \hat\rho}.
#
# \strong{Multivariate} (\code{include_covariates() == TRUE}): \code{parfm} becomes
# prohibitively slow with many covariates. Instead, the treatment PH coefficient
# \eqn{\hat\beta_\mathrm{PH}} and its SE are obtained from
# \code{survival::coxph()} with a \code{frailty(cluster, distribution = "gamma")}
# term on ALL data; the Weibull shape \eqn{\hat\rho} is estimated separately from a
# \code{survival::survreg()} fit on all data. The AFT conversion and delta-method SE
# follow identically to the IVWC multivariate subclass.
#
# The \pkg{parfm} package is required only for the univariate subclass.
#
# @keywords internal
DesignInferenceAbstractKKWeibullFrailtyCombinedLikelihood = R6::R6Class("DesignInferenceAbstractKKWeibullFrailtyCombinedLikelihood",
	inherit = DesignInferenceKKPassThrough,
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
			assertResponseType(seq_des_obj$get_response_type(), c("survival", "count"))
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		# @description
		# Returns the combined-likelihood estimate of the treatment effect (AFT log-time ratio).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		# @description
		# Returns a 1 - alpha confidence interval for the treatment effect.
		# @param alpha Significance level; default 0.05 gives a 95% CI.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Returns a 2-sided p-value for H0: beta_T = delta.
		# @param delta Null value; default 0.
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

			# For univariate, use the fast approximation instead of parfm
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			cluster_id = m_vec
			reservoir_idx = which(cluster_id == 0L)
			if (length(reservoir_idx) > 0L){
				cluster_id[reservoir_idx] = max(cluster_id) + seq_along(reservoir_idx)
			}

			dat = data.frame(
				time    = private$y,
				event   = private$dead,
				w       = private$w,
				cluster = factor(cluster_id)
			)

			res = private$fit_fast_approx(dat)
			if (!is.null(res) && is.finite(res$beta_aft)){
				return(res$beta_aft)
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
				# Resample reservoir
				i_res_b = sample(i_reservoir, n_reservoir, replace = TRUE)

				# Resample matched pairs
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
				} else {
					i_match_b = integer(0)
					cluster_b = integer(0)
				}

				i_all = c(i_res_b, i_match_b)
				# Unique cluster IDs for singletons
				cluster_id = c(m + seq_along(i_res_b), cluster_b)

				dat_b = data.frame(
					time = y[i_all],
					event = dead[i_all],
					w = w[i_all],
					cluster = factor(cluster_id)
				)
				if (!is.null(X)) dat_b = cbind(dat_b, X[i_all, , drop = FALSE])

				res = private$fit_fast_approx(dat_b, cov_str)
				if (!is.null(res) && is.finite(res$beta_aft)) res$beta_aft else NA_real_
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
			if (is.null(cox_mod) && !is.null(cov_str)){
				# Fallback: drop covariates
				cox_mod = tryCatch(
					survival::coxph(survival::Surv(time, event) ~ w +
					                survival::frailty(cluster, distribution = "gamma"),
					                data = dat),
					error = function(e) NULL
				)
			}
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
			if (is.null(survreg_mod) && !is.null(cov_str)){
				survreg_mod = tryCatch(survival::survreg(survival::Surv(time, event) ~ w,
				                                         data = dat, dist = "weibull"),
				                       error = function(e) NULL)
			}
			if (is.null(survreg_mod)) return(NULL)

			rho = 1 / survreg_mod$scale
			list(beta_aft = -beta_ph / rho, ssq_aft = (se_ph / rho)^2)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Weibull Frailty combined-likelihood: could not compute a finite standard error.")
			}
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Cluster IDs: matched subjects keep their pair ID;
			# reservoir subjects (m_vec == 0) each get a unique singleton ID.
			cluster_id = m_vec
			reservoir_idx = which(cluster_id == 0L)
			if (length(reservoir_idx) > 0L){
				cluster_id[reservoir_idx] = max(cluster_id) + seq_along(reservoir_idx)
			}

			dat = data.frame(
				time    = private$y,
				event   = private$dead,
				w       = private$w,
				cluster = factor(cluster_id)
			)

			if (private$include_covariates()){
				X_all = as.matrix(private$X)
				colnames(X_all) = paste0("x", seq_len(ncol(X_all)))
				dat = cbind(dat, X_all)
				cov_str = paste(colnames(X_all), collapse = " + ")

				res = private$fit_fast_approx(dat, cov_str)
				if (is.null(res)){
					private$cached_values$beta_hat_T   = NA_real_
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}
				alpha_aft     = res$beta_aft
				ssq_alpha_aft = res$ssq_aft
			} else {
				# Univariate: parfm parametric Weibull gamma-frailty on all data
				# (singleton clusters for reservoir subjects degenerate to standard Weibull).
				mod = tryCatch({
					parfm::parfm(
						survival::Surv(time, event) ~ w,
						cluster = "cluster",
						data    = dat,
						dist    = "weibull",
						frailty = "gamma"
					)
				}, error = function(e) NULL)

				if (is.null(mod)){
					private$cached_values$beta_hat_T   = NA_real_
					private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}

				rho           = as.numeric(mod["rho", "ESTIMATE"])
				beta_ph       = as.numeric(mod["w",   "ESTIMATE"])
				se_ph         = as.numeric(mod["w",   "SE"])
				alpha_aft     = -beta_ph / rho
				ssq_alpha_aft = (se_ph / rho)^2
			}

			private$cached_values$beta_hat_T   = if (is.finite(alpha_aft)) alpha_aft else NA_real_
			private$cached_values$s_beta_hat_T = if (is.finite(ssq_alpha_aft) && ssq_alpha_aft > 0)
			                                        sqrt(ssq_alpha_aft) else NA_real_
			private$cached_values$is_z         = TRUE
		}
	)
)
