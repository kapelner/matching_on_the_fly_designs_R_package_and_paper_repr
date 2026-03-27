# Abstract class for LWA-style Marginal Cox / Standard Cox Compound Inference
#
# @description
# This class implements a compound estimator for KK matching-on-the-fly designs with
# survival responses. For matched pairs, it uses a marginal Cox proportional hazards
# model with Lee-Wei-Amato style cluster-robust variance, treating each pair as a
# cluster of size two. For reservoir subjects, it uses standard Cox regression. The
# two estimates (both log-hazard ratios) are combined via a variance-weighted linear
# combination.
#
# @keywords internal
InferenceAbstractKKLWACoxIVWC = R6::R6Class("InferenceAbstractKKLWACoxIVWC",
	inherit = InferenceAsymp,
	public = list(

		# @description
		# Initialize the inference object.
		# @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "survival")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
		},

		# @description
		# Returns the estimated treatment effect (log-hazard ratio).
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

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			if (m > 0){
				private$lwa_cox_for_matched_pairs()
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m  = private$cached_values$ssq_beta_T_matched
			m_ok   = !is.null(beta_m) && is.finite(beta_m) &&
			         !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			if (nRT > 0 && nRC > 0){
				private$cox_for_reservoir()
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r  = private$cached_values$ssq_beta_T_reservoir
			r_ok   = !is.null(beta_r) && is.finite(beta_r) &&
			         !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0

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
				stop("LWA marginal Cox compound estimator: could not compute a finite standard error.")
			}
		},

		build_design_matrix = function(w, X){
			X_full = matrix(w, ncol = 1)
			colnames(X_full) = "w"
			if (private$include_covariates()){
				X_full = cbind(X_full, as.matrix(X))
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(1L %in% keep)) keep[r_full] = 1L
					keep = sort(keep)
					X_full = X_full[, keep, drop = FALSE]
				}
			}
			X_full
		},

		fit_cox_model = function(y, dead, w, X, cluster = NULL, robust = FALSE){
			if (length(y) == 0L || sum(dead) == 0L) return(NULL)

			X_full = private$build_design_matrix(w, X)
			dat = data.frame(y = y, dead = dead, w = X_full[, "w"])
			formula_str = "survival::Surv(y, dead) ~ w"

			X_covs = X_full[, colnames(X_full) != "w", drop = FALSE]
			if (ncol(X_covs) > 0){
				colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
				dat = cbind(dat, X_covs)
				formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
			}

			if (!is.null(cluster)){
				dat$cluster = factor(cluster)
				formula_str = paste(formula_str, "+ cluster(cluster)")
			}

			mod = tryCatch(
				suppressWarnings(survival::coxph(as.formula(formula_str), data = dat, robust = robust)),
				error = function(e) NULL
			)

			if (is.null(mod) && ncol(X_covs) > 0){
				fallback_formula = "survival::Surv(y, dead) ~ w"
				if (!is.null(cluster)){
					fallback_formula = paste(fallback_formula, "+ cluster(cluster)")
				}
				mod = tryCatch(
					suppressWarnings(survival::coxph(as.formula(fallback_formula), data = dat, robust = robust)),
					error = function(e) NULL
				)
			}

			if (is.null(mod)) return(NULL)

			coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
			vcv   = tryCatch(stats::vcov(mod), error = function(e) NULL)
			if (is.null(coefs) || is.null(vcv) || !("w" %in% names(coefs)) || !("w" %in% rownames(vcv))){
				return(NULL)
			}

			beta = as.numeric(coefs["w"])
			se   = sqrt(as.numeric(vcv["w", "w"]))
			if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)

			list(beta = beta, ssq = se^2)
		},

		lwa_cox_for_matched_pairs = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0)
			if (length(i_matched) == 0L) return(invisible(NULL))

			fit = private$fit_cox_model(
				y = private$y[i_matched],
				dead = private$dead[i_matched],
				w = private$w[i_matched],
				X = private$X[i_matched, , drop = FALSE],
				cluster = m_vec[i_matched],
				robust = TRUE
			)
			if (is.null(fit)) return(invisible(NULL))

			private$cached_values$beta_T_matched = fit$beta
			private$cached_values$ssq_beta_T_matched = fit$ssq
		},

		cox_for_reservoir = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_reservoir = which(m_vec == 0L)
			if (length(i_reservoir) == 0L) return(invisible(NULL))

			fit = private$fit_cox_model(
				y = private$y[i_reservoir],
				dead = private$dead[i_reservoir],
				w = private$w[i_reservoir],
				X = private$X[i_reservoir, , drop = FALSE],
				robust = FALSE
			)
			if (is.null(fit)) return(invisible(NULL))

			private$cached_values$beta_T_reservoir = fit$beta
			private$cached_values$ssq_beta_T_reservoir = fit$ssq
		}
	)
)
