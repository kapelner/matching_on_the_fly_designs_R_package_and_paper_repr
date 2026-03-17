# Abstract class for Stratified Cox / Standard Cox Compound Inference
#
# @description
# This class implements a compound estimator for KK matching-on-the-fly designs with
# survival responses. For matched pairs, it uses stratified Cox proportional hazards
# regression (each pair is a stratum). For reservoir subjects, it uses standard Cox
# regression. The two estimates (both log-hazard ratios) are combined via a
# variance-weighted linear combination.
#
# @keywords internal
SeqDesignInferenceAbstractKKStratCoxIVWC = R6::R6Class("SeqDesignInferenceAbstractKKStratCoxIVWC",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "survival")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
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

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

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

			# --- Matched pairs: Stratified Cox ---
			if (m > 0){
				private$strat_cox_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: Standard Cox ---
			if (nRT > 0 && nRC > 0){
				private$cox_for_reservoir()
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
				stop("Stratified Cox compound estimator: could not compute a finite standard error.")
			}
		},

		strat_cox_for_matched_pairs = function(){
			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			i_matched = which(match_indic > 0)
			y_m       = private$y[i_matched]
			dead_m    = private$dead[i_matched]
			w_m       = private$w[i_matched]
			strata_m  = match_indic[i_matched]

			# Filter strata that have no events (provide no information for Cox partial likelihood)
			strata_with_events = unique(strata_m[dead_m == 1])
			if (length(strata_with_events) == 0) return(invisible(NULL))

			i_valid = which(strata_m %in% strata_with_events)

			dat = data.frame(y = y_m[i_valid], dead = dead_m[i_valid], w = w_m[i_valid], strata = strata_m[i_valid])
			formula_str = "survival::Surv(y, dead) ~ w + strata(strata)"

			if (private$include_covariates()){
				X_m = as.matrix(private$X[i_matched[i_valid], , drop = FALSE])
				# QR-reduce to full rank while always preserving the treatment column
				X_full = cbind(w = dat$w, X_m)
				qr_full = qr(X_full)
				r_full  = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(1L %in% keep)) keep[r_full] = 1L
					keep    = sort(keep)
					X_full  = X_full[, keep, drop = FALSE]
				}
				# Remove 'w' from X_full as it's already in the formula
				X_covs = X_full[, colnames(X_full) != "w", drop = FALSE]
				if (ncol(X_covs) > 0){
					colnames(X_covs) = paste0("x", 1:ncol(X_covs))
					dat = cbind(dat, X_covs)
					formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
				}
			}

			mod = tryCatch(
				survival::coxph(as.formula(formula_str), data = dat),
				error = function(e) NULL
			)
			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(coef(mod)["w"])
			se   = sqrt(as.numeric(vcov(mod)["w", "w"]))
			private$cached_values$beta_T_matched     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(se) && se > 0) se^2 else NA_real_
		},

		cox_for_reservoir = function(){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			match_indic_safe = private$match_indic
			if (is.null(match_indic_safe)) match_indic_safe = rep(0L, private$n)
			match_indic_safe[is.na(match_indic_safe)] = 0L
			dead_r = private$dead[match_indic_safe == 0]
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)

			dat = data.frame(y = y_r, dead = dead_r, w = w_r)
			formula_str = "survival::Surv(y, dead) ~ w"

			if (private$include_covariates()){
				X_full = cbind(w = dat$w, X_r)
				qr_full = qr(X_full)
				r_full  = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(1L %in% keep)) keep[r_full] = 1L
					keep    = sort(keep)
					X_full  = X_full[, keep, drop = FALSE]
				}
				X_covs = X_full[, colnames(X_full) != "w", drop = FALSE]
				if (ncol(X_covs) > 0){
					colnames(X_covs) = paste0("x", 1:ncol(X_covs))
					dat = cbind(dat, X_covs)
					formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
				}
			}

			mod = tryCatch(
				survival::coxph(as.formula(formula_str), data = dat),
				error = function(e) NULL
			)

			# Fallback if multivariate failed
			if (is.null(mod) && private$include_covariates()){
				mod = tryCatch(
					survival::coxph(survival::Surv(y, dead) ~ w, data = data.frame(y = y_r, dead = dead_r, w = w_r)),
					error = function(e) NULL
				)
			}

			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(coef(mod)["w"])
			se   = sqrt(as.numeric(vcov(mod)["w", "w"]))
			private$cached_values$beta_T_reservoir     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(se) && se > 0) se^2 else NA_real_
		}
	)
)
