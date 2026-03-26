# Abstract class for Stratified Cox Combined-Likelihood Compound Inference
#
# @description
# Fits a single joint Cox model over all KK design data for survival responses.
# Reservoir subjects (m_vec == 0) share stratum 0 (common baseline hazard);
# matched subjects belong to their own pair stratum (strata 1..m), which conditions
# away each pair's baseline hazard identically to stratified Cox.
#
# The combined partial likelihood is therefore:
#   L = L_strat_pairs(beta_T, beta_xs) x L_cox_reservoir(beta_T, beta_xs)
# and is maximised by a single coxph call with strata(stratum).
#
# @keywords internal
DesignInferenceAbstractKKStratCoxCombinedLikelihood = R6::R6Class("DesignInferenceAbstractKKStratCoxCombinedLikelihood",
	inherit = DesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "survival")
			if (!is(des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
		},

		# @description
		# Returns the combined-likelihood estimate of the treatment effect (log-HR).
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

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Stratified Cox combined-likelihood: could not compute a finite standard error.")
			}
		},

		shared_combined_likelihood = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Reservoir subjects get stratum 0 (shared baseline);
			# matched subjects get their pair ID as stratum.
			dat = data.frame(
				y       = private$y,
				dead    = private$dead,
				w       = private$w,
				stratum = m_vec
			)

			formula_str = "survival::Surv(y, dead) ~ w + strata(stratum)"

			if (private$include_covariates()){
				X_full = cbind(w = private$w, as.matrix(private$X))
				qr_full = qr(X_full)
				r_full  = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(1L %in% keep)) keep[r_full] = 1L
					keep   = sort(keep)
					X_full = X_full[, keep, drop = FALSE]
				}
				X_covs = X_full[, colnames(X_full) != "w", drop = FALSE]
				if (ncol(X_covs) > 0){
					colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
					dat = cbind(dat, X_covs)
					formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
				}
			}

			mod = tryCatch(
				survival::coxph(as.formula(formula_str), data = dat),
				error = function(e) NULL
			)

			# Fallback: drop covariates if multivariate fit failed
			if (is.null(mod) && private$include_covariates()){
				mod = tryCatch(
					survival::coxph(survival::Surv(y, dead) ~ w + strata(stratum), data = dat),
					error = function(e) NULL
				)
			}

			if (is.null(mod) || !is.finite(coef(mod)["w"])){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(coef(mod)["w"])
			se = sqrt(as.numeric(vcov(mod)["w", "w"]))
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z         = TRUE
			invisible(NULL)
		}
	)
)
