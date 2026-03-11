# Abstract class for Conditional Logistic Combined-Likelihood Compound Inference
#
# @description
# Fits a single joint likelihood over all KK design data for incidence responses.
# The matched-pair component uses the conditional logistic likelihood on discordant
# pairs (signed within-pair differences); the reservoir component uses the marginal
# logistic likelihood. Both share the treatment effect beta_T and covariate slopes
# beta_xs.
#
# Combined design matrix (column layout: beta_0, beta_T, beta_xs):
#   Discordant pair rows: [0 | t_diff_k | X_diff_k]
#   Reservoir rows:       [1 | w_i      | X_i     ]
# The zero in column 1 for pair rows encodes that beta_0 is absent from the
# conditional likelihood (eliminated by conditioning on the pair sum).
# Fitting a single logistic regression on the stacked dataset maximises the
# combined log-likelihood L = L_cond(pairs) + L_marg(reservoir).
#
# @keywords internal
SeqDesignInferenceAbstractKKClogitCombinedLikelihood = R6::R6Class("SeqDesignInferenceAbstractKKClogitCombinedLikelihood",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
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
				stop("Clogit combined-likelihood: could not compute a finite standard error (possible perfect separation or insufficient data).")
			}
		},

		# Fit the combined logistic likelihood over discordant matched-pair differences
		# and reservoir observations with SHARED covariate effects beta_xs.
		#
		# Column layout of X_comb: [beta_0 | beta_T | beta_xs (p cols)]
		#   col 1       : beta_0   (0 for discordant pair rows; 1 for reservoir rows)
		#   col 2       : beta_T   (t_diff for pair rows; w_i for reservoir rows)
		#   cols 3..p+2 : beta_xs  (X_diff for pair rows; X_i for reservoir rows)
		#
		# Pairs-only case: beta_0 column is all-zero and dropped; layout [beta_T | beta_xs].
		# Reservoir-only case: standard logistic layout [beta_0 | beta_T | beta_xs].
		#
		# Discordant-pair processing and design-matrix construction are done in C++
		# (build_kk_combined_clogit_design_cpp / collect_discordant_pairs_cpp) to
		# avoid O(m^2) R vector reallocation and intermediate matrix allocations.
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

			# ---- Build combined design matrix ------------------------------------
			X_comb   = NULL
			y_comb   = NULL
			j_beta_T = 2L

			if (m > 0){
				match_indic = private$match_indic
				if (is.null(match_indic)) match_indic = rep(0L, private$n)
				match_indic[is.na(match_indic)] = 0L
				i_matched = which(match_indic > 0)
				y_m      = private$y[i_matched]
				w_m      = private$w[i_matched]
				strata_m = match_indic[i_matched]
				X_mat    = if (p > 0L) as.matrix(private$X[i_matched, , drop = FALSE]) else matrix(nrow = length(y_m), ncol = 0L)

				if (has_reservoir){
					# Combined: build [0|t_diff|X_diff; 1|w_r|X_r] in one C++ call.
					# If all pairs are concordant (nd=0) the result contains only
					# reservoir rows — still valid with j_beta_T = 2L.
					y_r    = KKstats$y_reservoir
					w_r    = KKstats$w_reservoir
					X_r    = if (p > 0L) as.matrix(KKstats$X_reservoir) else matrix(nrow = length(y_r), ncol = 0L)
					design = build_kk_combined_clogit_design_cpp(
						as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m),
						as.double(y_r), as.double(w_r), X_r
					)
					X_comb   = design$X_comb
					y_comb   = design$y_comb
					j_beta_T = 2L
				} else {
					# Pairs only: collect discordant pairs; drop beta_0 column.
					res = collect_discordant_pairs_cpp(
						as.double(y_m), as.double(w_m), X_mat, as.integer(strata_m)
					)
					if (res$nd > 0){
						X_comb   = if (p > 0L) cbind(res$t_diffs, res$X_diffs) else matrix(res$t_diffs, ncol = 1L)
						y_comb   = res$y_01
						j_beta_T = 1L
					}
				}
			} else if (has_reservoir){
				# Reservoir only: standard marginal logistic [beta_0 | beta_T | beta_xs]
				y_r    = KKstats$y_reservoir
				w_r    = KKstats$w_reservoir
				X_comb = if (p > 0L) cbind(1, w_r, as.matrix(KKstats$X_reservoir)) else cbind(1, w_r)
				y_comb = y_r
			}

			if (is.null(X_comb)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- QR-reduce to full rank, preserving beta_T column ---------------
			qr_X = qr(X_comb)
			if (qr_X$rank < ncol(X_comb)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_beta_T %in% keep)) keep[qr_X$rank] = j_beta_T
				keep     = sort(keep)
				X_comb   = X_comb[, keep, drop = FALSE]
				j_beta_T = which(keep == j_beta_T)
			}

			n_total  = nrow(X_comb)
			n_params = ncol(X_comb)
			if (n_total <= n_params){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- Fit combined logistic ------------------------------------------
			mod = tryCatch(
				fast_logistic_regression_with_var(X_comb, y_comb, j = j_beta_T),
				error = function(e) NULL
			)
			if (is.null(mod) || !is.finite(mod$b[j_beta_T])){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(mod$b[j_beta_T])
			private$cached_values$s_beta_hat_T = sqrt(mod$ssq_b_j)
			private$cached_values$is_z         = TRUE
			invisible(NULL)
		}
	)
)
