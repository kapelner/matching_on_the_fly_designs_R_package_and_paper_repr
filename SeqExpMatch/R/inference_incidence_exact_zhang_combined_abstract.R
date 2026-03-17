# Abstract base class for Zhang combined exact CIs
#
# @description
# Consolidates the shared logic for the Zhang (2026) combined test-inversion
# algorithm. Handles match data calculation, the bisection solver, and
# p-value combination rules.
#
# @keywords internal
SeqDesignInferenceAbstractZhangCombinedBase = R6::R6Class("SeqDesignInferenceAbstractZhangCombinedBase",
	inherit = SeqDesignInference,
	public = list(

		# @description
		# Initialize the object.
		# @param seq_des_obj  A SeqDesign object.
		# @param num_cores    Number of CPU cores for parallel processing.
		# @param verbose      Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			# For KK designs set up match data eagerly (mirrors SeqDesignInferenceKKPassThrough)
			if (is(seq_des_obj, "SeqDesignKK14")){
				private$match_indic = seq_des_obj$.__enclos_env__$private$match_indic
				private$compute_basic_match_data()
			}
		}
	),

	private = list(

		match_indic = NULL,

		# -----------------------------------------------------------------------
		# KK match-data setup (mirrors SeqDesignInferenceKKPassThrough)
		# -----------------------------------------------------------------------

		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			match_indic = private$match_indic
			if (is.null(match_indic)){
				match_indic = rep(0, private$n)
			}
			match_indic[is.na(match_indic)] = 0
			y = private$y
			w = private$w

			private$cached_values$KKstats = compute_zhang_match_data_cpp(w, match_indic, y, private$X)
		},

		# -----------------------------------------------------------------------
		# Core Zhang Bisection Algorithm
		# -----------------------------------------------------------------------

		ci_exact_zhang_combined = function(alpha, pval_epsilon, combination_method = "Fisher"){
			# Component sizes
			if (!is.null(private$cached_values$KKstats)){
				m   = private$cached_values$KKstats$m
				nRT = private$cached_values$KKstats$nRT
				nRC = private$cached_values$KKstats$nRC
			} else {
				m   = 0L
				nRT = sum(private$w == 1L, na.rm = TRUE)
				nRC = sum(private$w == 0L, na.rm = TRUE)
			}

			est = private$compute_treatment_estimate_internal()
			if (!is.finite(est)){
				stop("Cannot compute exact CI: point estimate is not finite.")
			}

			# Expand the MLE CI (at 2*alpha) by 50% on each side
			mle_ci   = private$compute_asymp_confidence_interval_internal(alpha * 2)
			ci_width = mle_ci[2] - mle_ci[1]
			lo_bound = mle_ci[1] - 0.5 * ci_width
			hi_bound = mle_ci[2] + 0.5 * ci_width

			p_fn = function(delta_0){
				p_M = if (m > 0)             private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
				p_R = if (nRT > 0 && nRC > 0) private$compute_exact_pval_reservoir(delta_0)    else NA_real_
				private$combine_exact_pvals(p_M, p_R, m, nRT, nRC, combination_method)
			}

			lower = private$bisect_ci_boundary(p_fn, inside = est, outside = lo_bound, pval_th = alpha, tol = pval_epsilon)
			upper = private$bisect_ci_boundary(p_fn, inside = est, outside = hi_bound, pval_th = alpha, tol = pval_epsilon)

			c(lower, upper)
		},

		# -----------------------------------------------------------------------
		# Abstract hooks for subclasses
		# -----------------------------------------------------------------------

		compute_treatment_estimate_internal      = function() stop("must implement"),
		compute_asymp_confidence_interval_internal = function(alpha) stop("must implement"),
		
		# Children must implement these to provide the per-component p-values
		compute_exact_pval_matched_pairs         = function(delta_0) stop("must implement"),
		compute_exact_pval_reservoir             = function(delta_0) stop("must implement"),

		# -----------------------------------------------------------------------
		# Combination Rules
		# -----------------------------------------------------------------------

		combine_exact_pvals = function(p_M, p_R, m, nRT, nRC, method){
			has_M = m > 0              && is.finite(p_M) && p_M > 0
			has_R = nRT > 0 && nRC > 0 && is.finite(p_R) && p_R > 0

			if (has_M && has_R){
				switch(method,
					Fisher = {
						stats::pchisq(-2 * (base::log(p_M) + base::log(p_R)), df = 4, lower.tail = FALSE)
					},
					Stouffer = {
						z_M = stats::qnorm(1 - p_M / 2)
						z_R = stats::qnorm(1 - p_R / 2)
						z_combined = (z_M + z_R) / base::sqrt(2)
						2 * stats::pnorm(-base::abs(z_combined))
					},
					min_p = {
						1 - (1 - base::min(p_M, p_R))^2
					}
				)
			} else if (has_M){
				p_M
			} else if (has_R){
				p_R
			} else {
				NA_real_
			}
		},

		combine_rand_pvals = function(p_M, p_R, m, nRT, nRC, method){
			private$combine_exact_pvals(p_M, p_R, m, nRT, nRC, method)
		},

		# -----------------------------------------------------------------------
		# Bisection Solver
		# -----------------------------------------------------------------------

		bisect_ci_boundary = function(p_fn, inside, outside, pval_th, tol){
			for (iter in seq_len(50L)){
				mid   = (inside + outside) / 2
				p_mid = tryCatch(p_fn(mid), error = function(e) NA_real_)
				if (is.na(p_mid)) p_mid = 0

				if (abs(p_mid - pval_th) < tol || abs(outside - inside) < 1e-8) break

				if (p_mid > pval_th){
					inside  = mid
				} else {
					outside = mid
				}
			}
			(inside + outside) / 2
		}
	)
)
