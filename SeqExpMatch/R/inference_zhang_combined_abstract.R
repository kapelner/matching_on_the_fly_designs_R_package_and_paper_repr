# Abstract base class for Zhang combined randomisation CIs
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
			m = max(match_indic, na.rm = TRUE)
			y = private$y
			w = private$w

			yTs_matched     = array(NA, m)
			yCs_matched     = array(NA, m)
			y_matched_diffs = array(NA, m)
			X_matched_diffs = matrix(NA, nrow = m, ncol = ncol(private$X))
			if (m > 0){
				match_data      = match_diffs_cpp(w, match_indic, y, private$X, m)
				yTs_matched     = match_data$yTs_matched
				yCs_matched     = match_data$yCs_matched
				X_matched_diffs = match_data$X_matched_diffs
				nonzero_cols    = apply(X_matched_diffs, 2, function(col) any(col != 0))
				X_matched_diffs = X_matched_diffs[, nonzero_cols, drop = FALSE]
				y_matched_diffs = yTs_matched - yCs_matched
			}
			w_reservoir = w[match_indic == 0]

			private$cached_values$KKstats = list(
				X_matched_diffs = X_matched_diffs,
				yTs_matched     = yTs_matched,
				yCs_matched     = yCs_matched,
				y_matched_diffs = y_matched_diffs,
				X_reservoir     = private$X[match_indic == 0, , drop = FALSE],
				y_reservoir     = y[match_indic == 0],
				w_reservoir     = w_reservoir,
				nRT             = sum(w_reservoir, na.rm = TRUE),
				nRC             = sum(w_reservoir == 0, na.rm = TRUE),
				m               = m
			)
		},

		# -----------------------------------------------------------------------
		# Core Zhang Bisection Algorithm
		# -----------------------------------------------------------------------

		ci_rand_zhang_combined = function(alpha, pval_epsilon, combination_method = "Fisher"){
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

			est = self$compute_treatment_estimate()
			if (!is.finite(est)){
				stop("Cannot compute randomisation CI: point estimate is not finite.")
			}

			# Expand the MLE CI (at 2*alpha) by 50% on each side
			mle_ci   = self$compute_mle_confidence_interval(alpha * 2)
			ci_width = mle_ci[2] - mle_ci[1]
			lo_bound = mle_ci[1] - 0.5 * ci_width
			hi_bound = mle_ci[2] + 0.5 * ci_width

			p_fn = function(delta_0){
				p_M = if (m > 0)             private$compute_rand_pval_matched_pairs(delta_0) else NA_real_
				p_R = if (nRT > 0 && nRC > 0) private$compute_rand_pval_reservoir(delta_0)    else NA_real_
				private$combine_rand_pvals(p_M, p_R, m, nRT, nRC, combination_method)
			}

			lower = private$bisect_ci_boundary(p_fn, inside = est, outside = lo_bound, pval_th = alpha, tol = pval_epsilon)
			upper = private$bisect_ci_boundary(p_fn, inside = est, outside = hi_bound, pval_th = alpha, tol = pval_epsilon)

			c(lower, upper)
		},

		# -----------------------------------------------------------------------
		# Abstract hooks for subclasses
		# -----------------------------------------------------------------------

		compute_rand_pval_matched_pairs = function(delta_0) stop("must implement"),
		compute_rand_pval_reservoir     = function(delta_0) stop("must implement"),

		# -----------------------------------------------------------------------
		# Combination Rules
		# -----------------------------------------------------------------------

		combine_rand_pvals = function(p_M, p_R, m, nRT, nRC, method){
			has_M = m > 0              && is.finite(p_M) && p_M > 0
			has_R = nRT > 0 && nRC > 0 && is.finite(p_R) && p_R > 0

			if (has_M && has_R){
				switch(method,
					Fisher = {
						pchisq(-2 * (log(p_M) + log(p_R)), df = 4, lower.tail = FALSE)
					},
					Stouffer = {
						z_M = qnorm(1 - p_M / 2)
						z_R = qnorm(1 - p_R / 2)
						z_combined = (z_M + z_R) / sqrt(2)
						2 * pnorm(-abs(z_combined))
					},
					min_p = {
						1 - (1 - min(p_M, p_R))^2
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
