#' Pattern-1 mixin: KK compound pass-through behavior
#'
#' Provides the five private methods shared by all KK compound estimators:
#' \code{supports_likelihood_tests}, \code{compute_basic_match_data},
#' \code{compute_estimate_from_matched_and_reservoir}, \code{only_matches},
#' \code{only_reservoir}, and \code{compute_reservoir_and_match_statistics}.
#'
#' Splice into a daughter class via
#' \code{private = c(InferenceMixinKKPassThroughCompound$private, list(...))}.
#' The capability flag \code{private$kk_passthrough_compound} is set to \code{TRUE}.
#'
#' @keywords internal
#' @noRd
InferenceMixinKKPassThroughCompound = list(
	public = list(),
	private = list(
		kk_passthrough_compound = TRUE,
		is_a_kk_compound_estimator = function() TRUE,

		supports_likelihood_tests = function() FALSE,

		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),

		# Caches only the column-keep decision from the QR rank check (keyed by
		# the input column count) and re-applies it to whatever X/j_treat is
		# passed on each call, rather than returning a frozen matrix.
		#
		# This matters because the combined-design estimators
		# (OLS/robust-regression, one-lik/IVWC) call this from fit_combined()
		# with a design matrix whose column count depends on which branch was
		# taken (matched-pairs-only vs reservoir-only vs the full combined
		# design), which in turn depends on the current nRT/nRC split. During
		# the fast per-permutation randomization path (compute_rand_two_sided_pval),
		# a single reused worker calls this many times with a fresh X/j_treat
		# each iteration while its cache is deliberately preserved across
		# iterations (m_mat is NULL for KK designs, so the reduction is assumed
		# stable) -- but the branch can still flip whenever a reservoir split
		# becomes one-arm-only (nRT == 0 or nRC == 0). Returning a frozen cached
		# matrix in that case reuses stale, wrong-shaped data against the
		# current iteration's y: at best a degenerate randomization
		# distribution, at worst a length mismatch that corrupts memory in the
		# C++ IRLS/OLS kernels since Eigen's bounds assertions are compiled out
		# (NDEBUG is force-defined ahead of every Eigen include). Keying the
		# cache on ncol(X) forces a fresh recomputation on a branch change
		# instead of reusing a stale reduction.
		reduce_design_matrix_once = function(X, j_treat, cache_key){
			cached = private$cached_values[[cache_key]]
			if (!is.null(cached) && identical(cached$ncol_in, ncol(X))){
				keep = cached$keep
				if (is.null(keep)) return(list(X = X, j_treat = j_treat))
				return(list(X = X[, keep, drop = FALSE], j_treat = which(keep == j_treat)))
			}
			qr_X = qr(X)
			keep = NULL
			if (qr_X$rank < ncol(X)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_treat %in% keep)) keep[qr_X$rank] = j_treat
				keep = sort(unique(keep))
			}
			private$cached_values[[cache_key]] = list(ncol_in = ncol(X), keep = keep)
			if (is.null(keep)) list(X = X, j_treat = j_treat) else list(X = X[, keep, drop = FALSE], j_treat = which(keep == j_treat))
		},

		compute_estimate_from_matched_and_reservoir = function(run_matched, run_reservoir){
			if (!isTRUE(private$has_match_structure)) {
				private$cache_nonestimable_estimate("kk_design_required")
				return(invisible(NULL))
			}
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			if (is.null(private$cached_values$KKstats)) {
				return(invisible(NULL))
			}
			if (private$only_matches()){
				run_matched()
			} else if (private$only_reservoir()){
				run_reservoir()
			} else {
				run_matched()
				run_reservoir()
			}
		},

		only_matches = function(){
			if (is.null(private$cached_values$KKstats)) return(FALSE)
			nRT = private$cached_values$KKstats$nRT
			nRC = private$cached_values$KKstats$nRC
			if (!is.finite(nRT) || !is.finite(nRC)) return(FALSE)
			nRT <= 1 || nRC <= 1
		},

		only_reservoir = function(){
			if (is.null(private$cached_values$KKstats)) return(FALSE)
			m = private$cached_values$KKstats$m
			is.finite(m) && m <= 1
		},

		compute_reservoir_and_match_statistics = function(){
			nRC = private$cached_values$KKstats$nRC
			nRT = private$cached_values$KKstats$nRT
			nR = nRT + nRC
			m = private$cached_values$KKstats$m

			y_reservoir_T = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 1]
			y_reservoir_C = private$cached_values$KKstats$y_reservoir[private$cached_values$KKstats$w_reservoir == 0]

			ssqD_bar = if (is.finite(m) && m > 1){
							var_cpp(private$cached_values$KKstats$y_matched_diffs) / m
						} else {
							NA_real_
						}
			ssqR = if (is.finite(nRT) && is.finite(nRC) && nRT > 1 && nRC > 1 && nR > 2){
						(var_cpp(y_reservoir_T) * (nRT - 1) + var_cpp(y_reservoir_C) * (nRC - 1)) /
							(nR - 2) * (1 / nRT + 1 / nRC)
					} else {
						NA_real_
					}

			private$cached_values$KKstats$d_bar = if (is.finite(m) && m > 0) mean_cpp(private$cached_values$KKstats$y_matched_diffs) else {
				NA_real_
			}
			private$cached_values$KKstats$ssqD_bar = ssqD_bar
			private$cached_values$KKstats$r_bar = if (is.finite(nRT) && is.finite(nRC) && nRT > 0 && nRC > 0) mean_cpp(y_reservoir_T) - mean_cpp(y_reservoir_C) else {
				NA_real_
			}
			private$cached_values$KKstats$ssqR = ssqR
			private$cached_values$KKstats$w_star = if (is.finite(ssqR) && is.finite(ssqD_bar) && (ssqR + ssqD_bar) > 0) {
				ssqR / (ssqR + ssqD_bar)
			} else {
				NA_real_
			}
		}
	)
)
