#' Abstract mixin: Zhang combined randomisation CI for quantile regression
#'
#' @description
#' Analogous to \code{SeqDesignInferenceIncidExactAbstract} but for quantile
#' regression estimators.  Provides \code{compute_confidence_interval_rand()}
#' via Zhang's combined test-inversion method for both CRD (\eqn{m = 0},
#' all subjects in the reservoir) and KK matching-on-the-fly designs
#' (\eqn{m > 0}).
#'
#' \strong{Extension points for subclasses:}
#' \itemize{
#'   \item \code{compute_rand_pval_matched_pairs(delta_0)}: two-sided p-value
#'         for \eqn{H_0 : Q_\tau(\text{pairs component}) = \delta_0}.
#'         Return \code{NA_real_} when there are no matched pairs.
#'   \item \code{compute_rand_pval_reservoir(delta_0)}: two-sided p-value for
#'         \eqn{H_0 : Q_\tau(\text{reservoir component}) = \delta_0}.
#'         Return \code{NA_real_} when either reservoir arm is empty.
#' }
#'
#' @keywords internal
SeqDesignInferenceAbstractQuantileRandCI = R6::R6Class("SeqDesignInferenceAbstractQuantileRandCI",
	inherit = SeqDesignInference,
	public = list(

		#' @param seq_des_obj		A SeqDesign object.
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			# For KK designs, set up match data eagerly (mirrors SeqDesignInferenceKKPassThrough)
			if (is(seq_des_obj, "SeqDesignKK14")){
				private$match_indic = seq_des_obj$.__enclos_env__$private$match_indic
				private$compute_basic_match_data()
			}
		},

		#' @description
		#' Computes a randomization-based confidence interval via Zhang's combined test.
		#' For matched pairs the sign-flip randomization test is used; for reservoir subjects
		#' the permutation test is used. The two component p-values are combined with
		#' Fisher's chi-squared statistic and the CI boundary is found by bisection.
		#'
		#' @param alpha					The confidence level is 1 - \code{alpha}. Default \code{0.05}.
		#' @param nsim_exact_test		Number of random sign-flips / permutations per boundary evaluation. Default \code{499}.
		#' @param pval_epsilon			Bisection convergence tolerance on the p-value scale. Default \code{0.005}.
		#' @param show_progress			Ignored (kept for interface compatibility).
		#'
		#' @return 	A length-2 numeric vector giving the lower and upper CI boundary.
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 499, pval_epsilon = 0.005, show_progress = TRUE){
			if (!is.null(private[["custom_randomization_statistic_function"]])){
				stop("Custom randomization statistic functions are not supported for the Zhang combined CI method used by ", class(self)[1], ". The method uses its own fixed QR-based test statistics.")
			}
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			assertCount(nsim_exact_test, positive = TRUE)
			private$nsim_rand = as.integer(nsim_exact_test)
			private$ci_rand_zhang_combined(alpha, pval_epsilon)
		}
	),

	private = list(

		nsim_rand   = 499L,
		match_indic = NULL,

		# KK match-data setup (mirrors SeqDesignInferenceKKPassThrough so that
		# KK-based subclasses work without going through that branch).
		# For CRD designs this is never called; ci_rand_zhang handles m = 0 via
		# private$w directly.
		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			private$cached_values$KKstats = .compute_kk_basic_match_data(
				X = private$X,
				n = private$n,
				y = private$y,
				w = private$w,
				match_indic = private$match_indic
			)
		},

		ci_rand_zhang_combined = function(alpha, pval_epsilon){
			# Derive component sizes from KKstats (KK) or treatment vector (CRD)
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

			mle_ci   = self$compute_mle_confidence_interval(alpha * 2)
			ci_width = mle_ci[2] - mle_ci[1]
			lo_bound = mle_ci[1] - 0.5 * ci_width
			hi_bound = mle_ci[2] + 0.5 * ci_width

			p_fn = function(delta_0){
				p_M = if (m > 0)             private$compute_rand_pval_matched_pairs(delta_0) else NA_real_
				p_R = if (nRT > 0 && nRC > 0) private$compute_rand_pval_reservoir(delta_0)    else NA_real_
				private$combine_rand_pvals(p_M, p_R, m, nRT, nRC)
			}

			lower = private$bisect_ci_boundary(p_fn, inside = est, outside = lo_bound, pval_th = alpha, tol = pval_epsilon)
			upper = private$bisect_ci_boundary(p_fn, inside = est, outside = hi_bound, pval_th = alpha, tol = pval_epsilon)
			c(lower, upper)
		},

		# Abstract stubs — subclasses MUST override these
		compute_rand_pval_matched_pairs = function(delta_0){
			stop(class(self)[1], " must implement compute_rand_pval_matched_pairs(delta_0)")
		},

		compute_rand_pval_reservoir = function(delta_0){
			stop(class(self)[1], " must implement compute_rand_pval_reservoir(delta_0)")
		},

		# Fisher's chi-squared combination of independent p-values
		combine_rand_pvals = function(p_M, p_R, m, nRT, nRC){
			has_M = m > 0              && is.finite(p_M) && p_M > 0
			has_R = nRT > 0 && nRC > 0 && is.finite(p_R) && p_R > 0
			if (has_M && has_R){
				pchisq(-2 * (log(p_M) + log(p_R)), df = 4, lower.tail = FALSE)
			} else if (has_M){
				p_M
			} else if (has_R){
				p_R
			} else {
				NA_real_
			}
		},

		# Bisection to find CI boundary where p_fn(delta_0) == pval_th.
		# Precondition: p_fn(inside) > pval_th, p_fn(outside) < pval_th.
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
