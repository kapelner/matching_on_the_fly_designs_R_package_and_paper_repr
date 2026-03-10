#' Abstract mixin: randomisation-based CI for binary (incidence) outcomes
#'
#' @description
#' Sits between \code{SeqDesignInference} and incidence-focused inference
#' classes.  The single public method it contributes,
#' \code{compute_confidence_interval_rand()}, intercepts calls where the
#' response type is \code{"incidence"} and dispatches to the chosen CI
#' algorithm; for all other response types it delegates straight to
#' \code{super$compute_confidence_interval_rand()} (i.e. the standard
#' continuous / count / proportion / survival implementation in
#' \code{SeqDesignInference}).
#'
#' Currently supported designs for the incidence path:
#' \itemize{
#'   \item \code{SeqDesignCRD} — all subjects treated as a Bernoulli reservoir
#'         (\eqn{m = 0}); \pkg{bclogit} not required.
#'   \item \code{SeqDesignKK14} and subclasses — matched pairs plus a Bernoulli
#'         reservoir; \pkg{bclogit} required when \eqn{m > 0}.
#' }
#' Other designs raise an informative error.  Support for additional designs
#' can be added in future \code{method} implementations without touching this
#' class.
#'
#' \strong{Extension points for subclasses} (used by \code{"zhang_combined"}):
#' \itemize{
#'   \item \code{compute_rand_pval_matched_pairs(delta_0)}: two-sided p-value
#'         for \eqn{H_0 : \mathrm{ATE}_\mathrm{pairs} = \delta_0}.  Return
#'         \code{NA_real_} when there are no matched pairs or no discordant
#'         pairs.
#'   \item \code{compute_rand_pval_reservoir(delta_0)}: two-sided p-value for
#'         \eqn{H_0 : \mathrm{ATE}_\mathrm{reservoir} = \delta_0}.  Return
#'         \code{NA_real_} when either reservoir arm is empty.
#'   \item \code{combine_rand_pvals_fisher / _stouffer / _minp}: combination
#'         rules selected via the \code{combination_method} argument of
#'         \code{compute_confidence_interval_rand()}.
#' }
#'
#' @keywords internal
SeqDesignInferenceIncidExactAbstract = R6::R6Class("SeqDesignInferenceIncidExactAbstract",
	inherit = SeqDesignInference,
	public = list(

		#' @description
		#' Initialize the object.
		#' @param seq_des_obj  A SeqDesign object.
		#' @param num_cores    Number of CPU cores for parallel processing.
		#' @param verbose      Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			# For KK designs set up match data eagerly (mirrors SeqDesignInferenceKKPassThrough)
			if (is(seq_des_obj, "SeqDesignKK14")){
				private$match_indic = seq_des_obj$.__enclos_env__$private$match_indic
				private$compute_basic_match_data()
			}
		},

		#' @description
		#' Compute the randomisation-based CI.
		#'
		#' For \code{response_type != "incidence"} the call is forwarded
		#' unchanged to \code{SeqDesignInference$compute_confidence_interval_rand()}.
		#' For \code{response_type == "incidence"} the chosen \code{method} is
		#' executed.
		#'
		#' @param alpha            Significance level; CI covers \eqn{1-\alpha}.
		#'                         Default \code{0.05}.
		#' @param nsim_exact_test  Passed to Monte Carlo methods; ignored by exact
		#'                         methods.  Default \code{501}.
		#' @param pval_epsilon     Bisection convergence tolerance on the p-value
		#'                         scale.  Default \code{0.005}.
		#' @param show_progress    Forwarded to the base-class method for
		#'                         non-incidence types; ignored by incidence methods.
		#' @param method             CI algorithm to use for incidence outcomes.
		#'   Currently supported:
		#'   \describe{
		#'     \item{\code{"zhang_combined"}}{(default) Exact test-inversion via
		#'       combination of matched-pair and reservoir p-values (Zhang 2026).
		#'       Requires subclasses to implement
		#'       \code{compute_rand_pval_matched_pairs()} and
		#'       \code{compute_rand_pval_reservoir()}.}
		#'   }
		#' @param combination_method  How to combine the matched-pair and reservoir
		#'   p-values.  One of \code{"Fisher"} (default), \code{"Stouffer"}, or
		#'   \code{"min_p"}.
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.005, show_progress = TRUE, method = "zhang_combined", combination_method = "Fisher"){
			if (private$seq_des_obj_priv_int$response_type != "incidence"){
				return(super$compute_confidence_interval_rand(alpha, nsim_exact_test, pval_epsilon, show_progress))
			}
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			method             = match.arg(method, choices = c("zhang_combined"))
			combination_method = match.arg(combination_method, c("Fisher", "Stouffer", "min_p"))

			switch(method,
				zhang_combined = private$ci_rand_combine_matched_pair_and_reservoir_p_values(alpha, pval_epsilon, combination_method)
			)
		}
	),

	private = list(

		match_indic = NULL,

		# -----------------------------------------------------------------------
		# KK match-data setup (mirrors SeqDesignInferenceKKPassThrough so that
		# KK-based subclasses work without going through that branch).
		# For CRD designs this is never called; ci_rand_zhang handles m = 0 via
		# private$w directly.
		# -----------------------------------------------------------------------

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

		# -----------------------------------------------------------------------
		# "zhang_combined" method
		#
		# Test-inversion via combination of independent matched-pair and
		# reservoir p-values, with bisection to find CI boundaries.
		# The combination rule is chosen by private$combination_method.
		# -----------------------------------------------------------------------

		ci_rand_combine_matched_pair_and_reservoir_p_values = function(alpha, pval_epsilon, combination_method){
			est = self$compute_treatment_estimate()
			if (!is.finite(est)){
				stop("Cannot compute randomisation CI: point estimate is not finite.")
			}

			# Expand the MLE CI (at 2*alpha) by 50% on each side for reliable bisection bounds
			mle_ci   = self$compute_mle_confidence_interval(alpha * 2)
			ci_width = mle_ci[2] - mle_ci[1]
			lo_bound = mle_ci[1] - 0.5 * ci_width
			hi_bound = mle_ci[2] + 0.5 * ci_width

			# Combined two-sided p-value: large near est, declines as delta_0 moves away.
			# Switch on combination_method to select Fisher / Stouffer / min_p.
			p_fn = function(delta_0) private$compute_combined_exact_pval(delta_0, combination_method)

			# Bisect to find the two CI boundaries.
			# est is always inside the CI; the expanded MLE bounds are outside.
			lower = private$bisect_ci_boundary(p_fn, inside = est, outside = lo_bound, pval_th = alpha, tol = pval_epsilon)
			upper = private$bisect_ci_boundary(p_fn, inside = est, outside = hi_bound, pval_th = alpha, tol = pval_epsilon)

			c(lower, upper)
		},

		compute_combined_exact_pval = function(delta_0, combination_method = "Fisher"){
			# Design validation: only CRD and KK supported by this method
			is_crd = is(private$seq_des_obj, "SeqDesignCRD")
			is_kk  = is(private$seq_des_obj, "SeqDesignKK14")
			if (!is_crd && !is_kk){
				stop("zhang_combined CI rand is only supported for CRD (SeqDesignCRD) and KK (SeqDesignKK14 or subclass) designs.")
			}

			if (!is.null(private$cached_values$KKstats)){
				m   = private$cached_values$KKstats$m
				nRT = private$cached_values$KKstats$nRT
				nRC = private$cached_values$KKstats$nRC
			} else {
				m   = 0L
				nRT = sum(private$w == 1L, na.rm = TRUE)
				nRC = sum(private$w == 0L, na.rm = TRUE)
			}

			p_M = if (m > 0)              private$compute_rand_pval_matched_pairs(delta_0) else NA_real_
			p_R = if (nRT > 0 && nRC > 0) private$compute_rand_pval_reservoir(delta_0)    else NA_real_

			switch(combination_method,
				Fisher   = private$combine_rand_pvals_fisher(p_M, p_R, m, nRT, nRC),
				Stouffer = private$combine_rand_pvals_stouffer(p_M, p_R, m, nRT, nRC),
				min_p    = private$combine_rand_pvals_minp(p_M, p_R, m, nRT, nRC)
			)
		},

		# Default: no matched-pair test (CRD designs have m = 0 always;
		# KK-specific exact classes override this method when needed).
		compute_rand_pval_matched_pairs = function(delta_0) NA_real_,

		# Exact Fisher test under H0: OR_reservoir = exp(delta_0).
		# For CRD designs (KKstats = NULL) all subjects are treated as the reservoir.
		compute_rand_pval_reservoir = function(delta_0){
			if (!is.null(private$cached_values$KKstats)){
				y_r = private$cached_values$KKstats$y_reservoir
				w_r = private$cached_values$KKstats$w_reservoir
				nRT = private$cached_values$KKstats$nRT
				nRC = private$cached_values$KKstats$nRC
			} else {
				y_r = private$y
				w_r = private$w
				nRT = sum(w_r == 1L, na.rm = TRUE)
				nRC = sum(w_r == 0L, na.rm = TRUE)
			}
			if (nRT == 0L || nRC == 0L) return(NA_real_)

			n11 = sum(y_r[w_r == 1L])              # successes in treatment
			n01 = sum(y_r[w_r == 0L])              # successes in control
			n10 = as.integer(nRT) - n11            # failures  in treatment
			n00 = as.integer(nRC) - n01            # failures  in control

			if (n11 + n01 == 0L || n10 + n00 == 0L) return(NA_real_)

			# 2x2 table: rows = Y (1/0), cols = arm (T/C)
			# OR = (n11 * n00) / (n10 * n01) — standard treatment-vs-control OR
			fisher.test(
				matrix(c(n11, n10, n01, n00), 2L, 2L),
				or          = exp(delta_0),
				alternative = "two.sided"
			)$p.value
		},

		# -----------------------------------------------------------------------
		# Combination rules
		#
		# Fisher:   -2(log p_M + log p_R) ~ chi^2(4)
		# Stouffer: signed-z combination (equal weights)
		# min_p:    Tippett's method — P(min <= t) = 1-(1-t)^2
		#
		# When only one component is available, all three return the raw p-value.
		# -----------------------------------------------------------------------

		combine_rand_pvals_fisher = function(p_M, p_R, m, nRT, nRC){
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

		combine_rand_pvals_stouffer = function(p_M, p_R, m, nRT, nRC){
			has_M = m > 0              && is.finite(p_M) && p_M > 0
			has_R = nRT > 0 && nRC > 0 && is.finite(p_R) && p_R > 0

			if (has_M && has_R){
				# Convert two-sided p-values to unsigned z-scores, combine with equal weights
				z_M = qnorm(1 - p_M / 2)
				z_R = qnorm(1 - p_R / 2)
				z_combined = (z_M + z_R) / sqrt(2)
				2 * pnorm(-abs(z_combined))
			} else if (has_M){
				p_M
			} else if (has_R){
				p_R
			} else {
				NA_real_
			}
		},

		combine_rand_pvals_minp = function(p_M, p_R, m, nRT, nRC){
			has_M = m > 0              && is.finite(p_M) && p_M > 0
			has_R = nRT > 0 && nRC > 0 && is.finite(p_R) && p_R > 0

			if (has_M && has_R){
				# Tippett's method: P(min(p1,p2) <= t) = 1 - (1-t)^2 under independence
				1 - (1 - min(p_M, p_R))^2
			} else if (has_M){
				p_M
			} else if (has_R){
				p_R
			} else {
				NA_real_
			}
		},

		# -----------------------------------------------------------------------
		# Bisection helper
		#
		# Locates the delta_0 where p_fn(delta_0) == pval_th.
		# Preconditions: p_fn(inside) > pval_th, p_fn(outside) < pval_th.
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
