#' Zhang Exact Inference Helpers
#'
#' @description Internal method.
#' Standalone functions to support Zhang (2026) exact test-inversion inference.
#' These functions handle the bisection solver, p-value combination rules,
#' and component-wise p-value logic.
#'
#' @name zhang_combine_exact_pvals
#' @param p_M Matched p-value.
#' @param p_R Reservoir p-value.
#' @param m Number of matches.
#' @param nRT Number of treated in reservoir.
#' @param nRC Number of control in reservoir.
#' @param method Combination method (Fisher or Stouffer).
#' @keywords internal
zhang_combine_exact_pvals = function(p_M, p_R, m, nRT, nRC, method){
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
}

zhang_bisect_ci_boundary = function(p_fn, inside, outside, pval_th, tol){
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

# Computes the point estimate used for bisection starting points
zhang_incid_treatment_estimate = function(stats){
	n11 = stats$n11
	n10 = stats$n10
	n01 = stats$n01
	n00 = stats$n00
	# Sample Odds Ratio (with continuity correction)
	log((n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)))
}

# Computes the MLE-based CI used for bisection search bounds
zhang_incid_mle_ci = function(stats, alpha){
	n11 = stats$n11
	n10 = stats$n10
	n01 = stats$n01
	n00 = stats$n00
	est = log((n11 + 0.5) * (n00 + 0.5) / ((n10 + 0.5) * (n01 + 0.5)))
	se  = sqrt(1/(n11+0.5) + 1/(n10+0.5) + 1/(n01+0.5) + 1/(n00+0.5))
	z = stats::qnorm(1 - alpha/2)
	c(est - z * se, est + z * se)
}
