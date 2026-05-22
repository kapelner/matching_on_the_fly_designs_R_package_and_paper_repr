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

zhang_normalize_exact_inference_args = function(type, args_for_type = NULL, pval_epsilon = NULL){
	if (should_run_asserts()) {
		assertChoice(type, c("Zhang"))
		assertList(args_for_type, null.ok = TRUE)
	}
	default_args = list(combination_method = "Fisher")
	if (!is.null(pval_epsilon)) default_args$pval_epsilon = pval_epsilon
	utils::modifyList(setNames(list(default_args), type), if (is.null(args_for_type)) list() else args_for_type)
}

zhang_assert_exact_inference_params = function(inf_obj, type, args_for_type){
	private = inf_obj$.__enclos_env__$private
	if (should_run_asserts()) {
		assertChoice(type, c("Zhang"))
		assertList(args_for_type)
		if (!(type %in% names(args_for_type))) stop("args_for_type must contain a list for ", type)
	}
	args = args_for_type[[type]]
	if (should_run_asserts()) {
		assertList(args)
	}
	is_bernoulli = is(private$des_obj, "DesignSeqOneByOneBernoulli") || is(private$des_obj, "DesignFixedBernoulli")
	if (should_run_asserts()) {
		if (!is_bernoulli && !private$has_match_structure) stop("Zhang incidence inference requires Bernoulli or matching designs.")
		assertResponseType(private$des_obj$get_response_type(), "incidence")
		assertNoCensoring(private$any_censoring)
		assertChoice(args$combination_method, c("Fisher", "Stouffer", "min_p"))
		if (!is.null(args$pval_epsilon)) assertNumeric(args$pval_epsilon, lower = .Machine$double.xmin, upper = 1)
	}
	invisible(args)
}

zhang_get_exact_stats = function(inf_obj){
	private = inf_obj$.__enclos_env__$private
	if (!is.null(private$cached_values$incid_exact_zhang_stats)) return(private$cached_values$incid_exact_zhang_stats)
	if (private$has_match_structure){
		if (is.null(private$cached_values$KKstats)){
			private$m = private$des_obj_priv_int$m
			m_vec = if (is.null(private$m)) rep(0, private$n) else private$m
			m_vec[is.na(m_vec)] = 0
			private$cached_values$KKstats = compute_zhang_match_data_cpp(private$get_X(), private$y, private$w, m_vec)
		}
		KKstats = private$cached_values$KKstats
		exact_stats = list(
			m = as.integer(KKstats$m),
			nRT = as.integer(KKstats$nRT),
			nRC = as.integer(KKstats$nRC),
			d_plus = as.integer(KKstats$d_plus),
			d_minus = as.integer(KKstats$d_minus),
			n11 = as.integer(KKstats$n11),
			n10 = as.integer(KKstats$n10),
			n01 = as.integer(KKstats$n01),
			n00 = as.integer(KKstats$n00)
		)
	} else {
		nRT = sum(private$w == 1L, na.rm = TRUE)
		nRC = sum(private$w == 0L, na.rm = TRUE)
		n11 = sum(private$y[private$w == 1L])
		n01 = sum(private$y[private$w == 0L])
		exact_stats = list(
			m = 0L,
			nRT = as.integer(nRT),
			nRC = as.integer(nRC),
			d_plus = 0L,
			d_minus = 0L,
			n11 = as.integer(n11),
			n10 = as.integer(nRT - n11),
			n01 = as.integer(n01),
			n00 = as.integer(nRC - n01)
		)
	}
	private$cached_values$incid_exact_zhang_stats = exact_stats
	exact_stats
}

zhang_compute_exact_pval_matched_pairs = function(inf_obj, delta_0) {
	private = inf_obj$.__enclos_env__$private
	if (!private$has_match_structure) return(NA_real_)
	exact_stats = zhang_get_exact_stats(inf_obj)
	if (exact_stats$m == 0L || exact_stats$d_plus + exact_stats$d_minus == 0L) return(NA_real_)
	zhang_exact_binom_pval_cpp(exact_stats$d_plus, exact_stats$d_minus, delta_0)
}

zhang_compute_exact_pval_reservoir = function(inf_obj, delta_0){
	exact_stats = zhang_get_exact_stats(inf_obj)
	if (exact_stats$nRT == 0L || exact_stats$nRC == 0L) return(NA_real_)
	if (exact_stats$n11 + exact_stats$n01 == 0L || exact_stats$n10 + exact_stats$n00 == 0L) return(NA_real_)
	zhang_exact_fisher_pval_cpp(exact_stats$n11, exact_stats$n10, exact_stats$n01, exact_stats$n00, delta_0)
}

zhang_pval_exact_combined = function(inf_obj, delta_0, combination_method = "Fisher"){
	exact_stats = zhang_get_exact_stats(inf_obj)
	p_M = if (exact_stats$m > 0) zhang_compute_exact_pval_matched_pairs(inf_obj, delta_0) else NA_real_
	p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) zhang_compute_exact_pval_reservoir(inf_obj, delta_0) else NA_real_
	zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
}

zhang_compute_ci_bounds_parallel = function(inf_obj, est, lo_bound, hi_bound, alpha, pval_epsilon, combination_method){
	private = inf_obj$.__enclos_env__$private
	bound_specs = list(list(inside = est, outside = lo_bound), list(inside = est, outside = hi_bound))
	n_cores_ci = min(2L, inf_obj$num_cores)
	child_budget = max(1L, as.integer(floor(inf_obj$num_cores / n_cores_ci)))
	inf_template = inf_obj$duplicate()
	results = private$par_lapply(bound_specs, function(spec){
		worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
		exact_stats = zhang_get_exact_stats(worker_inf)
		p_fn = function(delta_0){
			p_M = if (exact_stats$m > 0) zhang_compute_exact_pval_matched_pairs(worker_inf, delta_0) else NA_real_
			p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) zhang_compute_exact_pval_reservoir(worker_inf, delta_0) else NA_real_
			zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
		}
		zhang_bisect_ci_boundary(p_fn, inside = spec$inside, outside = spec$outside, pval_th = alpha, tol = pval_epsilon)
	}, n_cores = n_cores_ci, budget = child_budget,
	export_list = list(
		inf_template = inf_template,
		combination_method = combination_method,
		alpha = alpha,
		pval_epsilon = pval_epsilon
	))
	c(results[[1]], results[[2]])
}

zhang_ci_exact_combined = function(inf_obj, alpha, pval_epsilon, combination_method = "Fisher"){
	exact_stats = zhang_get_exact_stats(inf_obj)
	est = zhang_incid_treatment_estimate(exact_stats)
	if (should_run_asserts()) {
		if (!is.finite(est)) stop("Cannot compute exact CI: point estimate is not finite.")
	}
	mle_ci = zhang_incid_mle_ci(exact_stats, alpha * 2)
	ci_width = mle_ci[2] - mle_ci[1]
	lo_bound = mle_ci[1] - 0.5 * ci_width
	hi_bound = mle_ci[2] + 0.5 * ci_width
	if (inf_obj$num_cores > 1L) {
		zhang_compute_ci_bounds_parallel(inf_obj, est, lo_bound, hi_bound, alpha, pval_epsilon, combination_method)
	} else {
		p_fn = function(delta_0){
			p_M = if (exact_stats$m > 0) zhang_compute_exact_pval_matched_pairs(inf_obj, delta_0) else NA_real_
			p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) zhang_compute_exact_pval_reservoir(inf_obj, delta_0) else NA_real_
			zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
		}
		c(
			zhang_bisect_ci_boundary(p_fn, inside = est, outside = lo_bound, pval_th = alpha, tol = pval_epsilon),
			zhang_bisect_ci_boundary(p_fn, inside = est, outside = hi_bound, pval_th = alpha, tol = pval_epsilon)
		)
	}
}
