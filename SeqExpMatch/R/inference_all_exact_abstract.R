# Exact Inference via Zhang's Method
#
# @description
# Abstract class for Zhang (2026) exact test-inversion inference.
#
# @keywords internal
InferenceExact = R6::R6Class("InferenceExact",
	inherit = Inference,
	public = list(
		# @description
		# Computes an exact confidence interval via Zhang's combined test.
		# @param type             The type of exact inference. Currently only "Zhang" is supported.
		# @param alpha            Significance level.
		# @param args_for_type    A list of parameters including combination_method and pval_epsilon.
		compute_exact_confidence_interval = function(type = "Zhang", alpha = 0.05, args_for_type = list(Zhang = list(combination_method = "Fisher", pval_epsilon = 0.005))){
			private$assert_exact_inference_params(type, args_for_type)
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			if (type == "Zhang") private$ci_exact_zhang_combined(alpha, args_for_type[[type]]$pval_epsilon, args_for_type[[type]]$combination_method)
		},

		# @description
		# Computes an exact two-sided p-value via Zhang's combined test.
		# @param type             The type of exact inference.
		# @param delta            Null treatment effect (log-odds ratio).
		# @param args_for_type    A list containing combination_method.
		compute_exact_two_sided_pval_for_treatment_effect = function(type = "Zhang", delta = 0, args_for_type = list(Zhang = list(combination_method = "Fisher"))){
			private$assert_exact_inference_params(type, args_for_type); assertNumeric(delta)
			if (type == "Zhang") private$pval_exact_zhang_combined(delta, args_for_type[[type]]$combination_method)
		}
	),

	private = list(
		assert_exact_inference_params = function(type, args_for_type){
			assertChoice(type, c("Zhang")); assertList(args_for_type)
			if (!(type %in% names(args_for_type))) stop("args_for_type must contain a list for ", type)
			args = args_for_type[[type]]
			if (type == "Zhang") {
				is_bernoulli = is(private$des_obj, "DesignSeqOneByOneBernoulli") || is(private$des_obj, "FixedDesignBernoulli")
				if (!is_bernoulli && !private$is_KK) stop("Zhang inference requires Bernoulli or matching designs.")
				assertResponseType(private$des_obj$get_response_type(), "incidence"); assertNoCensoring(private$any_censoring)
				assertChoice(args$combination_method, c("Fisher", "Stouffer", "min_p"))
				if (!is.null(args$pval_epsilon)) assertNumeric(args$pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			}
		},

		pval_exact_zhang_combined = function(delta_0, combination_method = "Fisher"){
			exact_stats = private$get_exact_zhang_stats()
			p_M = if (exact_stats$m > 0) private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
			p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) private$compute_exact_pval_reservoir(delta_0) else NA_real_
			zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
		},

		compute_exact_pval_matched_pairs = function(delta_0) {
			if (!private$is_KK) return(NA_real_)
			exact_stats = private$get_exact_zhang_stats()
			if (exact_stats$m == 0L || exact_stats$d_plus + exact_stats$d_minus == 0L) return(NA_real_)
			zhang_exact_binom_pval_cpp(exact_stats$d_plus, exact_stats$d_minus, delta_0)
		},

		compute_exact_pval_reservoir = function(delta_0){
			exact_stats = private$get_exact_zhang_stats()
			if (exact_stats$nRT == 0L || exact_stats$nRC == 0L) return(NA_real_)
			if (exact_stats$n11 + exact_stats$n01 == 0L || exact_stats$n10 + exact_stats$n00 == 0L) return(NA_real_)
			zhang_exact_fisher_pval_cpp(exact_stats$n11, exact_stats$n10, exact_stats$n01, exact_stats$n00, delta_0)
		},

		get_exact_zhang_stats = function(){
			if (!is.null(private$cached_values$incid_exact_zhang_stats)) return(private$cached_values$incid_exact_zhang_stats)
			if (private$is_KK){
				if (is.null(private$cached_values$KKstats)){
					private$m = private$des_obj_priv_int$m
					m_vec = if (is.null(private$m)) rep(0, private$n) else private$m
					m_vec[is.na(m_vec)] = 0
					private$cached_values$KKstats = compute_zhang_match_data_cpp(private$w, m_vec, private$y, private$des_obj$get_X())
				}
				KKstats = private$cached_values$KKstats
				exact_stats = list(m = as.integer(KKstats$m), nRT = as.integer(KKstats$nRT), nRC = as.integer(KKstats$nRC), d_plus = as.integer(KKstats$d_plus), d_minus = as.integer(KKstats$d_minus), n11 = as.integer(KKstats$n11), n10 = as.integer(KKstats$n10), n01 = as.integer(KKstats$n01), n00 = as.integer(KKstats$n00))
			} else {
				nRT = sum(private$w == 1L, na.rm = TRUE); nRC = sum(private$w == 0L, na.rm = TRUE)
				n11 = sum(private$y[private$w == 1L]); n01 = sum(private$y[private$w == 0L])
				exact_stats = list(m = 0L, nRT = as.integer(nRT), nRC = as.integer(nRC), d_plus = 0L, d_minus = 0L, n11 = as.integer(n11), n10 = as.integer(nRT - n11), n01 = as.integer(n01), n00 = as.integer(nRC - n01))
			}
			private$cached_values$incid_exact_zhang_stats = exact_stats; exact_stats
		},

		ci_exact_zhang_combined = function(alpha, pval_epsilon, combination_method = "Fisher"){
			exact_stats = private$get_exact_zhang_stats(); est = zhang_incid_treatment_estimate(exact_stats)
			if (!is.finite(est)) stop("Cannot compute exact CI: point estimate is not finite.")
			mle_ci = zhang_incid_mle_ci(exact_stats, alpha * 2); ci_width = mle_ci[2] - mle_ci[1]
			lo_bound = mle_ci[1] - 0.5 * ci_width; hi_bound = mle_ci[2] + 0.5 * ci_width
			if (private$num_cores > 1L) private$compute_zhang_ci_bounds_parallel(est, lo_bound, hi_bound, alpha, pval_epsilon, combination_method)
			else {
				p_fn = function(delta_0){
					p_M = if (exact_stats$m > 0) private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
					p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) private$compute_exact_pval_reservoir(delta_0) else NA_real_
					zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
				}
				c(zhang_bisect_ci_boundary(p_fn, inside = est, outside = lo_bound, pval_th = alpha, tol = pval_epsilon), zhang_bisect_ci_boundary(p_fn, inside = est, outside = hi_bound, pval_th = alpha, tol = pval_epsilon))
			}
		},

		compute_zhang_ci_bounds_parallel = function(est, lo_bound, hi_bound, alpha, pval_epsilon, combination_method){
			bound_specs = list(list(inside = est, outside = lo_bound), list(inside = est, outside = hi_bound))
			child_budget = max(1L, as.integer(floor(private$num_cores / 2))); inf_template = self$duplicate()
			results = parallel::mclapply(bound_specs, function(spec){
				worker_inf = inf_template$duplicate(); worker_inf$.__enclos_env__$private$num_cores = child_budget; set_package_threads(child_budget)
				worker_private = worker_inf$.__enclos_env__$private; exact_stats = worker_private$get_exact_zhang_stats()
				p_fn = function(delta_0){
					p_M = if (exact_stats$m > 0) worker_private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
					p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) worker_private$compute_exact_pval_reservoir(delta_0) else NA_real_
					zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
				}
				zhang_bisect_ci_boundary(p_fn, inside = spec$inside, outside = spec$outside, pval_th = alpha, tol = pval_epsilon)
			}, mc.cores = min(2L, private$num_cores))
			c(results[[1]], results[[2]])
		}
	)
)
