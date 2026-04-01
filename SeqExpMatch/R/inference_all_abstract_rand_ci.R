#' Randomization-based Confidence Intervals
#'
#' Abstract class for randomization-based confidence interval inference.
#'
#' @keywords internal
InferenceRandCI = R6::R6Class("InferenceRandCI",
	lock_objects = FALSE,
	inherit = InferenceRand,
	public = list(
		#' @description
		#' Compute a randomization-based two-sided p-value for the treatment effect.
		#' @param r Number of randomization vectors.
		#' @param delta Null treatment effect value.
		#' @param transform_responses Response transformation to apply during the test.
		#' @param na.rm Whether to remove non-finite simulated statistics.
		#' @param show_progress Whether to show progress.
		#' @param permutations Optional pre-generated assignment draws.
		#' @param type Optional incidence-specific exact randomization type.
		#' @param args_for_type Optional arguments keyed by \code{type}.
		#' @return A two-sided p-value.
		compute_two_sided_pval_for_treatment_effect_rand = function(r = 501, delta = 0, transform_responses = "none", na.rm = TRUE, show_progress = TRUE, permutations = NULL, type = NULL, args_for_type = NULL){
			private$assert_design_supports_resampling("Randomization inference")
			assertLogical(na.rm)
			if (private$des_obj_priv_int$response_type == "incidence" && is.null(private$custom_randomization_statistic_function)){
				if (!identical(transform_responses, "none")) {
					stop("transform_responses is not supported for incidence randomization inference.")
				}
				rand_type = if (is.null(type)) "Zhang" else type
				exact_args = private$normalize_exact_inference_args(
					rand_type,
					args_for_type = args_for_type
				)
				return(private$compute_exact_two_sided_pval_rand(rand_type, delta, exact_args))
			}
			private$assert_no_incidence_only_randomization_args(private$des_obj_priv_int$response_type, type, args_for_type)
			super$compute_two_sided_pval_for_treatment_effect_rand(
				r = r,
				delta = delta,
				transform_responses = transform_responses,
				na.rm = na.rm,
				show_progress = show_progress,
				permutations = permutations
			)
		},

		#' @description
		#' Computes a randomization-based confidence interval.
		#' @param alpha					Significance level.
		#' @param r		Number of randomization vectors.
		#' @param pval_epsilon			Bisection tolerance.
		#' @param show_progress		Show progress.
		#' @param type Optional incidence-specific exact randomization type.
		#' @param args_for_type Optional arguments keyed by \code{type}.
		#' @return 	Randomization CI.
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, type = NULL, args_for_type = NULL){
			private$assert_design_supports_resampling("Randomization inference")
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(r, positive = TRUE); assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			assertLogical(show_progress); show_progress = isTRUE(show_progress) && self$num_cores == 1

			resp_type = private$des_obj_priv_int$response_type
			if (resp_type == "incidence" && is.null(private$custom_randomization_statistic_function)){
				rand_type = if (is.null(type)) "Zhang" else type
				exact_args = private$normalize_exact_inference_args(
					rand_type,
					args_for_type = args_for_type,
					pval_epsilon = pval_epsilon
				)
				return(private$compute_exact_confidence_interval_rand(rand_type, alpha, exact_args))
			}
			private$assert_no_incidence_only_randomization_args(resp_type, type, args_for_type)

			is_glm = inherits(self, "InferenceMLEorKMforGLMs") || 
			         inherits(self, "InferenceAbstractKKGEE") || 
			         inherits(self, "InferenceAbstractKKGLMM") ||
			         inherits(self, "InferenceKKPassThrough") ||
			         inherits(self, "InferencePropUniFractionalLogit") ||
			         inherits(self, "InferencePropZeroOneInflatedBetaAbstract") ||
			         inherits(self, "InferencePropGCompAbstract") ||
			         inherits(self, "InferenceCountZeroAugmentedPoissonAbstract") ||
			         inherits(self, "InferenceCountHurdleNegBinAbstract")
			temp_inf = if (resp_type %in% c("count", "proportion", "survival")) self$duplicate() else self
			transform_arg = "none"
			
			if (resp_type == "count"){
				transform_arg = if (is_glm) "log" else "already_transformed"
				if (!is_glm) temp_inf$.__enclos_env__$private$y = log1p(temp_inf$.__enclos_env__$private$y)
			} else if (resp_type == "proportion"){
				transform_arg = if (is_glm) "logit" else "already_transformed"
				y_clamped = pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, temp_inf$.__enclos_env__$private$y))
				temp_inf$.__enclos_env__$private$y = if (is_glm) y_clamped else logit(y_clamped)
			} else if (resp_type == "survival"){
				transform_arg = if (is_glm) "log" else "already_transformed"
				if (!is_glm) temp_inf$.__enclos_env__$private$y = log(pmax(.Machine$double.eps, temp_inf$.__enclos_env__$private$y))
			}
			if (resp_type %in% c("count", "proportion", "survival")) temp_inf$.__enclos_env__$private$cached_values = list()

			perms = temp_inf$.__enclos_env__$private$generate_permutations(r)
			bounds = private$build_randomization_ci_search_bounds(temp_inf, r, alpha, transform_arg, perms)

				# Run the lower and upper bounds sequentially and reserve all available
				# cores for the inner p-value computations inside each bisection step.
				ci = c(
					temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
						r, bounds$l, bounds$est, alpha / 2, pval_epsilon, transform_arg, TRUE, show_progress, perms
					),
					temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
						r, bounds$est, bounds$u, alpha / 2, pval_epsilon, transform_arg, FALSE, show_progress, perms
					)
				)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),

	private = list(
		assert_no_incidence_only_randomization_args = function(resp_type, type, args_for_type){
			if (!is.null(type)) {
				stop("Randomization type dispatch is only supported for incidence outcomes.")
			}
			if (!is.null(args_for_type)) {
				stop("args_for_type is only used for incidence randomization inference.")
			}
			invisible(NULL)
		},

		normalize_exact_inference_args = function(type, args_for_type = NULL, pval_epsilon = NULL){
			assertChoice(type, c("Zhang"))
			assertList(args_for_type, null.ok = TRUE)
			default_args = list(combination_method = "Fisher")
			if (!is.null(pval_epsilon)) default_args$pval_epsilon = pval_epsilon
			utils::modifyList(setNames(list(default_args), type), if (is.null(args_for_type)) list() else args_for_type)
		},

		assert_exact_inference_params = function(type, args_for_type){
			assertChoice(type, c("Zhang"))
			assertList(args_for_type)
			if (!(type %in% names(args_for_type))) stop("args_for_type must contain a list for ", type)
			args = args_for_type[[type]]
			assertList(args)
				is_bernoulli = is(private$des_obj, "DesignSeqOneByOneBernoulli") || is(private$des_obj, "FixedDesignBernoulli")
				if (!is_bernoulli && !private$has_match_structure) stop("Zhang randomization inference requires Bernoulli or matching designs.")
			assertResponseType(private$des_obj$get_response_type(), "incidence")
			assertNoCensoring(private$any_censoring)
			assertChoice(args$combination_method, c("Fisher", "Stouffer", "min_p"))
			if (!is.null(args$pval_epsilon)) assertNumeric(args$pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			invisible(args)
		},

		compute_exact_confidence_interval_rand = function(type, alpha, args_for_type){
			private$assert_design_supports_resampling("Randomization inference")
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$assert_exact_inference_params(type, args_for_type)
			switch(type,
				Zhang = private$ci_exact_zhang_combined(
					alpha = alpha,
					pval_epsilon = args_for_type[[type]]$pval_epsilon,
					combination_method = args_for_type[[type]]$combination_method
				)
			)
		},

		compute_exact_two_sided_pval_rand = function(type, delta, args_for_type){
			private$assert_design_supports_resampling("Randomization inference")
			assertNumeric(delta)
			private$assert_exact_inference_params(type, args_for_type)
			switch(type,
				Zhang = private$pval_exact_zhang_combined(
					delta_0 = delta,
					combination_method = args_for_type[[type]]$combination_method
				)
			)
		},

		pval_exact_zhang_combined = function(delta_0, combination_method = "Fisher"){
			exact_stats = private$get_exact_zhang_stats()
			p_M = if (exact_stats$m > 0) private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
			p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) private$compute_exact_pval_reservoir(delta_0) else NA_real_
			zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
		},

		compute_exact_pval_matched_pairs = function(delta_0) {
				if (!private$has_match_structure) return(NA_real_)
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
				if (private$has_match_structure){
				if (is.null(private$cached_values$KKstats)){
					private$m = private$des_obj_priv_int$m
					m_vec = if (is.null(private$m)) rep(0, private$n) else private$m
					m_vec[is.na(m_vec)] = 0
					private$cached_values$KKstats = compute_zhang_match_data_cpp(private$w, m_vec, private$y, private$get_X())
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
		},

		ci_exact_zhang_combined = function(alpha, pval_epsilon, combination_method = "Fisher"){
			exact_stats = private$get_exact_zhang_stats()
			est = zhang_incid_treatment_estimate(exact_stats)
			if (!is.finite(est)) stop("Cannot compute exact CI: point estimate is not finite.")
			mle_ci = zhang_incid_mle_ci(exact_stats, alpha * 2)
			ci_width = mle_ci[2] - mle_ci[1]
			lo_bound = mle_ci[1] - 0.5 * ci_width
			hi_bound = mle_ci[2] + 0.5 * ci_width
			if (self$num_cores > 1L) {
				private$compute_zhang_ci_bounds_parallel(est, lo_bound, hi_bound, alpha, pval_epsilon, combination_method)
			} else {
				p_fn = function(delta_0){
					p_M = if (exact_stats$m > 0) private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
					p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) private$compute_exact_pval_reservoir(delta_0) else NA_real_
					zhang_combine_exact_pvals(p_M, p_R, exact_stats$m, exact_stats$nRT, exact_stats$nRC, combination_method)
				}
				c(
					zhang_bisect_ci_boundary(p_fn, inside = est, outside = lo_bound, pval_th = alpha, tol = pval_epsilon),
					zhang_bisect_ci_boundary(p_fn, inside = est, outside = hi_bound, pval_th = alpha, tol = pval_epsilon)
				)
			}
		},

			compute_zhang_ci_bounds_parallel = function(est, lo_bound, hi_bound, alpha, pval_epsilon, combination_method){
				bound_specs = list(list(inside = est, outside = lo_bound), list(inside = est, outside = hi_bound))
				n_cores_ci = min(2L, self$num_cores)
				child_budget = max(1L, as.integer(floor(self$num_cores / n_cores_ci)))
				inf_template = self$duplicate()
				results = private$par_lapply(bound_specs, function(spec){
					worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
					worker_private = worker_inf$.__enclos_env__$private
					exact_stats = worker_private$get_exact_zhang_stats()
					p_fn = function(delta_0){
						p_M = if (exact_stats$m > 0) worker_private$compute_exact_pval_matched_pairs(delta_0) else NA_real_
						p_R = if (exact_stats$nRT > 0 && exact_stats$nRC > 0) worker_private$compute_exact_pval_reservoir(delta_0) else NA_real_
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
			},

		build_randomization_ci_search_bounds = function(inf_obj, r, alpha, transform_arg, permutations){
			obj_private = inf_obj$.__enclos_env__$private
			if (transform_arg == "none" && (is.null(obj_private$cached_values$t0s_rand) || length(obj_private$cached_values$t0s_rand) < r)) {
				inf_obj$compute_two_sided_pval_for_treatment_effect_rand(r, 0, transform_arg, FALSE, show_progress = FALSE, permutations = permutations)
			}
			est = as.numeric(inf_obj$compute_treatment_estimate())
			if (length(est) == 0L || !is.finite(est[1])) est = NA_real_ else est = est[1]
			asym_ci = tryCatch(as.numeric(inf_obj$compute_asymp_confidence_interval(alpha = alpha * 2)), error = function(e) c(NA_real_, NA_real_))
			if (length(asym_ci) < 2L || !all(is.finite(asym_ci[1:2]))) asym_ci = c(NA_real_, NA_real_) else asym_ci = sort(asym_ci[1:2])
			if (!is.finite(est) && all(is.finite(asym_ci))) est = mean(asym_ci)
			if (!all(is.finite(asym_ci))) {
				scale_guess = stats::sd(obj_private$y, na.rm = TRUE)
				if (!is.finite(scale_guess) || scale_guess <= 0) scale_guess = stats::IQR(obj_private$y, na.rm = TRUE) / 1.349
				if (!is.finite(scale_guess) || scale_guess <= 0) scale_guess = 1
				if (!is.finite(est)) est = stats::median(obj_private$y, na.rm = TRUE)
				if (!is.finite(est)) est = 0
				asym_ci = c(est - 2 * scale_guess, est + 2 * scale_guess)
			}
			l = asym_ci[1]; u = asym_ci[2]
			if (l >= est) l = est - max(abs(u - est), 1); if (u <= est) u = est + max(abs(est - l), 1)
			list(est = est, l = private$expand_bound(inf_obj, l, est, r, transform_arg, permutations, alpha / 2, TRUE), u = private$expand_bound(inf_obj, u, est, r, transform_arg, permutations, alpha / 2, FALSE))
		},

		expand_bound = function(inf_obj, bound, est, r, transform_arg, permutations, target_pval, lower){
			evaluate_pval = function(delta) {
				pval = tryCatch(as.numeric(inf_obj$compute_two_sided_pval_for_treatment_effect_rand(r, delta, transform_arg, FALSE, FALSE, permutations)), error = function(e) NA_real_)
				if (length(pval) == 0L) return(NA_real_) else pval[1]
			}
			pval_bound = evaluate_pval(bound)
			if (is.finite(pval_bound) && pval_bound < target_pval) return(bound)
			step = if (is.finite(abs(est - bound)) && abs(est - bound) > 0) abs(est - bound) else 1
			for (iter in seq_len(12L)) {
				step = step * 2; candidate = if (lower) est - step else est + step
				if (is.finite(evaluate_pval(candidate)) && evaluate_pval(candidate) < target_pval) { bound = candidate; break }
				bound = candidate
			}
			bound
		},

			compute_ci_by_inverting_the_randomization_test_iteratively = function(r, l, u, pval_th, tol, transform_responses, lower, show_progress = TRUE, permutations = NULL){
			pval_l = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(r, delta = l, transform_responses = transform_responses, show_progress = FALSE, permutations = permutations))
			pval_u = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(r, delta = u, transform_responses = transform_responses, show_progress = FALSE, permutations = permutations))
			if (length(pval_l) == 0) pval_l = NA_real_; if (length(pval_u) == 0) pval_u = NA_real_
			for (k in seq_len(30L)) {
				if (!is.na(pval_l) && !is.na(pval_u)) break
				if (is.na(pval_l)) { l = (l + u) / 2; pval_l = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(r, delta = l, transform_responses = transform_responses, show_progress = FALSE, permutations = permutations)) }
				if (is.na(pval_u)) { u = (l + u) / 2; pval_u = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(r, delta = u, transform_responses = transform_responses, show_progress = FALSE, permutations = permutations)) }
			}
			if (is.na(pval_l) || is.na(pval_u) || !all(is.finite(c(l, u)))) return(NA_real_)
			iter = 0; progress_label = if (lower) "CI lower" else "CI upper"
			repeat {
				pval_span = abs(pval_u - pval_l)
				if ((abs(u - l)) <= tol || pval_span <= tol) {
					if (isTRUE(show_progress)) cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g) done\n", progress_label, iter, pval_span, tol))
					return(if(lower) l else u)
				}
				m = (l + u) / 2.0; pval_m = as.numeric(self$compute_two_sided_pval_for_treatment_effect_rand(r, delta = m, transform_responses = transform_responses, show_progress = FALSE, permutations = permutations))
				if (is.na(pval_m)) { if (lower) { l = m; pval_l = 0 } else { u = m; pval_u = 0 }; iter = iter + 1; next }
				if (pval_m >= pval_th && lower) { u = m; pval_u = pval_m }
				else if (pval_m >= pval_th && !lower) { l = m; pval_l = pval_m }
				else if (lower) { l = m; pval_l = pval_m }
				else { u = m; pval_u = pval_m }
				iter = iter + 1
				if (isTRUE(show_progress)) cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g)", progress_label, iter, pval_span, tol))
			}
		}
	)
)


#' Exact Zhang Incidence Inference
#'
#' @export
InferenceIncidExactZhang = R6::R6Class("InferenceIncidExactZhang",
	lock_objects = FALSE,
	inherit = InferenceRandCI,
	public = list(
		#' @description
		#' Initialize exact Zhang incidence inference.
		#' @param des_obj A completed design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceIncidExactZhang} object.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Compute the Zhang incidence treatment estimate.
		#' @param estimate_only Ignored for this estimator.
		#' @return The treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			stats = private$get_exact_zhang_stats()
			zhang_incid_treatment_estimate(stats)
		},

		#' @description
		#' Compute an exact Zhang confidence interval.
		#' @param alpha Significance level.
		#' @param pval_epsilon Bisection tolerance for the inversion routine.
		#' @param args_for_type Optional arguments keyed by exact type.
		#' @return A confidence interval.
		compute_exact_confidence_interval = function(alpha = 0.05, pval_epsilon = 0.005, args_for_type = NULL){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			exact_args = private$normalize_exact_inference_args(
				"Zhang",
				args_for_type = args_for_type,
				pval_epsilon = pval_epsilon
			)
			private$compute_exact_confidence_interval_rand("Zhang", alpha, exact_args)
		},

		#' @description
		#' Compute the exact Zhang two-sided p-value for a treatment effect.
		#' @param delta Null treatment effect value.
		#' @param args_for_type Optional arguments keyed by exact type.
		#' @return A two-sided p-value.
		compute_exact_two_sided_pval_for_treatment_effect = function(delta = 0, args_for_type = NULL){
			assertNumeric(delta, len = 1)
			exact_args = private$normalize_exact_inference_args(
				"Zhang",
				args_for_type = args_for_type
			)
			private$compute_exact_two_sided_pval_rand("Zhang", delta, exact_args)
		}
	)
)
