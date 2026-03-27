#' Randomization-based Confidence Intervals
#'
#' @description
#' Abstract class for randomization-based confidence interval inference.
#'
#' @keywords internal
InferenceRandCI = R6::R6Class("InferenceRandCI",
	inherit = InferenceRand,
	public = list(
		# @description
		# Computes a randomization-based confidence interval.
		# @param alpha					Significance level.
		# @param r		Number of randomization vectors.
		# @param pval_epsilon			Bisection tolerance.
		# @param show_progress		Show progress.
		# @return 	Randomization CI.
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE){
			private$assert_design_supports_resampling("Randomization inference")
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(r, positive = TRUE); assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			assertLogical(show_progress); show_progress = isTRUE(show_progress) && private$num_cores == 1

			resp_type = private$des_obj_priv_int$response_type
			if (resp_type == "incidence") stop("Randomization CI not supported for incidence. Use InferenceIncidExactZhang.")

			is_glm = inherits(self, "InferenceMLEorKMforGLMs") || inherits(self, "InferenceAbstractKKGEE") || inherits(self, "InferenceAbstractKKGLMM")
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

			use_parallel_ci = private$num_cores > 1L && !(private$has_private_method("compute_fast_randomization_distr") && is.null(private[["custom_randomization_statistic_function"]]))

			if (use_parallel_ci) {
				ci = private$compute_ci_both_bounds_parallel(r, bounds$l, bounds$est, bounds$est, bounds$u, alpha / 2, pval_epsilon, transform_arg, perms, inf_obj = temp_inf)
			} else {
				ci = c(
					temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(r, bounds$l, bounds$est, alpha / 2, pval_epsilon, transform_arg, TRUE, show_progress, perms),
					temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(r, bounds$est, bounds$u, alpha / 2, pval_epsilon, transform_arg, FALSE, show_progress, perms)
				)
			}
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),

	private = list(
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

		compute_ci_both_bounds_parallel = function(r, l_lower, u_lower, l_upper, u_upper, pval_th, tol, transform_responses, permutations = NULL, inf_obj = self){
			bound_specs = list(list(l = l_lower, u = u_lower, lower = TRUE), list(l = l_upper, u = u_upper, lower = FALSE))
			inf_template = inf_obj$duplicate()
			results = parallel::mclapply(bound_specs, function(spec) {
				worker_inf = inf_template$duplicate(); worker_inf$.__enclos_env__$private$num_cores = 1L
				set_package_threads(1L)
				worker_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(r, l = spec$l, u = spec$u, pval_th = pval_th, tol = tol, transform_responses = transform_responses, lower = spec$lower, show_progress = FALSE, permutations = permutations)
			}, mc.cores = min(2L, private$num_cores))
			c(results[[1]], results[[2]])
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
