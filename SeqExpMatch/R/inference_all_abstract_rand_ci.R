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
		#' @param ci_search_control Optional control list for randomization-CI search. Supported
		#'   entries are \code{fallback}, \code{seed}, \code{max_radius_se_mult},
		#'   \code{max_radius_scale_mult}, \code{max_expansions}, \code{seed_boot_B},
		#'   Monte Carlo settings \code{mc_enable}, \code{mc_batch_size},
		#'   \code{mc_min_draws}, and \code{mc_conf_level}, midpoint-cache settings
		#'   \code{pval_cache_enable} and \code{pval_cache_resolution}, and model-fit
		#'   reuse settings \code{fit_warm_start_enable} and
		#'   \code{fit_reuse_factorizations}. Set \code{mc_enable = FALSE} to force
		#'   full enumeration of all requested randomization draws.
		#' @return 	Randomization CI.
		compute_confidence_interval_rand = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, type = NULL, args_for_type = NULL, ci_search_control = NULL){
			private$assert_design_supports_resampling("Randomization inference")
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertCount(r, positive = TRUE); assertNumeric(pval_epsilon, lower = .Machine$double.xmin, upper = 1)
			assertLogical(show_progress); show_progress = isTRUE(show_progress) && self$num_cores == 1
			ci_search_control = private$normalize_randomization_ci_search_control(ci_search_control, r, pval_epsilon)
			ci_search_control$mc_stop_threshold = alpha / 2

			dispatch_cores = private$effective_parallel_cores("rand_ci", self$num_cores)
			if (dispatch_cores != self$num_cores) {
				serial_inf = self$duplicate(verbose = private$verbose)
				serial_inf$num_cores = dispatch_cores
				return(serial_inf$compute_confidence_interval_rand(
					alpha = alpha,
					r = r,
					pval_epsilon = pval_epsilon,
					show_progress = show_progress,
					type = type,
					args_for_type = args_for_type,
					ci_search_control = ci_search_control
				))
			}

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
			old_mc_control = temp_inf$.__enclos_env__$private$randomization_mc_control
			temp_inf$.__enclos_env__$private$randomization_mc_control = ci_search_control
			on.exit({ temp_inf$.__enclos_env__$private$randomization_mc_control = old_mc_control }, add = TRUE)

			perms = temp_inf$.__enclos_env__$private$generate_permutations(r)
			ci_pval_cache = if (isTRUE(ci_search_control$pval_cache_enable)) new.env(parent = emptyenv()) else NULL
			bounds = private$build_randomization_ci_search_bounds(temp_inf, r, alpha, transform_arg, perms, ci_search_control, ci_pval_cache)
			if (!all(is.finite(c(bounds$l, bounds$u)))) {
				fallback_ci = as.numeric(bounds$fallback_ci)
				fallback_mode = ci_search_control$fallback
				if (identical(fallback_mode, "error")) {
					stop("Randomization CI search failed to bracket the target p-value within the configured search radius.")
				}
				if (identical(fallback_mode, "na") || length(fallback_ci) < 2L || !all(is.finite(fallback_ci[1:2]))) {
					ci = c(NA_real_, NA_real_)
				} else {
					ci = sort(fallback_ci[1:2])
				}
				names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
				return(ci)
			}

			# Run the lower and upper bounds sequentially and reserve all available
			# cores for the inner p-value computations inside each bisection step.
			ci = c(
				temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
					r, bounds$l, bounds$est, alpha / 2, pval_epsilon, transform_arg, TRUE, show_progress, perms, ci_search_control, ci_pval_cache
				),
				temp_inf$.__enclos_env__$private$compute_ci_by_inverting_the_randomization_test_iteratively(
					r, bounds$est, bounds$u, alpha / 2, pval_epsilon, transform_arg, FALSE, show_progress, perms, ci_search_control, ci_pval_cache
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

		compute_randomization_ci_pval_cached = function(inf_obj, r, delta, transform_responses, permutations, ci_search_control, ci_pval_cache){
			cache_enabled = isTRUE(ci_search_control$pval_cache_enable) && !is.null(ci_pval_cache)
			if (cache_enabled) {
				cache_key = private$normalize_delta_for_cache(delta, ci_search_control$pval_cache_resolution)
				if (!is.null(ci_pval_cache[[cache_key]])) return(ci_pval_cache[[cache_key]])
			}
			pval = tryCatch(
				as.numeric(inf_obj$compute_two_sided_pval_for_treatment_effect_rand(
					r,
					delta = delta,
					transform_responses = transform_responses,
					na.rm = FALSE,
					show_progress = FALSE,
					permutations = permutations
				)),
				error = function(e) NA_real_
			)
			if (length(pval) == 0L) pval = NA_real_ else pval = pval[1]
			if (cache_enabled) ci_pval_cache[[cache_key]] = pval
			pval
		},

		normalize_randomization_ci_search_control = function(ci_search_control, r, pval_epsilon){
			default_mc_batch = min(as.integer(r), max(25L, as.integer(ceiling(2 * sqrt(as.integer(r))))))
			defaults = list(
				fallback = "fallback",
				seed = "asymp_then_boot",
				max_radius_se_mult = 12,
				max_radius_scale_mult = 6,
				max_expansions = 5L,
				seed_boot_B = max(51L, min(as.integer(r), 201L)),
				pval_cache_enable = TRUE,
				pval_cache_resolution = pval_epsilon,
				mc_enable = as.integer(r) >= 200L,
				mc_batch_size = default_mc_batch,
				mc_min_draws = min(as.integer(r), max(100L, 2L * default_mc_batch)),
				mc_conf_level = 0.99,
				fit_warm_start_enable = TRUE,
				fit_reuse_factorizations = TRUE
			)
			assertList(ci_search_control, null.ok = TRUE)
			ctrl = utils::modifyList(defaults, if (is.null(ci_search_control)) list() else ci_search_control)
			assertChoice(ctrl$fallback, c("fallback", "na", "error"))
			assertChoice(ctrl$seed, c("asymp_then_boot", "boot_then_asymp", "asymp_only", "boot_only", "none"))
			assertNumber(ctrl$max_radius_se_mult, lower = 0, finite = TRUE)
			assertNumber(ctrl$max_radius_scale_mult, lower = 0, finite = TRUE)
			assertCount(as.integer(ctrl$max_expansions), positive = TRUE)
			assertCount(as.integer(ctrl$seed_boot_B), positive = TRUE)
			assertFlag(ctrl$pval_cache_enable)
			assertNumber(ctrl$pval_cache_resolution, lower = .Machine$double.eps, finite = TRUE)
			assertFlag(ctrl$mc_enable)
			assertCount(as.integer(ctrl$mc_batch_size), positive = TRUE)
			assertCount(as.integer(ctrl$mc_min_draws), positive = TRUE)
			assertNumber(ctrl$mc_conf_level, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			assertFlag(ctrl$fit_warm_start_enable)
			assertFlag(ctrl$fit_reuse_factorizations)
			ctrl$max_expansions = as.integer(ctrl$max_expansions)
			ctrl$seed_boot_B = as.integer(ctrl$seed_boot_B)
			ctrl$pval_cache_resolution = as.numeric(ctrl$pval_cache_resolution)
			ctrl$mc_batch_size = min(as.integer(r), as.integer(ctrl$mc_batch_size))
			ctrl$mc_min_draws = min(as.integer(r), max(as.integer(ctrl$mc_batch_size), as.integer(ctrl$mc_min_draws)))
			ctrl
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

		build_randomization_ci_search_bounds = function(inf_obj, r, alpha, transform_arg, permutations, ci_search_control, ci_pval_cache){
			normalize_ci = function(ci){
				ci = as.numeric(ci)
				if (length(ci) < 2L || !all(is.finite(ci[1:2]))) return(c(NA_real_, NA_real_))
				sort(ci[1:2])
			}
			obj_private = inf_obj$.__enclos_env__$private
			if (transform_arg == "none" && (is.null(obj_private$cached_values$t0s_rand) || length(obj_private$cached_values$t0s_rand) < r)) {
				private$compute_randomization_ci_pval_cached(inf_obj, r, 0, transform_arg, permutations, ci_search_control, ci_pval_cache)
			}
			est = as.numeric(inf_obj$compute_treatment_estimate())
			if (length(est) == 0L || !is.finite(est[1])) est = NA_real_ else est = est[1]
			asym_ci = normalize_ci(tryCatch(inf_obj$compute_asymp_confidence_interval(alpha = alpha * 2), error = function(e) c(NA_real_, NA_real_)))
			boot_ci = c(NA_real_, NA_real_)
			need_boot = ci_search_control$seed %in% c("boot_then_asymp", "boot_only") || !all(is.finite(asym_ci))
			if (need_boot) {
				boot_ci = normalize_ci(tryCatch(inf_obj$compute_bootstrap_confidence_interval(alpha = alpha * 2, B = ci_search_control$seed_boot_B, na.rm = TRUE, show_progress = FALSE), error = function(e) c(NA_real_, NA_real_)))
			}
			if (!is.finite(est) && all(is.finite(asym_ci))) est = mean(asym_ci)
			if (!is.finite(est) && all(is.finite(boot_ci))) est = mean(boot_ci)
			response_scale = stats::sd(obj_private$y, na.rm = TRUE)
			if (!is.finite(response_scale) || response_scale <= 0) response_scale = stats::IQR(obj_private$y, na.rm = TRUE) / 1.349
			if (!is.finite(response_scale) || response_scale <= 0) response_scale = 1
			if (!is.finite(est)) est = stats::median(obj_private$y, na.rm = TRUE)
			if (!is.finite(est)) est = 0
			se_guess = NA_real_
			if (all(is.finite(asym_ci))) se_guess = max(abs(asym_ci - est)) / max(stats::qnorm(1 - alpha), .Machine$double.eps)
			if (!is.finite(se_guess) && all(is.finite(boot_ci))) se_guess = max(abs(boot_ci - est)) / max(stats::qnorm(1 - alpha), .Machine$double.eps)
			if (!is.finite(se_guess) || se_guess <= 0) se_guess = response_scale / sqrt(max(1, obj_private$n))
			if (!is.finite(se_guess) || se_guess <= 0) se_guess = response_scale
			default_seed_ci = sort(c(est - 2 * se_guess, est + 2 * se_guess))
			if (!all(is.finite(asym_ci)) && !all(is.finite(boot_ci))) {
				asym_ci = default_seed_ci
			}
			fallback_ci = if (all(is.finite(asym_ci))) asym_ci else if (all(is.finite(boot_ci))) boot_ci else c(NA_real_, NA_real_)
			seed_ci = switch(ci_search_control$seed,
				asymp_then_boot = if (all(is.finite(asym_ci))) asym_ci else if (all(is.finite(boot_ci))) boot_ci else default_seed_ci,
				boot_then_asymp = if (all(is.finite(boot_ci))) boot_ci else if (all(is.finite(asym_ci))) asym_ci else default_seed_ci,
				asymp_only = if (all(is.finite(asym_ci))) asym_ci else c(NA_real_, NA_real_),
				boot_only = if (all(is.finite(boot_ci))) boot_ci else c(NA_real_, NA_real_),
				none = default_seed_ci
			)
			max_radius = max(ci_search_control$max_radius_se_mult * se_guess, ci_search_control$max_radius_scale_mult * response_scale, 1)
			min_radius = min(max(se_guess, 10 * .Machine$double.eps), max_radius)
			l = max(seed_ci[1], est - max_radius)
			u = min(seed_ci[2], est + max_radius)
			if (!is.finite(l) || l >= est) l = est - min_radius
			if (!is.finite(u) || u <= est) u = est + min_radius
			l = private$expand_bound(inf_obj, l, est, r, transform_arg, permutations, alpha / 2, TRUE, max_radius, ci_search_control$max_expansions, ci_search_control, ci_pval_cache)
			u = private$expand_bound(inf_obj, u, est, r, transform_arg, permutations, alpha / 2, FALSE, max_radius, ci_search_control$max_expansions, ci_search_control, ci_pval_cache)
			list(est = est, l = l, u = u, fallback_ci = fallback_ci)
		},

		expand_bound = function(inf_obj, bound, est, r, transform_arg, permutations, target_pval, lower, max_radius, max_expansions, ci_search_control, ci_pval_cache){
			evaluate_pval = function(delta) {
				private$compute_randomization_ci_pval_cached(inf_obj, r, delta, transform_arg, permutations, ci_search_control, ci_pval_cache)
			}
			if (!is.finite(bound) || !is.finite(est) || !is.finite(max_radius) || max_radius <= 0) return(NA_real_)
			bound = if (lower) max(bound, est - max_radius) else min(bound, est + max_radius)
			pval_bound = evaluate_pval(bound)
			if (is.finite(pval_bound) && pval_bound < target_pval) return(bound)
			step = abs(est - bound)
			if (!is.finite(step) || step <= 0) step = min(max_radius / 4, 1)
			for (iter in seq_len(max_expansions)) {
				step = min(step * 2, max_radius)
				candidate = if (lower) est - step else est + step
				pval_candidate = evaluate_pval(candidate)
				if (is.finite(pval_candidate) && pval_candidate < target_pval) return(candidate)
				if (step >= max_radius) break
			}
			NA_real_
		},

			compute_ci_by_inverting_the_randomization_test_iteratively = function(r, l, u, pval_th, tol, transform_responses, lower, show_progress = TRUE, permutations = NULL, ci_search_control = NULL, ci_pval_cache = NULL){
			evaluate_pval = function(delta) {
				private$compute_randomization_ci_pval_cached(self, r, delta, transform_responses, permutations, ci_search_control, ci_pval_cache)
			}
			pval_l = evaluate_pval(l)
			pval_u = evaluate_pval(u)
			if (length(pval_l) == 0) pval_l = NA_real_; if (length(pval_u) == 0) pval_u = NA_real_
			for (k in seq_len(30L)) {
				if (!is.na(pval_l) && !is.na(pval_u)) break
				if (is.na(pval_l)) { l = (l + u) / 2; pval_l = evaluate_pval(l) }
				if (is.na(pval_u)) { u = (l + u) / 2; pval_u = evaluate_pval(u) }
			}
			if (is.na(pval_l) || is.na(pval_u) || !all(is.finite(c(l, u)))) return(NA_real_)
			iter = 0; progress_label = if (lower) "CI lower" else "CI upper"
			repeat {
				pval_span = abs(pval_u - pval_l)
				if ((abs(u - l)) <= tol || pval_span <= tol) {
					if (isTRUE(show_progress)) cat(sprintf("\r%s iter=%d pval_span=%.6g (target<=%.6g) done\n", progress_label, iter, pval_span, tol))
					return(if(lower) l else u)
				}
				m = (l + u) / 2.0; pval_m = evaluate_pval(m)
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
