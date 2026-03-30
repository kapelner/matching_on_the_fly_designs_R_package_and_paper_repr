#' Randomization-based Inference
#'
#' Abstract class for randomization-based inference.
#'
#' @keywords internal
InferenceRand = R6::R6Class("InferenceRand",
	inherit = Inference,
	lock_objects = FALSE,
	public = list(
		#' @description
		#' Set Custom Randomization Statistic Computation
		#' @param custom_randomization_statistic_function	A function that returns a scalar value.
		set_custom_randomization_statistic_function = function(custom_randomization_statistic_function){
			assertFunction(custom_randomization_statistic_function, null.ok = TRUE)
			private[["custom_randomization_statistic_function"]] = custom_randomization_statistic_function
			if (!is.null(custom_randomization_statistic_function)){
				environment(private[["custom_randomization_statistic_function"]]) = environment(self$initialize)
			}
			private$cached_values$t0s_rand = NULL
			private$cached_values$rand_distr_cache = list()
			private$cached_values$custom_stat_analysis = NULL
		},

		#' @description
		#' Under the sharp null, computes the randomization distribution.
		#' @param r						The number of randomization vectors.
		#' @param delta					The null difference.
		#' @param transform_responses	Type of transformation.
		#' @param show_progress		Show progress bar.
		#' @param permutations		Pre-computed permutations.
		#' @return 	A vector of estimates.
		compute_beta_hat_T_randomization_distr_under_sharp_null = function(r = 501, delta = 0, transform_responses = "none", show_progress = TRUE, permutations = NULL){
			private$assert_design_supports_resampling("Randomization inference")
			assertNumeric(delta); assertCount(r, positive = TRUE)

			if (is.null(permutations)) permutations = private$generate_permutations(r)
			setup = private$setup_randomization_template_and_shifts(delta, transform_responses)
			
			if (!is.null(permutations) && private$has_private_method("compute_fast_randomization_distr")) {
				fast_distr = private$compute_fast_randomization_distr(setup$template$.__enclos_env__$private$y, permutations, delta, transform_responses)
				if (!is.null(fast_distr)) return(fast_distr)
			}

			custom_stat_analysis = private$analyze_custom_randomization_statistic()
			use_lightweight_custom_stat = isTRUE(custom_stat_analysis$can_use_lightweight_yw_only)
			use_perms = !is.null(permutations) && (!is.null(permutations$w_mat) || length(permutations) >= r)
			
			need_thread_objs = !(use_lightweight_custom_stat && use_perms)
			inf_template = if (need_thread_objs) self$duplicate() else NULL
			des_template = if (need_thread_objs) setup$template$duplicate() else NULL

			# Warm up the design template cache if it is a sequential design that uses covariates.
			if (!is.null(des_template) && isTRUE(des_template$.__enclos_env__$private$uses_covariates)) {
				if (private$verbose) cat("Warming up design cache... ")
				tryCatch({
					priv = des_template$.__enclos_env__$private
					old_t = priv$t
					if (is.null(priv$all_subject_data_cache)) priv$all_subject_data_cache = list()
					n_subjects = des_template$get_n()
					for (t_temp in 1 : n_subjects) {
						priv$t = t_temp
						priv$compute_all_subject_data()
					}
					priv$t = old_t
					if (private$verbose) cat("done.\n")
				}, error = function(e) {
					if (exists("old_t")) des_template$.__enclos_env__$private$t = old_t
					if (private$verbose) cat("failed.\n")
				})
			}

			if (!is.null(inf_template) && private$is_KK && private$object_has_private_method(inf_template, "compute_basic_match_data"))
				inf_template$.__enclos_env__$private$compute_basic_match_data()

			actual_rand_cores = private$num_cores
			if (actual_rand_cores > 1L && need_thread_objs) {
				do_warmup_iter = function() {
					w_des = if (!is.null(des_template)) des_template$duplicate() else NULL
					w_inf = if (!is.null(inf_template)) inf_template$duplicate() else NULL
					if (!is.null(w_inf)) w_inf$.__enclos_env__$private$num_cores = 1L
					private$run_randomization_iteration(w_des, w_inf, if(use_perms) 1L else NULL, permutations, delta, setup$y_delta, setup$base_template_y, setup$base_template_dead, custom_stat_analysis, setup$lightweight_custom_context)
				}
				# Run warmup TWICE and use the second timing. The first call often pays
				# cold-start penalties (C++ JIT, OS page-cache misses, R bytecode compilation)
				# that inflate the estimate 5–15× vs steady-state cost, causing the guard to
				# wrongly choose parallel for small r values like r = 19.
				system.time(do_warmup_iter())  # First call: discarded (cold-start overhead)
				t_rand_warmup = system.time(do_warmup_iter())[[3]]  # Second call: representative cost
				# For fork cluster: round-trip per task ~10ms, but first call also pays ~300ms
				# cluster-creation cost. For mclapply: per-fork cost ~500ms per worker.
				fork_overhead_estimate = if (isTRUE(private$make_fork_cluster)) 0.01 else 0.5
				cluster_create_overhead = if (isTRUE(private$make_fork_cluster) && is.null(private$fork_cluster)) 0.3 else 0.0

				# Better logic: if total estimated time is much greater than overhead, use all cores.
				# We multiply overhead by actual_rand_cores to account for total fork cost.
				if (t_rand_warmup * r < fork_overhead_estimate * actual_rand_cores * 2.0 + cluster_create_overhead) {
					actual_rand_cores = 1L
				}
			} else if (actual_rand_cores > 1L && !need_thread_objs) {
				# Lightweight stats (yw only) are so fast that fork overhead usually dominates 
				# unless r is very large. Default to serial for lightweight unless num_cores is small.
				if (r < 2000) actual_rand_cores = 1L
			}

			beta_hat_T_diff_ws = unlist(private$par_lapply(1:r, function(idx) {
				set_package_threads(1L)
				suppressWarnings({
					worker_des = if (!is.null(des_template)) des_template$duplicate() else NULL
					worker_inf = if (!is.null(inf_template)) inf_template$duplicate() else NULL
					if (!is.null(worker_inf)) worker_inf$.__enclos_env__$private$num_cores = 1L
					private$run_randomization_iteration(worker_des, worker_inf, if(use_perms) idx else NULL, permutations, delta, setup$y_delta, setup$base_template_y, setup$base_template_dead, custom_stat_analysis, setup$lightweight_custom_context)
				})
			}, n_cores = actual_rand_cores, show_progress = show_progress))

			if (!is.numeric(beta_hat_T_diff_ws)) beta_hat_T_diff_ws = as.numeric(beta_hat_T_diff_ws)
			beta_hat_T_diff_ws
		},

		#' @description
		#' Computes a randomization-based p-value.
		#' @param r		Number of randomization vectors.
		#' @param delta					Null difference.
		#' @param transform_responses	Transformation.
		#' @param na.rm 				Remove NAs.
		#' @param show_progress		Show progress.
		#' @param permutations		Pre-computed permutations.
		#' @return 	Randomization p-value.
		compute_two_sided_pval_for_treatment_effect_rand = function(r = 501, delta = 0, transform_responses = "none", na.rm = TRUE, show_progress = TRUE, permutations = NULL){
			private$assert_design_supports_resampling("Randomization inference")
			assertLogical(na.rm)
			if (private$des_obj_priv_int$response_type == "incidence") stop("Randomization tests are not supported for incidence. Use Zhang method.")
			if (is.null(permutations)) permutations = private$generate_permutations(r)

			cache_key = private$build_randomization_distribution_cache_key(r, delta, transform_responses, permutations)

			if (transform_responses == "none" && is.null(private[["custom_randomization_statistic_function"]]) && !is.null(private$cached_values$t0s_rand) && length(private$cached_values$t0s_rand) >= r) {
				t0s = private$cached_values$t0s_rand[seq_len(r)] + delta
				t = private$compute_treatment_estimate_during_randomization_inference()
				if (length(t) != 1 || !is.finite(t)) return(NA_real_)
				na_t0s = !is.finite(t0s)
				nsim_adj = sum(!na_t0s)
				if (nsim_adj == 0L) return(NA_real_)
				return(min(1, max(2 / nsim_adj, 2 * min(sum(t0s >= t, na.rm = TRUE) / nsim_adj, sum(t0s <= t, na.rm = TRUE) / nsim_adj))))
			}

			if (is.null(private$cached_values$rand_distr_cache)) private$cached_values$rand_distr_cache = list()

			if (!is.null(cache_key) && !is.null(private$cached_values$rand_distr_cache[[cache_key]]) && length(private$cached_values$rand_distr_cache[[cache_key]]) >= r) {
				t0s = private$cached_values$rand_distr_cache[[cache_key]][seq_len(r)]
			} else {
				t0s = self$compute_beta_hat_T_randomization_distr_under_sharp_null(r, delta, transform_responses, show_progress, permutations)
				if (!is.null(cache_key)) private$cached_values$rand_distr_cache[[cache_key]] = t0s
			}

			if (delta == 0 && transform_responses == "none" && is.null(private[["custom_randomization_statistic_function"]])) private$cached_values$t0s_rand = t0s

			t = private$compute_treatment_estimate_during_randomization_inference()
			if (length(t) != 1 || !is.finite(t)) return(NA_real_)
			na_t0s = !is.finite(t0s)
			r_adjusted = sum(!na_t0s)
			if (r_adjusted == 0) return(NA_real_)
			min(1, max(2 / r_adjusted, 2 * min(sum(t0s >= t, na.rm = TRUE) / r_adjusted, sum(t0s <= t, na.rm = TRUE) / r_adjusted)))
		}
	),

	private = list(
		custom_randomization_statistic_function = NULL,

		setup_randomization_template_and_shifts = function(delta, transform_responses){
			bypass_checks = FALSE
			if (transform_responses == "already_transformed"){
				transform_responses = "none"; bypass_checks = TRUE
			}

			template = private$des_obj$duplicate()

			if (transform_responses == "log"){
				if (private$des_obj_priv_int$response_type != "count") template$.__enclos_env__$private$y = log(copy(private$y))
			} else if (transform_responses == "logit"){
				template$.__enclos_env__$private$y = logit(copy(private$y))
			} else if (transform_responses == "log1p"){
				template$.__enclos_env__$private$y = log1p(copy(private$y))
			}

			if (delta != 0 && !bypass_checks){
				if (private$des_obj_priv_int$response_type == "incidence") stop("randomization tests with delta nonzero not supported for incidence")
				if (private$des_obj_priv_int$response_type == "count" && !(transform_responses %in% c("log1p", "already_transformed", "log"))) stop("delta nonzero requires log1p for counts")
				if (private$des_obj_priv_int$response_type == "proportion" && transform_responses != "logit") stop("delta nonzero requires logit for proportions")
				if (private$des_obj_priv_int$response_type == "survival" && transform_responses != "log") stop("delta nonzero requires log for survival")
				
				if (transform_responses == "log" && private$des_obj_priv_int$response_type == "count") {
					template$.__enclos_env__$private$y[private$w == 1] = template$.__enclos_env__$private$y[private$w == 1] * exp(-delta)
				} else {
					template$.__enclos_env__$private$y[private$w == 1] = template$.__enclos_env__$private$y[private$w == 1] - delta
					if (transform_responses == "log" && private$des_obj_priv_int$response_type != "count") template$.__enclos_env__$private$y = exp(template$.__enclos_env__$private$y)
					else if (transform_responses == "logit") template$.__enclos_env__$private$y = inv_logit(template$.__enclos_env__$private$y)
					else if (transform_responses == "log1p") template$.__enclos_env__$private$y = expm1(template$.__enclos_env__$private$y)
				}
			}

			base_template_y = template$.__enclos_env__$private$y
			base_template_dead = template$.__enclos_env__$private$dead
			y_delta = base_template_y
			if (delta != 0) {
				if (transform_responses == "log") y_delta = y_delta * exp(delta)
				else if (transform_responses == "logit") y_delta = inv_logit(logit(y_delta) + delta)
				else if (transform_responses == "log1p") y_delta = (y_delta + 1) * exp(delta) - 1
				else y_delta = y_delta + delta
			}

			lightweight_custom_context = NULL
			if (isTRUE(private$analyze_custom_randomization_statistic()$can_use_lightweight_yw_only)) lightweight_custom_context = private$build_lightweight_custom_randomization_context()

			list(template = template, y_delta = y_delta, base_template_y = base_template_y, base_template_dead = base_template_dead, lightweight_custom_context = lightweight_custom_context)
		},

		run_randomization_iteration = function(thread_des_obj, thread_inf_obj, perm_idx, permutations, delta, y_delta, base_template_y, base_template_dead, custom_stat_analysis, lightweight_custom_context){
			use_perms = !is.null(perm_idx)
			get_perm_data = if (!is.null(permutations$w_mat)) {
				n_avail = ncol(permutations$w_mat)
				function(i) { j = ((i - 1L) %% n_avail) + 1L; list(w = permutations$w_mat[, j], m_vec = if (!is.null(permutations$m_mat)) permutations$m_mat[, j] else NULL) }
			} else function(i) permutations[[i]]

			if (isTRUE(custom_stat_analysis$can_use_lightweight_yw_only) && use_perms) {
				perm_data = get_perm_data(perm_idx); w_sim = perm_data$w; y_sim = base_template_y
				if (delta != 0) y_sim[w_sim == 1] = y_delta[w_sim == 1]
				return(private$evaluate_lightweight_custom_randomization_statistic(lightweight_custom_context, y_sim, w_sim, base_template_dead))
			}

			if (use_perms) {
				perm_data = get_perm_data(perm_idx)
				thread_des_obj$.__enclos_env__$private$w = perm_data$w
				if (private$is_KK && private$object_has_private_method(thread_des_obj, "m")) thread_des_obj$.__enclos_env__$private$m = perm_data$m_vec
				thread_inf_obj$.__enclos_env__$private$w = perm_data$w
				if (private$is_KK && private$object_has_private_method(thread_inf_obj, "m")) thread_inf_obj$.__enclos_env__$private$m = perm_data$m_vec
				thread_inf_obj$.__enclos_env__$private$des_obj_priv_int = thread_des_obj$.__enclos_env__$private
				thread_inf_obj$.__enclos_env__$private$y = thread_des_obj$.__enclos_env__$private$y
				thread_inf_obj$.__enclos_env__$private$dead = thread_des_obj$.__enclos_env__$private$dead
			} else {
				thread_des_obj$.__enclos_env__$private$draw_one_w()
				thread_inf_obj$.__enclos_env__$private$des_obj_priv_int = thread_des_obj$.__enclos_env__$private
				thread_inf_obj$.__enclos_env__$private$y = thread_des_obj$.__enclos_env__$private$y
				thread_inf_obj$.__enclos_env__$private$w = thread_des_obj$.__enclos_env__$private$w
				thread_inf_obj$.__enclos_env__$private$dead = thread_des_obj$.__enclos_env__$private$dead
			}

			thread_inf_obj$.__enclos_env__$private$cached_values = list()
			if (delta != 0) {
				w_sim = thread_inf_obj$.__enclos_env__$private$w; y_sim = thread_inf_obj$.__enclos_env__$private$y
				y_sim[w_sim == 1] = y_delta[w_sim == 1]; thread_inf_obj$.__enclos_env__$private$y = y_sim
			}

			if (private$is_KK && (is.null(private$custom_randomization_statistic_function) || isTRUE(custom_stat_analysis$needs_match_data)) && private$object_has_private_method(thread_inf_obj, "compute_basic_match_data")){
				w_key = if(use_perms) private$stable_signature(perm_data$w) else NULL
				cached_match = if(!is.null(w_key)) private$cached_values$m_cache[[w_key]] else NULL
				if (!is.null(cached_match)) thread_inf_obj$.__enclos_env__$private$m = cached_match
				else {
					if (!use_perms) thread_inf_obj$.__enclos_env__$private$m = thread_des_obj$.__enclos_env__$private$m
					thread_inf_obj$.__enclos_env__$private$compute_basic_match_data()
					if (!is.null(w_key)) private$cached_values$m_cache[[w_key]] = thread_inf_obj$.__enclos_env__$private$m
				}
				if (private$object_has_private_method(thread_inf_obj, "compute_reservoir_and_match_statistics")) thread_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
			}

			estimate = tryCatch(thread_inf_obj$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference(estimate_only = TRUE), error = function(e) NA_real_)
			if (is.list(estimate) && "b" %in% names(estimate)) return(as.numeric(estimate$b[1]))
			as.numeric(estimate)
		},

		build_randomization_distribution_cache_key = function(r, delta, transform_responses, permutations){
			if (is.null(permutations)) return(NULL)
			perm_sig = attr(permutations, "sig")
			if (is.null(perm_sig)) {
				perm_sig = private$permutations_signature(permutations, r)
				attr(permutations, "sig") = perm_sig
			}
			stat_sig = private$custom_randomization_statistic_signature()
			paste("nsim", as.integer(r), "delta", formatC(delta, digits = 17, format = "fg", flag = "#"), "transform", transform_responses, "stat", stat_sig, "perm", perm_sig, sep = "|")
		},

		custom_randomization_statistic_signature = function(){
			if (is.null(private$custom_randomization_statistic_function)) return("default")
			private$stable_signature(list(formals = formals(private$custom_randomization_statistic_function), body = body(private$custom_randomization_statistic_function)))
		},

		permutations_signature = function(permutations, r){
			private$stable_signature(list(r = as.integer(r), permutations = permutations[seq_len(min(length(permutations), r))]))
		},

		next_generated_permutation_signature = function(cache_key){
			counter = private$cached_values$generated_permutation_sig_counter %||% 0L
			counter = counter + 1L; private$cached_values$generated_permutation_sig_counter = counter
			paste("generated", cache_key, counter, sep = ":")
		},

		analyze_custom_randomization_statistic = function(){
			if (!is.null(private$cached_values$custom_stat_analysis)) return(private$cached_values$custom_stat_analysis)
			if (is.null(private$custom_randomization_statistic_function)) {
				analysis = list(can_use_lightweight_yw_only = FALSE, needs_match_data = TRUE)
				private$cached_values$custom_stat_analysis = analysis; return(analysis)
			}
			dollar_paths = private$extract_dollar_paths(body(private$custom_randomization_statistic_function))
			path_strings = vapply(dollar_paths, paste, character(1), collapse = "$")
			allowed_lightweight_paths = c("private$y", "private$w", "private$dead", "private$des_obj_priv_int", "private$des_obj_priv_int$y", "private$des_obj_priv_int$w", "private$des_obj_priv_int$dead")
			references_self = any(vapply(dollar_paths, function(path) length(path) > 0L && identical(path[1], "self"), logical(1)))
			can_use_lightweight_yw_only = !references_self && all(path_strings %in% allowed_lightweight_paths)
			match_tokens = c("match", "reservoir", "pair", "stratum", "strata", "matched", "discordant", "concordant")
			needs_match_data = TRUE
			if (can_use_lightweight_yw_only) needs_match_data = FALSE
			else if (length(path_strings) > 0L) needs_match_data = any(vapply(match_tokens, function(token) any(grepl(token, path_strings, fixed = TRUE)), logical(1)))
			analysis = list(can_use_lightweight_yw_only = can_use_lightweight_yw_only, needs_match_data = needs_match_data)
			private$cached_values$custom_stat_analysis = analysis; return(analysis)
		},

		build_lightweight_custom_randomization_context = function(){
			if (is.null(private$custom_randomization_statistic_function)) return(NULL)
			orig_env = environment(private$custom_randomization_statistic_function)
			eval_env = new.env(parent = orig_env); private_proxy = new.env(parent = emptyenv()); seq_priv_proxy = new.env(parent = emptyenv())
			private_proxy$des_obj_priv_int = seq_priv_proxy; eval_env$private = private_proxy
			custom_fun = private$custom_randomization_statistic_function; environment(custom_fun) = eval_env
			list(fun = custom_fun, private_proxy = private_proxy, des_obj_priv_int_proxy = seq_priv_proxy)
		},

		evaluate_lightweight_custom_randomization_statistic = function(context, y, w, dead = NULL){
			if (is.null(context)) return(NA_real_)
			context$private_proxy$y = y; context$private_proxy$w = w; context$private_proxy$dead = dead
			context$des_obj_priv_int_proxy$y = y; context$des_obj_priv_int_proxy$w = w; context$des_obj_priv_int_proxy$dead = dead
			tryCatch(context$fun(), error = function(e) NA_real_)
		},

		generate_permutations = function(r = 501){
			private$assert_design_supports_resampling("Randomization inference")
			cache_key = as.character(as.integer(r))
			if (is.null(private$cached_values$permutations_cache)) private$cached_values$permutations_cache = list()
			if (!is.null(private$cached_values$permutations_cache[[cache_key]])) return(private$cached_values$permutations_cache[[cache_key]])
			perms = private$des_obj$draw_ws_according_to_design(r)
			if (is.matrix(perms)) perms = list(w_mat = perms, match_indic_mat = NULL)
			attr(perms, "sig") = private$next_generated_permutation_signature(cache_key)
			private$cached_values$permutations_cache[[cache_key]] = perms
			perms
		},

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$custom_randomization_statistic_function)) self$compute_treatment_estimate(estimate_only = estimate_only)
			else private$custom_randomization_statistic_function()
		}
	)
)
