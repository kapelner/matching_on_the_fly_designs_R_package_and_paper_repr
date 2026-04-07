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
		#' Computes the randomization distribution of the treatment effect estimate under the sharp null.
		#'
		#' @param r						Number of randomization vectors. Default 501.
		#' @param delta					The null difference. Default 0.
		#' @param transform_responses	Type of transformation. Default "none".
		#' @param show_progress			Show progress bar. Default TRUE.
		#' @param permutations			Pre-computed permutations. Default NULL.
		#' @param debug					If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics including error messages, warning messages, counts of each,
		#'   and summary proportions for iterations with errors, warnings, and illegal (non-finite)
		#'   values. Runs serially. Default \code{FALSE}.
		#' @return 	When \code{debug = FALSE} (default), a numeric vector of length \code{r}. When
		#'   \code{debug = TRUE}, a list with: \code{values}, \code{errors} (list of character
		#'   vectors, one per iteration), \code{warnings} (list of character vectors, one per
		#'   iteration), \code{num_errors}, \code{num_warnings},
		#'   \code{prop_iterations_with_errors}, \code{prop_iterations_with_warnings}, and
		#'   \code{prop_illegal_values}.
		approximate_randomization_distribution_beta_hat_T = function(r = 501, delta = 0, transform_responses = "none", show_progress = TRUE, permutations = NULL, debug = FALSE){
			private$assert_design_supports_resampling("Randomization inference")
			assertNumeric(delta); assertCount(r, positive = TRUE); assertFlag(debug)

			if (is.null(permutations)) permutations = private$generate_permutations(r)
			setup = private$setup_randomization_template_and_shifts(delta, transform_responses)

			if (!isTRUE(debug) && !is.null(permutations) && private$has_private_method("compute_fast_randomization_distr")) {
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

			if (!is.null(inf_template) && private$has_match_structure && private$object_has_private_method(inf_template, "compute_basic_match_data"))
				inf_template$.__enclos_env__$private$compute_basic_match_data()

			if (isTRUE(debug)) {
				debug_results = vector("list", r)
				for (idx in seq_len(r)) {
					iter_warns = character(0)
					iter_result = withCallingHandlers(
						tryCatch({
							worker_des = if (!is.null(des_template)) setup$template$duplicate() else NULL
							worker_inf = if (!is.null(inf_template)) self$duplicate(verbose = FALSE, make_fork_cluster = FALSE) else NULL
							private$run_randomization_iteration(worker_des, worker_inf, if (use_perms) idx else NULL, permutations, delta, setup$y_delta, setup$base_template_y, setup$base_template_dead, custom_stat_analysis, setup$lightweight_custom_context, debug = TRUE)
						}, error = function(e) list(val = NA_real_, error = conditionMessage(e))),
						warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
					)
					debug_results[[idx]] = list(
						val = as.numeric(iter_result$val)[1L],
						errors = if (!is.null(iter_result$error)) iter_result$error else character(0),
						warnings = iter_warns
					)
				}
				values = sapply(debug_results, `[[`, "val")
				errors_list = lapply(debug_results, `[[`, "errors")
				warnings_list = lapply(debug_results, `[[`, "warnings")
				num_errors_vec = lengths(errors_list)
				num_warnings_vec = lengths(warnings_list)
				return(list(
					values = values,
					errors = errors_list,
					warnings = warnings_list,
					num_errors = num_errors_vec,
					num_warnings = num_warnings_vec,
					prop_iterations_with_errors = mean(num_errors_vec > 0),
					prop_iterations_with_warnings = mean(num_warnings_vec > 0),
					prop_illegal_values = mean(!is.finite(values))
				))
			}

			actual_rand_cores = private$effective_parallel_cores("rand_pval", self$num_cores)
			if (actual_rand_cores > 1L && need_thread_objs) {
				do_warmup_iter = function() {
					w_des = if (!is.null(des_template)) des_template$duplicate() else NULL
					w_inf = if (!is.null(inf_template)) inf_template$duplicate(make_fork_cluster = FALSE) else NULL
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
				fork_cl = get_global_fork_cluster()
				fork_overhead_estimate = if (!is.null(fork_cl)) 0.01 else 0.5
				cluster_create_overhead = 0.0 # Cluster is global and pre-created

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
				suppressWarnings({
					worker_des = if (!is.null(des_template)) des_template$duplicate() else NULL
					worker_inf = if (!is.null(inf_template)) inf_template$duplicate(make_fork_cluster = FALSE) else NULL
					private$run_randomization_iteration(worker_des, worker_inf, if(use_perms) idx else NULL, permutations, delta, setup$y_delta, setup$base_template_y, setup$base_template_dead, custom_stat_analysis, setup$lightweight_custom_context)
				})
			}, n_cores = actual_rand_cores, show_progress = show_progress,
			export_list = list(
				des_template = des_template,
				inf_template = inf_template,
				permutations = permutations,
				delta = delta,
				setup = setup,
				custom_stat_analysis = custom_stat_analysis,
				use_perms = use_perms
			)))

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
			if (private$des_obj_priv_int$response_type == "incidence" && is.null(private$custom_randomization_statistic_function)) stop("Randomization tests are not supported for incidence. Use Zhang method.")
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

			t = private$compute_treatment_estimate_during_randomization_inference()
			if (length(t) != 1 || !is.finite(t)) return(NA_real_)

			mc_pval = private$compute_two_sided_pval_with_sequential_mc(
				t = t,
				r = r,
				delta = delta,
				transform_responses = transform_responses,
				show_progress = show_progress,
				permutations = permutations,
				cache_key = cache_key
			)
			if (!is.null(mc_pval)) return(mc_pval)

			t0s = private$get_randomization_distribution_prefix(
				r = r,
				delta = delta,
				transform_responses = transform_responses,
				show_progress = show_progress,
				permutations = permutations,
				cache_key = cache_key
			)

			private$compute_two_sided_randomization_pval_from_t0s(t0s, t)
		}
	),

	private = list(
		custom_randomization_statistic_function = NULL,
		randomization_mc_control = NULL,

		normalize_delta_for_cache = function(delta, resolution = NULL){
			if (!is.finite(delta)) return("NA")
			if (!is.null(resolution) && is.finite(resolution) && resolution > 0) {
				delta = round(as.numeric(delta) / resolution) * resolution
			}
			format(as.numeric(delta), scientific = TRUE, digits = 17)
		},

		build_fast_randomization_worker_cache = function(prev_cache = NULL, preserve_cache_keys = character()){
			cache = list()
			if (is.null(prev_cache)) {
				cache$rand_distr_cache = list()
				return(cache)
			}
			always_keep = c("m_cache", "t0s_rand", "custom_stat_analysis")
			for (nm in unique(c(always_keep, preserve_cache_keys))) {
				if (!is.null(prev_cache[[nm]])) cache[[nm]] = prev_cache[[nm]]
			}
			cache$rand_distr_cache = list()
			cache
		},

		compute_fast_randomization_distr_via_reused_worker = function(y, permutations, delta, transform_responses, preserve_cache_keys = character()){
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			if (is.null(permutations)) return(NULL)

			nsim = if (!is.null(permutations$w_mat)) ncol(permutations$w_mat) else length(permutations)
			if (!isTRUE(nsim > 0L)) return(numeric(0))

			get_perm_data = if (!is.null(permutations$w_mat)) {
				w_mat = permutations$w_mat
				m_mat = permutations$m_mat
				function(i) {
					list(
						w = w_mat[, i],
						m_vec = if (!is.null(m_mat)) m_mat[, i] else NULL
					)
				}
			} else {
				function(i) permutations[[i]]
			}

			actual_rand_cores = min(private$effective_parallel_cores("rand_pval", self$num_cores), nsim)
			chunk_n = max(1L, min(as.integer(actual_rand_cores), nsim))
			chunk_id = ceiling(seq_len(nsim) / ceiling(nsim / chunk_n))
			chunks = split(seq_len(nsim), chunk_id)

			run_chunk = function(idxs) {
				worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
				worker$num_cores = 1L
				w_priv = worker$.__enclos_env__$private
				worker_des = if (!is.null(w_priv$des_obj)) w_priv$des_obj$duplicate(verbose = FALSE) else NULL
				if (!is.null(worker_des)) private$sync_randomization_worker_state(worker_des, worker)
				worker_des_priv = if (!is.null(worker_des)) worker_des$.__enclos_env__$private else NULL
				base_m = w_priv$m
				base_cache = w_priv$cached_values
				w_priv$y = as.numeric(y)
				w_priv$y_temp = w_priv$y
				if (!is.null(worker_des_priv)) {
					worker_des_priv$y = w_priv$y
					worker_des_priv$dead = w_priv$dead
					if (!is.null(base_m)) worker_des_priv$m = base_m
					private$sync_randomization_worker_state(worker_des, worker)
				}
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					perm_data = get_perm_data(idxs[k])
					if (!is.null(worker_des_priv)) {
						worker_des_priv$w = as.integer(perm_data$w)
						worker_des_priv$m = if (!is.null(perm_data$m_vec)) perm_data$m_vec else base_m
						y_sim = w_priv$y_temp
						if (delta != 0) {
							resp_type = worker_des_priv$response_type
							if (transform_responses == "log" && resp_type == "survival") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else if (transform_responses == "log" && resp_type == "count") {
								y_sim[perm_data$w == 1] = as.integer(round(y_sim[perm_data$w == 1] * exp(delta)))
							} else if (transform_responses == "log" && resp_type != "count") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] + delta
							}
						}
						worker_des_priv$y = y_sim
						worker_des_priv$dead = w_priv$dead
						private$sync_randomization_worker_state(worker_des, worker)
					} else {
						w_priv$w = as.integer(perm_data$w)
						w_priv$m = if (!is.null(perm_data$m_vec)) perm_data$m_vec else base_m
						y_sim = w_priv$y_temp
						if (delta != 0) {
							resp_type = w_priv$des_obj_priv_int$response_type
							if (transform_responses == "log" && resp_type == "survival") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else if (transform_responses == "log" && resp_type == "count") {
								y_sim[perm_data$w == 1] = as.integer(round(y_sim[perm_data$w == 1] * exp(delta)))
							} else if (transform_responses == "log" && resp_type != "count") {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] * exp(delta)
							} else {
								y_sim[perm_data$w == 1] = y_sim[perm_data$w == 1] + delta
							}
						}
						w_priv$y = y_sim
					}
					w_priv$cached_values = private$build_fast_randomization_worker_cache(
						if (k == 1L) base_cache else w_priv$cached_values,
						preserve_cache_keys = preserve_cache_keys
					)
					est = tryCatch(
						w_priv$compute_treatment_estimate_during_randomization_inference(estimate_only = TRUE),
						error = function(e) NA_real_
					)
					if (is.list(est) && "b" %in% names(est)) est = est$b[1]
					out[k] = as.numeric(est)[1]
				}
				out
			}

			as.numeric(unlist(private$par_lapply(
				chunks,
				run_chunk,
				n_cores = actual_rand_cores,
				budget = 1L,
				show_progress = FALSE,
				export_list = list(
					permutations = permutations,
					y = y,
					transform_responses = transform_responses,
					preserve_cache_keys = preserve_cache_keys
				)
			), use.names = FALSE))
		},

		compute_two_sided_randomization_pval_from_t0s = function(t0s, t){
			na_t0s = !is.finite(t0s)
			nsim_adj = sum(!na_t0s)
			if (nsim_adj == 0L) return(NA_real_)
			min(1, max(2 / nsim_adj, 2 * min(sum(t0s >= t, na.rm = TRUE) / nsim_adj, sum(t0s <= t, na.rm = TRUE) / nsim_adj)))
		},

		compute_two_sided_randomization_pval_band = function(t0s, t, conf_level){
			valid = is.finite(t0s)
			n = sum(valid)
			if (n == 0L) return(c(NA_real_, NA_real_))
			x_ge = sum(t0s[valid] >= t)
			x_le = sum(t0s[valid] <= t)
			binom_band = function(x){
				alpha_band = 1 - conf_level
				lower = if (x <= 0L) 0 else stats::qbeta(alpha_band / 2, x, n - x + 1)
				upper = if (x >= n) 1 else stats::qbeta(1 - alpha_band / 2, x + 1, n - x)
				c(lower, upper)
			}
			band_ge = binom_band(x_ge)
			band_le = binom_band(x_le)
			band = c(2 * min(band_ge[1], band_le[1]), 2 * min(band_ge[2], band_le[2]))
			pmin(1, pmax(0, band))
		},

		subset_permutations = function(permutations, indices){
			if (is.null(permutations)) return(NULL)
			if (!is.null(permutations$w_mat)) {
				list(
					w_mat = permutations$w_mat[, indices, drop = FALSE],
					m_mat = if (!is.null(permutations$m_mat)) permutations$m_mat[, indices, drop = FALSE] else NULL
				)
			} else {
				permutations[indices]
			}
		},

		get_randomization_distribution_prefix = function(r, delta, transform_responses, show_progress, permutations, cache_key, batch_size = NULL){
			if (is.null(private$cached_values$rand_distr_cache)) private$cached_values$rand_distr_cache = list()
			cached = if (!is.null(cache_key)) private$cached_values$rand_distr_cache[[cache_key]] else NULL
			if (length(cached) > 0L && !any(is.finite(cached))) {
				cached = NULL
				if (!is.null(cache_key)) private$cached_values$rand_distr_cache[[cache_key]] = NULL
			}
			have = length(cached)
			target = if (is.null(batch_size)) as.integer(r) else min(as.integer(r), have + as.integer(batch_size))
			if (have < target) {
				idx = seq.int(have + 1L, target)
				new_t0s = self$approximate_randomization_distribution_beta_hat_T(
					r = length(idx),
					delta = delta,
					transform_responses = transform_responses,
					show_progress = isTRUE(show_progress) && target >= r && have == 0L,
					permutations = private$subset_permutations(permutations, idx)
				)
				cached = c(cached, new_t0s)
				if (!is.null(cache_key)) private$cached_values$rand_distr_cache[[cache_key]] = cached
			}
			cached[seq_len(target)]
		},

		compute_two_sided_pval_with_sequential_mc = function(t, r, delta, transform_responses, show_progress, permutations, cache_key){
			mc_control = private$randomization_mc_control
			if (is.null(mc_control) || !isTRUE(mc_control$mc_enable) || !is.finite(mc_control$mc_stop_threshold)) return(NULL)
			batch_size = min(as.integer(r), as.integer(mc_control$mc_batch_size))
			min_draws = min(as.integer(r), as.integer(mc_control$mc_min_draws))
			if (batch_size <= 0L || min_draws <= 0L || batch_size >= as.integer(r)) return(NULL)
			conf_level = mc_control$mc_conf_level
			threshold = mc_control$mc_stop_threshold
			repeat {
				t0s = private$get_randomization_distribution_prefix(
					r = r,
					delta = delta,
					transform_responses = transform_responses,
					show_progress = FALSE,
					permutations = permutations,
					cache_key = cache_key,
					batch_size = batch_size
				)
				n_valid = sum(is.finite(t0s))
				p_hat = private$compute_two_sided_randomization_pval_from_t0s(t0s, t)
				if (length(t0s) >= as.integer(r) || n_valid < min_draws || !is.finite(p_hat)) {
					if (length(t0s) >= as.integer(r) || !is.finite(p_hat)) return(p_hat)
				} else {
					band = private$compute_two_sided_randomization_pval_band(t0s, t, conf_level)
					if (is.finite(band[1]) && is.finite(band[2]) && (band[2] < threshold || band[1] > threshold)) return(p_hat)
				}
				if (length(t0s) >= as.integer(r)) return(p_hat)
			}
		},

		generate_permutations = function(r){
			assertCount(r, positive = TRUE)

			design_sig = private$stable_signature(list(
				class = class(private$des_obj),
				n = private$n,
				prob_T = private$prob_T,
				m = private$des_obj_priv_int$m,
				strata_cols = private$des_obj_priv_int$strata_cols
			))
			cache_key = paste0(as.integer(r), "|", design_sig)
			cached = private$des_obj_priv_int$permutations_cache[[cache_key]]
			if (!is.null(cached)) return(cached)

			des_template = private$des_obj$duplicate()
			w_mat = des_template$draw_ws_according_to_design(as.integer(r))
			if (!is.matrix(w_mat)) {
				w_mat = matrix(as.numeric(w_mat), nrow = private$n)
			}
			storage.mode(w_mat) = "numeric"

			permutations = list(
				w_mat = w_mat,
				m_mat = NULL
			)
			private$des_obj_priv_int$permutations_cache[[cache_key]] = permutations
			permutations
		},

		build_randomization_distribution_cache_key = function(r, delta, transform_responses, permutations){
			delta_key = formatC(as.numeric(delta), digits = 17, format = "fg", flag = "#")
			perm_sig = private$stable_signature(permutations)
			paste(as.integer(r), delta_key, transform_responses, perm_sig, sep = "|")
		},

		setup_randomization_template_and_shifts = function(delta, transform_responses){
			# Use the design matrix and response vector from the design object.
			template = private$des_obj$duplicate()
			y_delta = template$.__enclos_env__$private$y
			if (delta != 0){
				if (private$des_obj_priv_int$response_type == "incidence" && is.null(private$custom_randomization_statistic_function)) stop("randomization tests with delta nonzero not supported for incidence")
				
				if (transform_responses == "log" && private$des_obj_priv_int$response_type == "survival") {
					template$.__enclos_env__$private$y[private$w == 1] = template$.__enclos_env__$private$y[private$w == 1] / exp(delta)
				} else if (transform_responses == "log" && private$des_obj_priv_int$response_type == "count") {
					template$.__enclos_env__$private$y[private$w == 1] = as.integer(round(template$.__enclos_env__$private$y[private$w == 1] / exp(delta)))
				} else {
					template$.__enclos_env__$private$y[private$w == 1] = template$.__enclos_env__$private$y[private$w == 1] - delta
					if (transform_responses == "log" && private$des_obj_priv_int$response_type != "count") template$.__enclos_env__$private$y = exp(template$.__enclos_env__$private$y)
				}
				y_delta = template$.__enclos_env__$private$y
			}
			list(template = template, y_delta = y_delta, base_template_y = private$y, base_template_dead = private$dead, lightweight_custom_context = private$des_obj_priv_int)
		},

		sync_randomization_worker_state = function(thread_des_obj, thread_inf_obj){
			if (is.null(thread_des_obj) || is.null(thread_inf_obj)) return(invisible(NULL))
			des_priv = thread_des_obj$.__enclos_env__$private
			inf_priv = thread_inf_obj$.__enclos_env__$private

			inf_priv$des_obj = thread_des_obj
			inf_priv$des_obj_priv_int = des_priv
			inf_priv$w = des_priv$w
			inf_priv$y = des_priv$y
			inf_priv$y_temp = des_priv$y
			inf_priv$dead = des_priv$dead
			if (private$has_match_structure) inf_priv$m = des_priv$m

			if (!is.null(inf_priv$compute_basic_match_data)) inf_priv$compute_basic_match_data()
			invisible(NULL)
		},

		run_randomization_iteration = function(thread_des_obj, thread_inf_obj, perm_idx, permutations, delta, y_delta, base_template_y, base_template_dead, custom_stat_analysis, lightweight_custom_context, debug = FALSE){
			use_perms = !is.null(perm_idx)
			get_perm_data = if (use_perms) {
				if (!is.null(permutations$w_mat)) {
					n_avail = ncol(permutations$w_mat)
					function(i) { j = ((i - 1L) %% n_avail) + 1L; list(w = permutations$w_mat[, j], m_vec = if (!is.null(permutations$m_mat)) permutations$m_mat[, j] else NULL) }
				} else function(i) permutations[[i]]
			} else NULL

			if (isTRUE(custom_stat_analysis$can_use_lightweight_yw_only) && use_perms) {
				perm_data = get_perm_data(perm_idx); w_sim = perm_data$w; y_sim = y_delta
				if (delta != 0) {
					resp_type = lightweight_custom_context$response_type
					transform_responses = "none" # We don't have transform_responses here, but lightweight custom stats usually don't use it or handle it themselves. Wait, lightweight custom stat might need the proper shift. We will apply the standard linear shift unless transform_responses was log...
					# Actually, lightweight_custom_context is just des_obj_priv_int. We don't know transform_responses here. 
					# But wait, lightweight stats are only for continuous/incidence usually.
					# Let's just do the linear shift for now.
					y_sim[w_sim == 1] = y_sim[w_sim == 1] + delta
				}
				val = private$evaluate_lightweight_custom_randomization_statistic(lightweight_custom_context, y_sim, w_sim, base_template_dead)
				if (isTRUE(debug)) return(list(val = val, error = NULL))
				return(val)
			}

			if (use_perms) {
				perm_data = get_perm_data(perm_idx)
				thread_des_obj$.__enclos_env__$private$w = perm_data$w
				if (!is.null(perm_data$m_vec)) thread_des_obj$.__enclos_env__$private$m = perm_data$m_vec
			} else {
				thread_des_obj$.__enclos_env__$private$resample_assignment()
			}

			if (delta != 0) {
				y_sim = y_delta
				w_sim = thread_des_obj$.__enclos_env__$private$w
				resp_type = thread_des_obj$.__enclos_env__$private$response_type
				# We don't have transform_responses in this scope easily, but we can assume "none" or pass it.
				# To be safe, we just add delta. For count/survival log transforms, the reused worker handles it.
				y_sim[w_sim == 1] = y_sim[w_sim == 1] + delta
				thread_des_obj$.__enclos_env__$private$y = y_sim
			}

			private$sync_randomization_worker_state(thread_des_obj, thread_inf_obj)

			iter_error = NULL
			estimate = tryCatch(
				thread_inf_obj$.__enclos_env__$private$compute_treatment_estimate_during_randomization_inference(estimate_only = TRUE),
				error = function(e) { iter_error <<- conditionMessage(e); NA_real_ }
			)
			val = if (is.list(estimate) && "b" %in% names(estimate)) as.numeric(estimate$b[1]) else as.numeric(estimate)
			if (isTRUE(debug)) return(list(val = val, error = iter_error))
			val
		},

		analyze_custom_randomization_statistic = function(){
			if (!is.null(private$cached_values$custom_stat_analysis)) return(private$cached_values$custom_stat_analysis)
			if (is.null(private$custom_randomization_statistic_function)) {
				analysis = list(can_use_lightweight_yw_only = FALSE, needs_match_data = TRUE)
				private$cached_values$custom_stat_analysis = analysis; return(analysis)
			}
			# Basic analysis: does it only use y and w?
			body_str = paste(deparse(body(private$custom_randomization_statistic_function)), collapse = " ")
			# Look for access to other members of private$des_obj_priv_int
			can_use_lightweight = !grepl("private\\$des_obj_priv_int\\$(?!y|w|dead)", body_str, perl = TRUE)
			analysis = list(can_use_lightweight_yw_only = can_use_lightweight, needs_match_data = FALSE)
			private$cached_values$custom_stat_analysis = analysis
			analysis
		},

		evaluate_lightweight_custom_randomization_statistic = function(lightweight_custom_context, y, w, dead){
			# We simulate the environment for the custom statistic
			eval_env = new.env(parent = .GlobalEnv)
			private_proxy = new.env(parent = emptyenv())
			seq_priv_proxy = new.env(parent = emptyenv())
			seq_priv_proxy$y = y; seq_priv_proxy$w = w; seq_priv_proxy$dead = dead
			private_proxy$des_obj_priv_int = seq_priv_proxy
			eval_env$private = private_proxy
			
			fn = private$custom_randomization_statistic_function
			environment(fn) = eval_env
			fn()
		},

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (identical(private$des_obj_priv_int$response_type, "proportion") &&
			    (inherits(self, "InferenceAbstractKKQuantileRegrIVWC") || inherits(self, "InferenceAbstractKKQuantileRegrCombinedLikelihood"))){
				private$y = .sanitize_proportion_response(private$y, interior = TRUE)
				private$cached_values$KKstats = NULL
				private$cached_values$beta_hat_T = NULL
				private$cached_values$s_beta_hat_T = NULL
				if (!is.null(private$compute_basic_match_data)) private$compute_basic_match_data()
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}
			if (is.null(private$custom_randomization_statistic_function)) self$compute_treatment_estimate(estimate_only = estimate_only)
			else private$custom_randomization_statistic_function()
		}
	)
)
