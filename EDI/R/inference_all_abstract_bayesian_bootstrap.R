#' Bayesian Bootstrap-capable Inference
#'
#' Abstract class for Dirichlet-weight Bayesian bootstrap inference layered on
#' top of the existing nonparametric bootstrap infrastructure.
#'
#' @keywords internal
InferenceBayesianBootstrap = R6::R6Class("InferenceBayesianBootstrap",
	lock_objects = FALSE,
	inherit = InferenceNonParamBootstrap,
	public = list(
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			stop(class(self)[1], " must implement weighted bootstrap estimation.")
		},
		approximate_bayesian_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Bayesian bootstrap inference")
				private$assert_valid_bootstrap_type(bootstrap_type)
				assertCount(B, positive = TRUE)
				assertFlag(debug)
			}
			cache_key = private$bayesian_bootstrap_cache_key(B = B, bootstrap_type = bootstrap_type)
			if (!isTRUE(debug) && !is.null(private$cached_values$bayes_boot_distr_cache[[cache_key]])) {
				return(private$cached_values$bayes_boot_distr_cache[[cache_key]])
			}
			inf_template = self$duplicate()
			run_one_iter = function(worker_inf) {
				draw = private$bayesian_bootstrap_sample_weights(bootstrap_type = bootstrap_type)
				worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
				as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
					subject_or_block_weights = draw$subject_or_block_weights,
					estimate_only = TRUE
				))[1L]
			}
			if (isTRUE(debug)) {
				run_debug_iter = function(worker_inf = NULL, worker_state = NULL) {
					iter_warns = character(0)
					iter_errs = character(0)
					iter_val = withCallingHandlers(
						tryCatch({
							if (!is.null(worker_state)) {
								draw = private$bayesian_bootstrap_sample_weights(bootstrap_type = bootstrap_type)
								private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
								private$compute_bayesian_bootstrap_worker_estimate(worker_state)
							} else {
								run_one_iter(worker_inf)
							}
						}, error = function(e) { iter_errs <<- c(iter_errs, conditionMessage(e)); NA_real_ }),
						warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
					)
					list(val = as.numeric(iter_val)[1L], errors = iter_errs, warnings = iter_warns)
				}
				actual_debug_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
				chunk_n = max(1L, min(as.integer(actual_debug_cores), as.integer(B)))
				chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
				chunks = split(seq_len(B), chunk_id)
				run_debug_chunk = if (isTRUE(private$use_reusable_bootstrap_worker())) {
					function(idxs) {
						worker_state = private$create_bootstrap_worker_state()
						lapply(idxs, function(idx) run_debug_iter(worker_state = worker_state))
					}
				} else {
					function(idxs) {
						lapply(idxs, function(idx) {
							worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
							run_debug_iter(worker_inf = worker_inf)
						})
					}
				}
				debug_results = if (actual_debug_cores <= 1L) {
					run_debug_chunk(seq_len(B))
				} else {
					unlist(private$par_lapply(
						chunks,
						run_debug_chunk,
						n_cores = actual_debug_cores,
						budget = 1L,
						show_progress = show_progress
					), recursive = FALSE, use.names = FALSE)
				}
				debug_results = Filter(function(x) is.list(x) && !is.null(x$val), debug_results)
				values = vapply(debug_results, function(x) as.numeric(x$val)[1L], numeric(1))
				errors_list = lapply(debug_results, `[[`, "errors")
				warnings_list = lapply(debug_results, `[[`, "warnings")
				num_errors_vec = lengths(errors_list)
				num_warnings_vec = lengths(warnings_list)
				if (is.null(private$cached_values$bayes_boot_distr_cache)) private$cached_values$bayes_boot_distr_cache = list()
				private$cached_values$bayes_boot_distr_cache[[cache_key]] = values
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
			actual_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
			if (actual_cores > 1L) {
				do_warmup_iter = if (isTRUE(private$use_reusable_bootstrap_worker())) {
					function() {
						worker_state = private$create_bootstrap_worker_state()
						draw = private$bayesian_bootstrap_sample_weights(bootstrap_type = bootstrap_type)
						private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
						tryCatch(private$compute_bayesian_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
					}
				} else {
					function() {
						worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
						tryCatch(run_one_iter(worker_inf), error = function(e) NA_real_)
					}
				}
				system.time(do_warmup_iter())
				t_boot_warmup = system.time(do_warmup_iter())[[3]]
				fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.5
				if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores)) actual_cores = 1L
			}
			boot_distr = if (isTRUE(private$use_reusable_bootstrap_worker())) {
				private$compute_bayesian_bootstrap_distribution_with_reused_workers(
					B = B,
					actual_cores = actual_cores,
					show_progress = show_progress,
					bootstrap_type = bootstrap_type
				)
			} else {
				unlist(private$par_lapply(
					1:B,
					function(idx) {
						worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
						tryCatch(run_one_iter(worker_inf), error = function(e) NA_real_)
					},
					n_cores = actual_cores,
					show_progress = show_progress,
					export_list = list(inf_template = inf_template, run_one_iter = run_one_iter)
				))
			}
			boot_distr = as.numeric(boot_distr)
			if (is.null(private$cached_values$bayes_boot_distr_cache)) private$cached_values$bayes_boot_distr_cache = list()
			private$cached_values$bayes_boot_distr_cache[[cache_key]] = boot_distr
			boot_distr
		},
		compute_bayesian_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = FALSE, show_progress = TRUE, min_number_usable_samples = 5L, bootstrap_type = NULL){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertFlag(na.rm)
			}
			type = tolower(type %||% "percentile")
			if (should_run_asserts()) {
				assertChoice(type, c("percentile", "symmetric", "wald"))
			}
			boot_distr = self$approximate_bayesian_bootstrap_distribution_beta_hat_T(
				B = B,
				show_progress = show_progress,
				bootstrap_type = bootstrap_type
			)
			if (isTRUE(na.rm)) boot_distr = boot_distr[is.finite(boot_distr)]
			else if (any(!is.finite(boot_distr))) return(NA_real_)
			if (length(boot_distr) < as.integer(min_number_usable_samples)) return(NA_real_)
			est = as.numeric(self$compute_estimate())[1L]
			if (!is.finite(est)) return(NA_real_)
			if (type == "percentile") {
				boot_null = boot_distr - mean(boot_distr) + delta
				n_bs = length(boot_null)
				return(min(1, max(2 / n_bs, 2 * min(
					sum(boot_null >= est) / n_bs,
					sum(boot_null <= est) / n_bs
				))))
			}
			if (type == "symmetric") {
				D_obs = abs(est - delta)
				D_boot = abs(boot_distr - mean(boot_distr))
				n_bs = length(D_boot)
				return(min(1, max(1 / n_bs, mean(D_boot >= D_obs))))
			}
			se_boot = stats::sd(boot_distr)
			if (!is.finite(se_boot) || se_boot <= 0) return(NA_real_)
			2 * stats::pnorm(-abs((est - delta) / se_boot))
		},
		compute_bayesian_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, bootstrap_type = NULL){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
			}
			type = tolower(type %||% "percentile")
			if (should_run_asserts()) {
				assertChoice(type, c("percentile", "basic", "wald"))
			}
			est = as.numeric(self$compute_estimate())[1L]
			ci = c(NA_real_, NA_real_)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			if (!is.finite(est)) return(ci)
			boot_distr = self$approximate_bayesian_bootstrap_distribution_beta_hat_T(
				B = B,
				show_progress = show_progress,
				bootstrap_type = bootstrap_type
			)
			boot_distr = boot_distr[is.finite(boot_distr)]
			if (length(boot_distr) < as.integer(min_number_usable_samples)) return(ci)
			if (type == "wald") {
				se_boot = stats::sd(boot_distr)
				if (!is.finite(se_boot) || se_boot <= 0) return(ci)
				z = stats::qnorm(1 - alpha / 2)
				ci[] = c(est - z * se_boot, est + z * se_boot)
				return(ci)
			}
			ci[] = private$ci_from_boot_distribution(boot_distr, alpha, type, est = est)
			ci
		}
	),
	private = list(
		current_bayesian_bootstrap_context = NULL,
		bayesian_bootstrap_cache_key = function(B, bootstrap_type = NULL){
			paste0(as.integer(B), "::", bootstrap_type %||% "default")
		},
		build_bayesian_bootstrap_context = function(bootstrap_type = NULL){
			n = private$des_obj$get_n()
			design_obj = private$des_obj
			is_matching_design = is(design_obj, "DesignMatching") &&
				isTRUE(tryCatch(design_obj$is_matching_design(), error = function(e) FALSE))
			is_blocking_design = is(design_obj, "DesignBlocking") &&
				isTRUE(tryCatch(design_obj$is_blocking_design(), error = function(e) FALSE))
			if (is_matching_design) {
				private$des_obj_priv_int$ensure_matching_structure_computed()
				cluster_ids = as.integer(design_obj$get_matching_cluster_ids(private$m))
				row_to_unit = match(cluster_ids, unique(cluster_ids))
				unit_group_id = rep(1L, max(row_to_unit))
				return(list(
					row_to_unit = as.integer(row_to_unit),
					unit_group_id = as.integer(unit_group_id),
					n_units = length(unit_group_id)
				))
			}
			if (is_blocking_design) {
				block_ids = as.integer(design_obj$get_block_ids())
				if (identical(bootstrap_type, "resample_blocks")) {
					row_to_unit = match(block_ids, unique(block_ids))
					unit_group_id = rep(1L, max(row_to_unit))
					return(list(
						row_to_unit = as.integer(row_to_unit),
						unit_group_id = as.integer(unit_group_id),
						n_units = length(unit_group_id)
					))
				}
				if (is(design_obj, "DesignFixedBlockedCluster")) {
					cluster_ids = as.character(private$des_obj_priv_int$Xraw[seq_len(n), ][[private$des_obj_priv_int$cluster_col]])
					unit_keys = paste(block_ids, cluster_ids, sep = "::")
					row_to_unit = match(unit_keys, unique(unit_keys))
					unit_group_id = block_ids[match(unique(unit_keys), unit_keys)]
					unit_group_id = match(unit_group_id, unique(unit_group_id))
					return(list(
						row_to_unit = as.integer(row_to_unit),
						unit_group_id = as.integer(unit_group_id),
						n_units = length(unit_group_id)
					))
				}
				row_to_unit = seq_len(n)
				unit_group_id = match(block_ids, unique(block_ids))
				return(list(
					row_to_unit = as.integer(row_to_unit),
					unit_group_id = as.integer(unit_group_id),
					n_units = n
				))
			}
			list(
				row_to_unit = seq_len(n),
				unit_group_id = rep(1L, n),
				n_units = n
			)
		},
		bayesian_bootstrap_sample_weights = function(bootstrap_type = NULL){
			ctx = private$build_bayesian_bootstrap_context(bootstrap_type = bootstrap_type)
			subject_or_block_weights = numeric(ctx$n_units)
			for (group_id in unique(ctx$unit_group_id)) {
				idx = which(ctx$unit_group_id == group_id)
				draw = stats::rgamma(length(idx), shape = 1, rate = 1)
				subject_or_block_weights[idx] = draw / sum(draw) * length(idx)
			}
			list(
				subject_or_block_weights = as.numeric(subject_or_block_weights),
				context = ctx
			)
		},
		expand_subject_or_block_weights_to_row_weights = function(subject_or_block_weights){
			ctx = private$current_bayesian_bootstrap_context
			if (is.null(ctx)) {
				stop("No Bayesian-bootstrap context is installed on this inference object.", call. = FALSE)
			}
			if (should_run_asserts()) {
				assertNumeric(subject_or_block_weights, len = ctx$n_units, lower = 0, any.missing = FALSE)
			}
			as.numeric(subject_or_block_weights[ctx$row_to_unit])
		},
		load_bayesian_bootstrap_weights_into_worker = function(worker_state, draw){
			worker_state$current_subject_or_block_weights = as.numeric(draw$subject_or_block_weights)
			worker_state$current_bayesian_bootstrap_context = draw$context
			worker_state$worker$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
		},
		compute_bayesian_bootstrap_worker_estimate = function(worker_state){
			theta = as.numeric(worker_state$worker$compute_estimate_with_bootstrap_weights(
				subject_or_block_weights = worker_state$current_subject_or_block_weights,
				estimate_only = TRUE
			))[1L]
			if (is.function(worker_state$worker$is_nonestimable) &&
			    isTRUE(worker_state$worker$is_nonestimable("estimate"))) {
				return(NA_real_)
			}
			theta
		},
		compute_bayesian_bootstrap_distribution_with_reused_workers = function(B, actual_cores, show_progress = FALSE, bootstrap_type = NULL){
			chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
			chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
			chunks = split(seq_len(B), chunk_id)
			run_chunk = function(idxs) {
				worker_state = private$create_bootstrap_worker_state()
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					draw = private$bayesian_bootstrap_sample_weights(bootstrap_type = bootstrap_type)
					out[k] = tryCatch({
						private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
						private$compute_bayesian_bootstrap_worker_estimate(worker_state)
					}, error = function(e) NA_real_)
				}
				out
			}
			if (actual_cores <= 1L) {
				return(as.numeric(run_chunk(seq_len(B))))
			}
			as.numeric(unlist(private$par_lapply(
				chunks,
				run_chunk,
				n_cores = actual_cores,
				budget = 1L,
				show_progress = show_progress
			), use.names = FALSE))
		}
	)
)
