#' Mixin for KK Matching-on-the-Fly Pass-Through Inference
#'
#' A Pattern-1 mixin (plain list with \code{$public} and \code{$private} slots)
#' providing KK design validation, bootstrap distribution logic, and match-data
#' helpers. Splice into a daughter class via
#' \code{public = c(InferenceMixinKKPassThrough$public, list(...))} and
#' \code{private = c(InferenceMixinKKPassThrough$private, list(...))}.
#'
#' Capability flag: \code{private$kk_passthrough == TRUE}.
#'
#' @keywords internal
InferenceMixinKKPassThrough = list(
	public = list(
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect
		#'
		#' @param B  					Number of bootstrap samples. The default is 501.
		#'
		#' @return A vector of length \code{B} with the bootstrap values of the estimates of the
		#'   treatment effect
		#'
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics. Default \code{FALSE}.
		#' @param bootstrap_type Optional bootstrap-resampling scheme. Legal public values are:
		#'   \describe{
		#'     \item{\code{NULL}}{Use the design's default row-resampling bootstrap.}
		#'     \item{\code{"within_blocks"}}{Only legal for blocking-style designs.}
		#'     \item{\code{"resample_blocks"}}{Only legal for the same blocking-style designs.}
		#'   }
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
				if (should_run_asserts()) {
					private$assert_valid_bootstrap_type(bootstrap_type)
				}
				if (!private$has_match_structure){
					super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug = debug, bootstrap_type = bootstrap_type)
				} else {
					if (should_run_asserts()) {
						assertCount(B, positive = TRUE)
					}
					if (is.null(private$cached_values$KKstats)){
						private$compute_basic_match_data()
					}
					n = private$n
					y = private$y
					dead = private$dead
					w = private$w
					X = private$get_X()
					# Let Design initialise and own the bootstrap pair structure
					des_priv = private$des_obj_priv_int
					des_priv$init_matching_bootstrap_structure()
					n_reservoir = des_priv$boot_n_reservoir
					m            = nrow(des_priv$boot_pair_rows)
					# For the C++ fast-path we still need i_reservoir / m_vec in Inference scope
					i_reservoir  = des_priv$boot_i_reservoir
					m_vec = private$m
					if (is.null(m_vec)) m_vec = rep(0L, n)
					m_vec = as.integer(m_vec)
					m_vec[is.na(m_vec)] = 0L
					# Check if subclass provides a C++ OpenMP dispatcher to bypass the slow R loop
					if (!isTRUE(debug) && private$has_private_method("compute_fast_bootstrap_distr")) {
						fast_distr = private$compute_fast_bootstrap_distr(B, i_reservoir, n_reservoir, m, y, w, m_vec)
						if (!is.null(fast_distr)) {
							return(fast_distr)
						}
					}
					kk_boot_context = private$create_kk_bootstrap_context(
						y = y, dead = dead, w = w, X = X,
						m = m, n_reservoir = n_reservoir
					)
					kk_boot_draws = replicate(
						as.integer(B),
						des_priv$draw_bootstrap_indices(),
						simplify = FALSE
					)
					if (isTRUE(private$use_reusable_kk_bootstrap_worker())) {
						if (isTRUE(debug)) {
							return(private$compute_kk_bootstrap_debug_with_reused_worker(kk_boot_draws, kk_boot_context))
						}
						actual_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
						if (actual_cores > 1L) {
							do_warmup_iter = function() {
								worker_state = private$create_kk_bootstrap_worker_state(kk_boot_context)
								sample_info = kk_boot_draws[[1L]]
								private$load_kk_bootstrap_sample_into_worker(worker_state, sample_info)
								tryCatch(private$compute_kk_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
							}
							system.time(do_warmup_iter())
							t_boot_warmup = system.time(do_warmup_iter())[[3]]
							fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.5
							if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores)) {
								actual_cores = 1L
							}
						}
						return(private$compute_kk_bootstrap_distribution_with_reused_workers(
							kk_boot_draws = kk_boot_draws,
							kk_boot_context = kk_boot_context,
							actual_cores = actual_cores,
							show_progress = show_progress
						))
					}
					if (isTRUE(debug)) {
						debug_results = vector("list", B)
						has_res_stat_debug = private$has_private_method("compute_reservoir_and_match_statistics")
						for (b in seq_len(B)) {
							boot_sample = kk_boot_draws[[b]]
							i_b   = boot_sample$i_b
							m_vec_b = boot_sample$m_vec_b
							iter_warns = character(0)
							iter_val = withCallingHandlers(
								tryCatch({
									boot_inf_obj = self$duplicate()
									boot_inf_obj$.__enclos_env__$private$y = y[i_b]
									boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
									boot_inf_obj$.__enclos_env__$private$w = w[i_b]
									boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
									boot_inf_obj$.__enclos_env__$private$cached_values = list()
									boot_inf_obj$.__enclos_env__$private$m = m_vec_b
									boot_inf_obj$.__enclos_env__$private$compute_basic_match_data()
									if (has_res_stat_debug) boot_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
									boot_inf_obj$compute_estimate(estimate_only = TRUE)
								}, error = function(e) list(val = NA_real_, error = conditionMessage(e))),
								warning = function(wrn) { iter_warns <<- c(iter_warns, conditionMessage(wrn)); invokeRestart("muffleWarning") }
							)
							debug_results[[b]] = list(
								val = if (is.list(iter_val) && !is.null(iter_val$val)) as.numeric(iter_val$val)[1L] else as.numeric(iter_val)[1L],
								errors = if (is.list(iter_val) && !is.null(iter_val$error)) iter_val$error else character(0),
								warnings = iter_warns
							)
						}
						debug_results = debug_results[!vapply(debug_results, is.null, logical(1))]
						if (length(debug_results) == 0L) {
							stop("All bootstrap iterations failed or returned invalid results. Check for worker crashes or out-of-memory issues.")
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
					# Pure R KK bootstrap implementation
					if (self$num_cores == 1) {
						pb = NULL
						if (private$verbose) {
							pb = utils::txtProgressBar(min = 0, max = B, style = 3)
							on.exit(close(pb), add = TRUE)
						}
						beta_hat_T_bs = rep(NA_real_, B)
						has_res_stat_serial = private$has_private_method("compute_reservoir_and_match_statistics")
						for (b in 1:B) {
							boot_sample = kk_boot_draws[[b]]
							i_b    = boot_sample$i_b
							m_vec_b = boot_sample$m_vec_b
							boot_inf_obj = self$duplicate()
							boot_inf_obj$.__enclos_env__$private$y = y[i_b]
							boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
							boot_inf_obj$.__enclos_env__$private$w = w[i_b]
							boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
							boot_inf_obj$.__enclos_env__$private$cached_values = list()
							beta_hat_T_bs[b] = tryCatch({
								boot_inf_obj$.__enclos_env__$private$m = m_vec_b
								boot_inf_obj$.__enclos_env__$private$compute_basic_match_data()
								if (has_res_stat_serial) boot_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
								boot_inf_obj$compute_estimate(estimate_only = TRUE)
							}, error = function(e) { NA_real_ })
							if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
						}
						return(beta_hat_T_bs)
					} else {
						# Parallel bootstrap execution for KK designs
						cores_to_use = self$num_cores
						if (cores_to_use > 1L) {
							boot_sample_w = des_priv$draw_bootstrap_indices()
							i_b_w   = boot_sample_w$i_b
							m_vec_b_w = boot_sample_w$m_vec_b
							t_warmup_kk = system.time({
								boot_w = self$duplicate()
								boot_w$.__enclos_env__$private$y = y[i_b_w]
								boot_w$.__enclos_env__$private$dead = dead[i_b_w]
								boot_w$.__enclos_env__$private$w = w[i_b_w]
								boot_w$.__enclos_env__$private$X = X[i_b_w, , drop = FALSE]
								boot_w$.__enclos_env__$private$cached_values = list()
								boot_w$.__enclos_env__$private$m = m_vec_b_w
								boot_w$.__enclos_env__$private$compute_basic_match_data()
								if (private$object_has_private_method(boot_w, "compute_reservoir_and_match_statistics"))
									boot_w$.__enclos_env__$private$compute_reservoir_and_match_statistics()
								tryCatch(boot_w$compute_estimate(estimate_only = TRUE), error = function(e) NA_real_)
							})[[3]]
							fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.5
							if (!(t_warmup_kk * B > fork_overhead_estimate * cores_to_use))
								cores_to_use = 1L
						}
						kk_template = self$duplicate()
						has_res_stat = private$object_has_private_method(kk_template, "compute_reservoir_and_match_statistics")
						# Use private$par_lapply which handles both persistent fork clusters and other strategies.
						beta_hat_T_bs = unlist(private$par_lapply(seq_along(kk_boot_draws), function(b) {
							boot_sample = kk_boot_draws[[b]]
							i_b    = boot_sample$i_b
							m_vec_b = boot_sample$m_vec_b
							worker_inf = kk_template$duplicate()
							worker_inf$.__enclos_env__$private$y = y[i_b]
							worker_inf$.__enclos_env__$private$dead = dead[i_b]
							worker_inf$.__enclos_env__$private$w = w[i_b]
							worker_inf$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
							worker_inf$.__enclos_env__$private$cached_values = list()
							tryCatch({
								worker_inf$.__enclos_env__$private$m = m_vec_b
								worker_inf$.__enclos_env__$private$compute_basic_match_data()
								if (has_res_stat) worker_inf$.__enclos_env__$private$compute_reservoir_and_match_statistics()
								worker_inf$compute_estimate(estimate_only = TRUE)
							}, error = function(e) NA_real_)
						}, n_cores = cores_to_use, show_progress = private$verbose,
						export_list = list(
							y = y, dead = dead, w = w, X = X, kk_template = kk_template,
							kk_boot_draws = kk_boot_draws,
							has_res_stat = has_res_stat
						)))
						return(beta_hat_T_bs)
					}
				}
		}
	),
	private = list(
		m = NULL,
		kk_passthrough = TRUE,
		supports_information_preference = function(){
			FALSE
		},
		supports_observed_information = function(){
			FALSE
		},
		get_supported_testing_types_impl = function(){
			"wald"
		},
		get_supported_information_preferences_impl = function(){
			"auto"
		},
		use_reusable_kk_bootstrap_worker = function(){
			FALSE
		},
		init_kk_passthrough = function(des_obj){
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")) {
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14) or DesignFixedBinaryMatch.")
				}
			}
			if (private$has_match_structure){
				if (inherits(des_obj, "DesignFixedBinaryMatch")){
					des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
				}
				private$m = des_obj$.__enclos_env__$private$m
				private$compute_basic_match_data()
			}
		},
		create_kk_bootstrap_context = function(y, dead, w, X, m, n_reservoir){
			X_mat = if (is.null(X)) {
				matrix(numeric(0), nrow = length(y), ncol = 0L)
			} else {
				as.matrix(X)
			}
			list(
				y = as.numeric(y),
				dead = dead,
				w = as.integer(w),
				X = X_mat,
				m = as.integer(m),
				n_reservoir = as.integer(n_reservoir)
			)
		},
		create_kk_bootstrap_worker_state = function(kk_boot_context){
			worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			worker$num_cores = 1L
			worker_priv = worker$.__enclos_env__$private
			list(
				worker = worker,
				worker_priv = worker_priv,
				base_y = kk_boot_context$y,
				base_dead = kk_boot_context$dead,
				base_w = kk_boot_context$w,
				base_X = kk_boot_context$X,
				n_reservoir = kk_boot_context$n_reservoir,
				has_res_stat = private$object_has_private_method(worker, "compute_reservoir_and_match_statistics")
			)
		},
		load_kk_bootstrap_sample_into_worker = function(worker_state, sample_info){
			worker_priv = worker_state$worker_priv
			i_b = sample_info$i_b
			worker_priv$y = worker_state$base_y[i_b]
			worker_priv$y_temp = worker_priv$y
			worker_priv$dead = worker_state$base_dead[i_b]
			worker_priv$w = worker_state$base_w[i_b]
			worker_priv$X = worker_state$base_X[i_b, , drop = FALSE]
			worker_priv$cached_values = list(
				KKstats = compute_bootstrap_matching_stats_cpp(
					X = worker_state$base_X,
					y = worker_state$base_y,
					w = worker_state$base_w,
					i_b = i_b,
					n_reservoir = worker_state$n_reservoir
				)
			)
			worker_priv$best_X_colnames = NULL
			worker_priv$best_Xmm_colnames = NULL
			worker_priv$fit_warm_coefficients = NULL
			worker_priv$cached_mod = NULL
			worker_priv$m = sample_info$m_vec_b
			if (isTRUE(worker_state$has_res_stat)) {
				worker_priv$compute_reservoir_and_match_statistics()
			}
		},
		compute_kk_bootstrap_worker_estimate = function(worker_state){
			as.numeric(worker_state$worker$compute_estimate(estimate_only = TRUE))[1L]
		},
		compute_kk_bootstrap_debug_with_reused_worker = function(kk_boot_draws, kk_boot_context){
			B = length(kk_boot_draws)
			worker_state = private$create_kk_bootstrap_worker_state(kk_boot_context)
			debug_results = vector("list", B)
			for (b in seq_len(B)) {
				iter_warns = character(0)
				iter_val = withCallingHandlers(
					tryCatch({
						sample_info = kk_boot_draws[[b]]
						private$load_kk_bootstrap_sample_into_worker(worker_state, sample_info)
						private$compute_kk_bootstrap_worker_estimate(worker_state)
					}, error = function(e) list(val = NA_real_, error = conditionMessage(e))),
					warning = function(wrn) { iter_warns <<- c(iter_warns, conditionMessage(wrn)); invokeRestart("muffleWarning") }
				)
				debug_results[[b]] = list(
					val = if (is.list(iter_val) && !is.null(iter_val$val)) as.numeric(iter_val$val)[1L] else as.numeric(iter_val)[1L],
					errors = if (is.list(iter_val) && !is.null(iter_val$error)) iter_val$error else character(0),
					warnings = iter_warns
				)
			}
			debug_results = debug_results[!vapply(debug_results, is.null, logical(1))]
			if (should_run_asserts()) {
				if (length(debug_results) == 0L) {
					stop("All bootstrap iterations failed or returned invalid results. Check for worker crashes or out-of-memory issues.")
				}
			}
			values = sapply(debug_results, `[[`, "val")
			errors_list = lapply(debug_results, `[[`, "errors")
			warnings_list = lapply(debug_results, `[[`, "warnings")
			num_errors_vec = lengths(errors_list)
			num_warnings_vec = lengths(warnings_list)
			list(
				values = values,
				errors = errors_list,
				warnings = warnings_list,
				num_errors = num_errors_vec,
				num_warnings = num_warnings_vec,
				prop_iterations_with_errors = mean(num_errors_vec > 0),
				prop_iterations_with_warnings = mean(num_warnings_vec > 0),
				prop_illegal_values = mean(!is.finite(values))
			)
		},
		compute_kk_bootstrap_distribution_with_reused_workers = function(kk_boot_draws, kk_boot_context, actual_cores, show_progress = FALSE){
			B = length(kk_boot_draws)
			chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
			chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
			chunks = split(seq_len(B), chunk_id)
			run_chunk = function(idxs) {
				worker_state = private$create_kk_bootstrap_worker_state(kk_boot_context)
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					sample_info = kk_boot_draws[[idxs[[k]]]]
					out[k] = tryCatch({
						private$load_kk_bootstrap_sample_into_worker(worker_state, sample_info)
						private$compute_kk_bootstrap_worker_estimate(worker_state)
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
		},
		compute_basic_kk_match_data_impl = function(){
			if (!isTRUE(private$has_match_structure)) {
				private$cache_nonestimable_estimate("kk_design_required")
				return(invisible(NULL))
			}
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv     = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},
		compute_concordant_and_discordant_match_statistics = function(){
			m = private$cached_values$KKstats$m
			y_matched_diffs = private$cached_values$KKstats$y_matched_diffs
			i_m_conc = which(y_matched_diffs == 0)
			i_m_disc = setdiff(seq_len(m), i_m_conc)
			private$cached_values$KKstats$i_m_conc = i_m_conc
			private$cached_values$KKstats$i_m_disc = i_m_disc
			private$cached_values$KKstats$n_m_conc = length(i_m_conc)
			private$cached_values$KKstats$n_m_disc = length(i_m_disc)
			private$cached_values$KKstats$y_matched_diffs_disc = y_matched_diffs[i_m_disc]
			private$cached_values$KKstats$X_matched_diffs_disc =
				private$cached_values$KKstats$X_matched_diffs[i_m_disc, drop = FALSE]
			i_conc = which(private$m %in% i_m_conc)
			private$cached_values$KKstats$X_conc = private$get_X()[i_conc, , drop = FALSE]
			private$cached_values$KKstats$y_conc = private$y[i_conc]
			private$cached_values$KKstats$w_conc = private$w[i_conc]
		},
		compute_model_matrix_with_matching_dummies = function(){
			if (is.null(private$cached_values$data_frame_with_matching_dummies)){
				if (!is.null(private$m) & uniqueN(private$m) > 1){
					mm = model.matrix(~ 0 + factor(private$m))
					mm = mm[, 2 : (ncol(mm) - 1)]
				} else {
					mm = NULL
				}
				private$cached_values$data_frame_with_matching_dummies = cbind(data.frame(w = private$w), mm)
			}
			private$cached_values$data_frame_with_matching_dummies
		}
	)
)
