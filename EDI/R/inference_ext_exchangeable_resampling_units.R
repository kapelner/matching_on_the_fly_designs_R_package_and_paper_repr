#' Shared Exchangeable Resampling Unit Helpers
#'
#' Internal helpers for resampling methods that operate on empirical units.
#' m-out-of-n bootstrap and PRW subsampling share unit discovery but differ in
#' replacement semantics and calibration.
#'
#' Pattern-1 file-split extension, spliced into exactly one class
#' (\code{InferenceNonParamBootstrap}) -- not a reusable mixin.
#'
#' @keywords internal
#' @noRd
InferenceExtExchangeableResamplingUnits = list(
	public = list(),
	private = list(
		resolve_resampling_unit = function(unit = "auto"){
			unit = tolower(as.character(unit)[1L])
			valid = c("auto", "observation", "cluster", "block", "pair", "matched_set")
			if (should_run_asserts()) assertChoice(unit, valid)
			if (unit != "auto") return(unit)
			design_obj = private$des_obj
			is_matching_design = isTRUE(design_obj$is_matching_design())
			is_cluster_design = design_obj$is_a_cluster_capable()
			is_blocking_design = isTRUE(design_obj$is_blocking_design())
			if (is_matching_design) {
				if (isTRUE(private$is_KK)) "matched_set" else "pair"
			} else if (is_cluster_design) {
				"cluster"
			} else if (is_blocking_design) {
				"observation"
			} else {
				"observation"
			}
		},
		get_exchangeable_units = function(unit = "auto", resampling_type = NULL){
			design_obj = private$des_obj
			n = if (!is.null(design_obj)) design_obj$get_n() else private$n
			n = as.integer(n)
			if (!is.finite(n) || n <= 0L) {
				stop("Exchangeable resampling units are unavailable.", call. = FALSE)
			}
			unit = private$resolve_resampling_unit(unit)
			des_priv = private$des_obj_priv_int
			all_idx = seq_len(n)
			make_info = function(units, unit_type, strata_ids = NULL, unit_kind = NULL, m_vec_full = NULL){
				list(
					units = units,
					unit_type = unit_type,
					strata_ids = strata_ids,
					unit_kind = unit_kind,
					m_vec_full = m_vec_full,
					n_units = length(units)
				)
			}
			if (unit %in% c("pair", "matched_set")) {
				if (is.null(des_priv) || !is.function(des_priv$init_matching_bootstrap_structure)) {
					stop("Matching structure is unavailable for resampling.", call. = FALSE)
				}
				des_priv$init_matching_bootstrap_structure()
				pair_rows = des_priv$boot_pair_rows
				i_reservoir = des_priv$boot_i_reservoir
				units = list()
				unit_kind = character(0)
				if (!is.null(pair_rows) && nrow(pair_rows) > 0L) {
					for (i in seq_len(nrow(pair_rows))) {
						units[[length(units) + 1L]] = as.integer(pair_rows[i, ])
						unit_kind = c(unit_kind, "pair")
					}
				}
				if (!is.null(i_reservoir) && length(i_reservoir) > 0L) {
					for (i in seq_along(i_reservoir)) {
						units[[length(units) + 1L]] = as.integer(i_reservoir[i])
						unit_kind = c(unit_kind, "reservoir")
					}
				}
				if (!length(units)) stop("Matching structure has no exchangeable units.", call. = FALSE)
				return(make_info(units, unit, unit_kind = unit_kind, m_vec_full = private$m))
			}
			if (unit == "cluster") {
				cluster_ids = private$get_resampling_cluster_ids()
				if (is.null(cluster_ids)) stop("Cluster IDs are unavailable for resampling.", call. = FALSE)
				group_vals = unique(cluster_ids[!is.na(cluster_ids)])
				units = lapply(group_vals, function(g) all_idx[cluster_ids == g])
				strata_ids = private$get_resampling_strata_ids()
				unit_strata = NULL
				if (!is.null(strata_ids)) {
					unit_strata = vapply(units, function(idx) as.character(strata_ids[idx[1L]]), character(1))
				}
				return(make_info(units, "cluster", strata_ids = unit_strata))
			}
			if (unit == "block" || identical(resampling_type, "resample_blocks")) {
				block_ids = private$get_resampling_block_ids()
				if (is.null(block_ids)) stop("Block IDs are unavailable for resampling.", call. = FALSE)
				block_vals = sort(unique(block_ids[is.finite(block_ids) & block_ids > 0L]))
				units = lapply(block_vals, function(g) all_idx[block_ids == g])
				return(make_info(units, "block"))
			}
			strata_ids = private$get_resampling_strata_ids()
			make_info(lapply(all_idx, function(i) i), "observation", strata_ids = strata_ids)
		},
		get_resampling_cluster_ids = function(){
			des_priv = private$des_obj_priv_int
			if (is.null(des_priv) || is.null(des_priv$cluster_col) || is.null(des_priv$Xraw)) return(NULL)
			n = private$des_obj$get_n()
			as.character(des_priv$Xraw[seq_len(n), ][[des_priv$cluster_col]])
		},
		get_resampling_strata_ids = function(){
			des_priv = private$des_obj_priv_int
			if (is.null(des_priv) || !is.function(des_priv$get_strata_keys)) return(NULL)
			tryCatch(as.character(des_priv$get_strata_keys()), error = function(e) NULL)
		},
		get_resampling_block_ids = function(){
			if (is.null(private$des_obj) || !is.function(private$des_obj$get_block_ids)) return(NULL)
			tryCatch(as.integer(private$des_obj$get_block_ids()), error = function(e) NULL)
		},
		resampling_effective_p = function(){
			X = tryCatch(private$get_X(), error = function(e) private$X)
			if (is.null(X)) return(1L)
			max(1L, as.integer(NCOL(X)))
		},
		resolve_resampling_size = function(size, n_units, label, deterministic_exponent = 0.7){
			if (is.list(size)) stop(label, " selection must be resolved before drawing.", call. = FALSE)
			if (is.null(size)) {
				size = floor(as.numeric(n_units) ^ deterministic_exponent)
			}
			size = as.integer(size)[1L]
			min_size = max(5L, private$resampling_effective_p() + 2L)
			max_size = floor(as.integer(n_units) / 2L)
			if (!is.finite(size) || size < min_size || size > max_size) {
				stop(label, " must satisfy ", min_size, " <= ", label, " <= ", max_size, ".", call. = FALSE)
			}
			size
		},
		allocate_resampling_sizes_by_stratum = function(total_size, strata_ids, replace){
			strata_ids = as.character(strata_ids)
			strata_vals = unique(strata_ids)
			counts = tabulate(match(strata_ids, strata_vals), nbins = length(strata_vals))
			raw = total_size * counts / sum(counts)
			alloc = floor(raw)
			remaining = as.integer(total_size) - sum(alloc)
			if (remaining > 0L) {
				ord = order(raw - alloc, decreasing = TRUE)
				alloc[ord[seq_len(remaining)]] = alloc[ord[seq_len(remaining)]] + 1L
			}
			if (!isTRUE(replace) && any(alloc > counts)) {
				stop("Subsampling stratum is too small for the requested subsample size.", call. = FALSE)
			}
			names(alloc) = strata_vals
			alloc
		},
		sample_exchangeable_unit_ids = function(unit_info, size, replace, stratified = TRUE){
			n_units = length(unit_info$units)
			if (n_units == 0L) stop("No exchangeable units are available.", call. = FALSE)
			if (isTRUE(stratified) && !is.null(unit_info$strata_ids)) {
				alloc = private$allocate_resampling_sizes_by_stratum(size, unit_info$strata_ids, replace = replace)
				out = integer(0)
				for (stratum in names(alloc)) {
					k = as.integer(alloc[[stratum]])
					if (k <= 0L) next
					ids = which(unit_info$strata_ids == stratum)
					out = c(out, sample(ids, k, replace = replace))
				}
				return(as.integer(out))
			}
			sample.int(n_units, size, replace = replace)
		},
		build_resampling_draw_from_units = function(unit_info, selected_unit_ids, size_label, preserve_order = FALSE){
			selected_unit_ids = as.integer(selected_unit_ids)
			if (!length(selected_unit_ids)) return(list(i_b = integer(0), m_vec_b = NULL))
			rows = as.integer(unlist(unit_info$units[selected_unit_ids], use.names = FALSE))
			m_vec_b = NULL
			if (!is.null(unit_info$unit_kind)) {
				m_out = integer(0)
				pair_id = 0L
				for (unit_id in selected_unit_ids) {
					unit_rows = unit_info$units[[unit_id]]
					if (identical(unit_info$unit_kind[[unit_id]], "pair")) {
						pair_id = pair_id + 1L
						m_out = c(m_out, rep(pair_id, length(unit_rows)))
					} else {
						m_out = c(m_out, rep(0L, length(unit_rows)))
					}
				}
				m_vec_b = as.integer(m_out)
			} else if (!is.null(unit_info$m_vec_full)) {
				m_vec_b = private$renumber_match_ids(unit_info$m_vec_full[rows])
			}
			if (isTRUE(preserve_order)) {
				ord = order(rows)
				rows = rows[ord]
				if (!is.null(m_vec_b)) m_vec_b = m_vec_b[ord]
			}
			draw = list(i_b = rows, m_vec_b = m_vec_b)
			draw$unit_ids = selected_unit_ids
			draw$unit_type = unit_info$unit_type
			draw$n_units = unit_info$n_units
			draw[[size_label]] = length(selected_unit_ids)
			draw
		},
		generate_exchangeable_resampling_draws = function(unit_info, B, size, replace, stratified = TRUE, preserve_order = FALSE, size_label){
			B = as.integer(B)
			size = as.integer(size)
			draws = tryCatch(
				exchangeable_resampling_draws_cpp(
					units = unit_info$units,
					strata_ids_sexp = unit_info$strata_ids %||% NULL,
					unit_kind_sexp = unit_info$unit_kind %||% NULL,
					m_vec_full_sexp = unit_info$m_vec_full %||% NULL,
					B = B,
					size = size,
					replace = isTRUE(replace),
					stratified = isTRUE(stratified),
					preserve_order = isTRUE(preserve_order),
					unit_type = unit_info$unit_type,
					size_label = size_label
				),
				error = function(e) NULL
			)
			if (!is.null(draws)) return(draws)
			replicate(B, {
				selected = private$sample_exchangeable_unit_ids(
					unit_info = unit_info,
					size = size,
					replace = replace,
					stratified = stratified
				)
				private$build_resampling_draw_from_units(unit_info, selected, size_label = size_label, preserve_order = preserve_order)
			}, simplify = FALSE)
		},
		compute_resampling_draw_distribution = function(draws, operation, show_progress = TRUE, debug = FALSE, cache_key = NULL){
			private$ensure_resampling_distribution_cache(operation)
			if (!isTRUE(debug) && !is.null(cache_key)) {
				cached = private$get_cached_resampling_distribution(operation, cache_key)
				if (!is.null(cached)) return(cached)
			}
			n_draws = length(draws)
			if (n_draws == 0L) return(numeric(0))
			if (isTRUE(debug)) {
				debug_results = lapply(draws, function(draw) {
					iter_warns = character(0)
					iter_errs = character(0)
					iter_val = withCallingHandlers(
						tryCatch({
							sub_inf = private$bootstrap_subset_inference(draw, smooth = FALSE)
							if (is.null(sub_inf)) NA_real_ else as.numeric(sub_inf$compute_estimate(estimate_only = TRUE))[1L]
						}, error = function(e) { iter_errs <<- c(iter_errs, conditionMessage(e)); NA_real_ }),
						warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
					)
					list(val = as.numeric(iter_val)[1L], errors = iter_errs, warnings = iter_warns)
				})
				values = vapply(debug_results, function(x) x$val, numeric(1))
				errors_list = lapply(debug_results, `[[`, "errors")
				warnings_list = lapply(debug_results, `[[`, "warnings")
				num_errors_vec = lengths(errors_list)
				num_warnings_vec = lengths(warnings_list)
				if (!is.null(cache_key)) private$set_cached_resampling_distribution(operation, cache_key, values)
				return(list(
					values = values,
					errors = errors_list,
					warnings = warnings_list,
					num_errors = num_errors_vec,
					num_warnings = num_warnings_vec,
					prop_iterations_with_errors = mean(num_errors_vec > 0),
					prop_iterations_with_warnings = mean(num_warnings_vec > 0),
					prop_illegal_values = mean(!is.finite(values)),
					finite_fraction = mean(is.finite(values))
				))
			}
			actual_cores = private$effective_parallel_cores(operation, self$num_cores)
			values = if (isTRUE(private$use_reusable_bootstrap_worker())) {
				private$compute_reusable_bootstrap_worker_distribution(
					draws = draws,
					actual_cores = actual_cores,
					show_progress = show_progress,
					operation = operation
				)
			} else {
				as.numeric(unlist(private$par_lapply(seq_along(draws), function(i) {
					sub_inf = private$bootstrap_subset_inference(draws[[i]], smooth = FALSE)
					if (is.null(sub_inf)) return(NA_real_)
					tryCatch({
						theta = as.numeric(sub_inf$compute_estimate(estimate_only = TRUE))[1L]
						if (is.finite(theta)) theta else NA_real_
					}, error = function(e) NA_real_)
				}, n_cores = actual_cores, show_progress = show_progress), use.names = FALSE))
			}
			if (!is.numeric(values)) values = as.numeric(values)
			if (!is.null(cache_key)) private$set_cached_resampling_distribution(operation, cache_key, values)
			values
		},
		get_cached_centered_resampling_pivot = function(operation, cache_key){
			cache_name = paste0(operation, "_centered_pivot_cache")
			cache = private$cached_values[[cache_name]]
			if (is.null(cache)) return(NULL)
			cache[[cache_key]]
		},
		set_cached_centered_resampling_pivot = function(operation, cache_key, value){
			cache_name = paste0(operation, "_centered_pivot_cache")
			if (is.null(private$cached_values[[cache_name]])) {
				private$cached_values[[cache_name]] = list()
			}
			private$cached_values[[cache_name]][[cache_key]] = value
			invisible(value)
		},
		resampling_ci_from_centered_distribution = function(centered_scaled, alpha, est, n_units, scaling = "sqrt_n"){
			centered_scaled = centered_scaled[is.finite(centered_scaled)]
			if (!length(centered_scaled)) return(c(NA_real_, NA_real_))
			tau_n = private$resampling_scaling_factor(n_units, scaling)
			q = stats::quantile(centered_scaled, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8)
			ci = c(est - q[2L] / tau_n, est - q[1L] / tau_n)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		resampling_centered_pval = function(centered_scaled, est, delta, n_units, scaling = "sqrt_n"){
			centered_scaled = centered_scaled[is.finite(centered_scaled)]
			if (!length(centered_scaled)) return(NA_real_)
			tau_n = private$resampling_scaling_factor(n_units, scaling)
			t_obs = tau_n * (as.numeric(est)[1L] - as.numeric(delta)[1L])
			p_raw = 2 * min(mean(centered_scaled <= t_obs), mean(centered_scaled >= t_obs))
			min(1, max(2 / length(centered_scaled), p_raw))
		},
		resampling_scaling_factor = function(size, scaling = "sqrt_n"){
			if (identical(scaling, "sqrt_n")) return(sqrt(as.numeric(size)[1L]))
			if (is.list(scaling) && !is.null(scaling$rate_exponent)) {
				return(as.numeric(size)[1L] ^ as.numeric(scaling$rate_exponent)[1L])
			}
			if (is.function(scaling)) return(as.numeric(scaling(size))[1L])
			stop("Unsupported resampling scaling.", call. = FALSE)
		}
	)
)
