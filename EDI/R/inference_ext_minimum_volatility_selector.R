#' Shared Minimum-Volatility Resample-Size Selector
#'
#' Internal implementation of the Politis/Romano/Wolf and Bickel/Sakov
#' minimum-volatility tuning heuristic. Method-specific extensions provide the
#' evaluator; this extension only owns grid evaluation and selection.
#'
#' Pattern-1 file-split extension, spliced into exactly one class
#' (\code{InferenceNonParamBootstrap}) -- not a reusable mixin.
#'
#' @keywords internal
#' @noRd
InferenceExtMinimumVolatilitySelector = list(
	public = list(),
	private = list(
		select_optimal_resample_size = function(size_grid, evaluator, objective = "ci_width",
				volatility_window = 3L, min_finite_fraction = 0.8, select = TRUE,
				size_name = "size", exponent_grid = NULL, selection_class = "EDIResampleSizeSelection"){
			if (should_run_asserts()) {
				assertIntegerish(size_grid, lower = 1, min.len = 1, any.missing = FALSE)
				assertFunction(evaluator)
				assertChoice(objective, c("ci_width", "pval_stability", "estimate_stability", "finite_fraction_penalized_ci_width"))
				assertCount(volatility_window, positive = TRUE)
				assertNumber(min_finite_fraction, lower = 0, upper = 1)
			}
			if (volatility_window < 3L || volatility_window %% 2L == 0L) {
				stop("volatility_window must be an odd integer >= 3.", call. = FALSE)
			}
			size_grid_input = as.integer(size_grid)
			exponent_input = if (!is.null(exponent_grid)) as.numeric(exponent_grid) else NULL
			size_grid = sort(unique(size_grid_input))
			K = length(size_grid)
			if (K < volatility_window) {
				stop("size_grid must have at least volatility_window unique values.", call. = FALSE)
			}
			rows = vector("list", K)
			for (k in seq_along(size_grid)) {
				size_k = size_grid[[k]]
				t0 = proc.time()[["elapsed"]]
				res_k = tryCatch(evaluator(size_k), error = function(e) {
					list(
						estimate = NA_real_,
						ci = c(NA_real_, NA_real_),
						pval = NA_real_,
						n_finite = 0L,
						finite_fraction = 0,
						dominant_failure_reason = conditionMessage(e),
						status = "error"
					)
				})
				ci = as.numeric(res_k$ci %||% c(NA_real_, NA_real_))
				if (length(ci) < 2L) ci = c(NA_real_, NA_real_)
				ci_width = if (all(is.finite(ci[1:2]))) abs(diff(ci[1:2])) else NA_real_
				objective_value = switch(
					objective,
					ci_width = ci_width,
					pval_stability = as.numeric(res_k$pval %||% NA_real_)[1L],
					estimate_stability = as.numeric(res_k$estimate %||% NA_real_)[1L],
					finite_fraction_penalized_ci_width = if (is.finite(ci_width)) ci_width / max(.Machine$double.eps, as.numeric(res_k$finite_fraction %||% 0)[1L]) else NA_real_
				)
				rows[[k]] = data.frame(
					size = size_k,
					estimate = as.numeric(res_k$estimate %||% NA_real_)[1L],
					ci_lower = ci[1L],
					ci_upper = ci[2L],
					ci_width = ci_width,
					pval = as.numeric(res_k$pval %||% NA_real_)[1L],
					objective_value = objective_value,
					n_finite = as.integer(res_k$n_finite %||% 0L),
					finite_fraction = as.numeric(res_k$finite_fraction %||% 0)[1L],
					dominant_failure_reason = as.character(res_k$dominant_failure_reason %||% NA_character_)[1L],
					status = as.character(res_k$status %||% "ok")[1L],
					elapsed_sec = unname(proc.time()[["elapsed"]] - t0),
					stringsAsFactors = FALSE
				)
			}
			grid_table = do.call(rbind, rows)
			names(grid_table)[names(grid_table) == "size"] = size_name
			if (!is.null(exponent_input)) {
				grid_table[[paste0(size_name, "_pow_of_n")]] = NA_real_
				for (size_i in size_grid) {
					first_idx = match(size_i, size_grid_input)
					grid_table[[paste0(size_name, "_pow_of_n")]][grid_table[[size_name]] == size_i] = exponent_input[first_idx]
				}
			}
			half_window = (volatility_window - 1L) %/% 2L
			volatility = rep(NA_real_, K)
			eligible_idx = (half_window + 1L):(K - half_window)
			for (k in eligible_idx) {
				win = grid_table$objective_value[(k - half_window):(k + half_window)]
				if (all(is.finite(win))) volatility[k] = stats::sd(win)
			}
			grid_table$volatility = volatility
			grid_table$eligible = is.finite(grid_table$volatility) &
				is.finite(grid_table$finite_fraction) &
				grid_table$finite_fraction >= min_finite_fraction &
				grid_table$status == "ok"
			result = list(
				grid_table = grid_table,
				objective = objective,
				volatility_window = as.integer(volatility_window),
				min_finite_fraction = as.numeric(min_finite_fraction)
			)
			if (isTRUE(select)) {
				eligible = which(grid_table$eligible)
				if (!length(eligible)) {
					result[[paste0(size_name, "_optimal")]] = NA_integer_
					result$status = "nonestimable"
					result$reason = paste0(size_name, "_selection_failed")
				} else {
					min_vol = min(grid_table$volatility[eligible])
					tied = eligible[abs(grid_table$volatility[eligible] - min_vol) <= sqrt(.Machine$double.eps)]
					k_star = tied[which.min(grid_table[[size_name]][tied])]
					result[[paste0(size_name, "_optimal")]] = as.integer(grid_table[[size_name]][k_star])
					pow_col = paste0(size_name, "_pow_of_n")
					if (pow_col %in% names(grid_table)) {
						result[[paste0(size_name, "_pow_of_n_optimal")]] = grid_table[[pow_col]][k_star]
					}
					result$status = "ok"
					result$reason = NULL
					result$tie_rule = "smallest_size"
				}
			}
			class(result) = c(selection_class, "list")
			result
		}
	)
)
