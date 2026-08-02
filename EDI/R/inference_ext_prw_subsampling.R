#' Politis/Romano/Wolf Subsampling Extension
#'
#' Public and private methods for PRW subsampling inference. Pattern-1
#' file-split extension, spliced into exactly one class
#' (\code{InferenceNonParamBootstrap}) -- not a reusable mixin.
#'
#' @keywords internal
#' @noRd
InferenceExtPRWSubsampling = list(
	public = list(),
	private = list(
		#' @description Creates the Politis/Romano/Wolf subsampling distribution of
		#'   the treatment-effect estimate.
		#'
		#' @param B Number of subsamples. Default 501.
		#' @param b Number of exchangeable units drawn without replacement. If
		#'   \code{NULL} (default), use the deterministic intermediate-size rule
		#'   \code{floor(n_units^0.7)}, where \code{n_units} is the number of
		#'   exchangeable units used by the design (observations, clusters, pairs, or
		#'   matched sets). The resolved value must satisfy
		#'   \code{max(5, p_eff + 2) <= b <= floor(n_units / 2)}.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param debug If \code{TRUE}, return distribution diagnostics in addition
		#'   to the subsampled estimates.
		#' @param subsampling_type Optional empirical-resampling scheme. See
		#'   \code{approximate_bootstrap_distribution_beta_hat_T()}.
		#' @param scaling Scaling sequence for centered subsampling pivots. The
		#'   default \code{"sqrt_n"} uses \code{sqrt(b)} for the subsample
		#'   distribution and converts back to the full-sample scale using
		#'   \code{sqrt(n_units)}.
		#' @param center Centering convention for diagnostics and cache keys.
		#'
		#' @details The default \code{b = NULL} rule is a cheap deterministic
		#'   intermediate sequence: \eqn{b \to \infty} and \eqn{b / n \to 0}, as
		#'   required by the Politis, Romano, and Wolf subsampling framework. The
		#'   exponent 0.7 is a pragmatic interior point in \eqn{(0, 1)}; it is not a
		#'   universal optimum. Use \code{select_optimal_b_subsampling()} for
		#'   data-adaptive minimum-volatility selection.
		#'
		#' @return A numeric vector of subsampled estimates, or when
		#'   \code{debug = TRUE}, a diagnostic list.
		#'
		#' @references Politis, D. N., Romano, J. P., and Wolf, M. (1999).
		#'   \emph{Subsampling}. Springer.
		approximate_subsampling_distribution_beta_hat_T_impl = function(B = 501, b = NULL, show_progress = TRUE, debug = FALSE, subsampling_type = NULL, scaling = "sqrt_n", center = "full_estimate"){
			private$active_resampling_operation = "subsampling"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("PRW subsampling inference")
				private$assert_valid_bootstrap_type(subsampling_type)
				assertCount(B, positive = TRUE)
				assertFlag(show_progress)
				assertFlag(debug)
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = subsampling_type)
			if (is.list(b)) {
				selection = do.call(self$select_optimal_b_subsampling, c(list(
					B = as.integer(b$B %||% min(as.integer(B), 251L)),
					alpha = as.numeric(b$alpha %||% 0.05),
					subsampling_type = subsampling_type,
					scaling = scaling,
					show_progress = show_progress
				), b[setdiff(names(b), c("method", "B", "alpha"))]))
				b = selection$b_optimal
				if (!is.finite(b)) stop("PRW subsampling b selection failed.", call. = FALSE)
			}
			b = private$resolve_resampling_size(b, unit_info$n_units, "b")
			cache_key = private$subsampling_cache_key(B = B, b = b, subsampling_type = subsampling_type, scaling = scaling, center = center)
			if (!isTRUE(debug)) {
				private$ensure_resampling_distribution_cache("subsampling")
				cached = private$get_cached_resampling_distribution("subsampling", cache_key)
				if (!is.null(cached)) return(cached)
			}
			if (!is.null(private$seed)) {
				had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				if (had_seed) old_seed = get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
				on.exit(
					if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else rm(".Random.seed", envir = .GlobalEnv),
					add = TRUE
				)
				set.seed(private$seed)
			}
			draws = private$generate_exchangeable_resampling_draws(
				unit_info = unit_info,
				B = B,
				size = b,
				replace = FALSE,
				stratified = !identical(subsampling_type, "resample_blocks"),
				preserve_order = TRUE,
				size_label = "b"
			)
			res = private$compute_resampling_draw_distribution(
				draws = draws,
				operation = "subsampling",
				show_progress = show_progress,
				debug = debug,
				cache_key = cache_key
			)
			if (isTRUE(debug)) {
				est = tryCatch(as.numeric(self$compute_estimate(estimate_only = TRUE))[1L], error = function(e) NA_real_)
				res$centered_scaled_values = private$resampling_scaling_factor(b, scaling) * (res$values - est)
				res$full_estimate = est
				res$b = b
				res$n_units = unit_info$n_units
				res$scaling = scaling
				res$unit_type = unit_info$unit_type
			}
			res
		},
		#' @description Computes a centered PRW subsampling two-sided p-value.
		#'
		#' @param delta Null treatment effect. Default 0.
		#' @param B Number of subsamples. Default 501.
		#' @param b Number of exchangeable units drawn without replacement. If
		#'   \code{NULL} (default), use \code{floor(n_units^0.7)} subject to the
		#'   validation bounds documented for
		#'   \code{approximate_subsampling_distribution_beta_hat_T()}.
		#' @param type P-value type. Currently only \code{"centered"} is supported.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param min_number_usable_samples Minimum number of finite subsampled
		#'   estimates required after filtering. Default 5.
		#' @param subsampling_type Optional empirical-resampling scheme.
		#' @param scaling Scaling sequence for centered subsampling pivots.
		#'
		#' @details The default \code{b = NULL} follows the intermediate-sequence
		#'   convention from the Politis/Romano/Wolf subsampling literature:
		#'   \eqn{b \to \infty} and \eqn{b / n \to 0}. The deterministic exponent
		#'   0.7 is a first-pass default; for unstable paths prefer the
		#'   minimum-volatility selector.
		#'
		#' @return A numeric two-sided p-value, or \code{NA_real_} if the path is
		#'   non-estimable.
		compute_subsampling_two_sided_pval_impl = function(delta = 0, B = 501, b = NULL, type = "centered", show_progress = TRUE, min_number_usable_samples = 5L, subsampling_type = NULL, scaling = "sqrt_n"){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertChoice(type, c("centered"))
			}
			if (is.list(b)) {
				selection = do.call(self$select_optimal_b_subsampling, c(list(
					B = as.integer(b$B %||% min(as.integer(B), 251L)),
					alpha = as.numeric(b$alpha %||% 0.05),
					subsampling_type = subsampling_type,
					scaling = scaling,
					show_progress = show_progress
				), b[setdiff(names(b), c("method", "B", "alpha"))]))
				b = selection$b_optimal
				if (!is.finite(b)) {
					if (isTRUE(private$harden)) private$cache_nonestimable_estimate("subsampling_b_selection_failed")
					return(NA_real_)
				}
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = subsampling_type)
			b = private$resolve_resampling_size(b, unit_info$n_units, "b")
			pivot = private$subsampling_centered_pivot(
				B = B,
				b = b,
				unit_info = unit_info,
				show_progress = show_progress,
				subsampling_type = subsampling_type,
				scaling = scaling
			)
			if (!isTRUE(pivot$ok)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate(pivot$reason)
				return(NA_real_)
			}
			if (length(pivot$finite) < as.integer(min_number_usable_samples)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate("subsampling_too_few_finite_estimates")
				return(NA_real_)
			}
			private$resampling_centered_pval(pivot$centered_scaled, est = pivot$est, delta = delta, n_units = pivot$n_units, scaling = scaling)
		},
		#' @description Computes a basic PRW subsampling confidence interval.
		#'
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of subsamples. Default 501.
		#' @param b Number of exchangeable units drawn without replacement. If
		#'   \code{NULL} (default), use \code{floor(n_units^0.7)} subject to the
		#'   validation bounds documented for
		#'   \code{approximate_subsampling_distribution_beta_hat_T()}.
		#' @param type Confidence-interval type. Currently only \code{"basic"} is
		#'   supported.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param min_number_usable_samples Minimum number of finite subsampled
		#'   estimates required after filtering. Default 5.
		#' @param subsampling_type Optional empirical-resampling scheme.
		#' @param scaling Scaling sequence for centered subsampling pivots.
		#'
		#' @details The \code{NULL} default is grounded in the standard
		#'   Politis/Romano/Wolf asymptotic condition \eqn{b \to \infty} and
		#'   \eqn{b / n \to 0}. The minimum-volatility selector is available when a
		#'   fixed deterministic exponent is too brittle for a specific
		#'   estimator/design path.
		#'
		#' @return A length-2 numeric confidence interval, or
		#'   \code{c(NA_real_, NA_real_)} if the path is non-estimable.
		compute_subsampling_confidence_interval_impl = function(alpha = 0.05, B = 501, b = NULL, type = "basic", show_progress = TRUE, min_number_usable_samples = 5L, subsampling_type = NULL, scaling = "sqrt_n"){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertChoice(type, c("basic"))
			}
			if (is.list(b)) {
				selection = do.call(self$select_optimal_b_subsampling, c(list(
					B = as.integer(b$B %||% min(as.integer(B), 251L)),
					alpha = alpha,
					subsampling_type = subsampling_type,
					scaling = scaling,
					show_progress = show_progress
				), b[setdiff(names(b), c("method", "B", "alpha"))]))
				b = selection$b_optimal
				if (!is.finite(b)) {
					if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, "subsampling_b_selection_failed", stage = "estimate"))
					return(c(NA_real_, NA_real_))
				}
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = subsampling_type)
			b = private$resolve_resampling_size(b, unit_info$n_units, "b")
			pivot = private$subsampling_centered_pivot(
				B = B,
				b = b,
				unit_info = unit_info,
				show_progress = show_progress,
				subsampling_type = subsampling_type,
				scaling = scaling
			)
			if (!isTRUE(pivot$ok)) {
				if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, pivot$reason, stage = "estimate"))
				return(c(NA_real_, NA_real_))
			}
			if (length(pivot$finite) < as.integer(min_number_usable_samples)) {
				if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, "subsampling_too_few_finite_estimates", stage = "estimate"))
				return(c(NA_real_, NA_real_))
			}
			ci = private$resampling_ci_from_centered_distribution(pivot$centered_scaled, alpha = alpha, est = pivot$est, n_units = pivot$n_units, scaling = scaling)
			if (!all(is.finite(ci[1:2])) || private$bootstrap_confidence_interval_extreme(ci, est = pivot$est)) {
				if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, "subsampling_extreme_confidence_interval", stage = "estimate"))
			}
			ci
		},
		#' @description Selects a PRW subsampling size by minimum volatility.
		#'
		#' @param B Number of subsamples per candidate size. Default 251.
		#' @param alpha Significance level for interval-width objectives.
		#' @param b_pow_of_n_grid Candidate exponent grid used when
		#'   \code{b_grid = NULL}. Defaults to \code{seq(0.5, 0.9, by = 0.05)}.
		#' @param b_grid Optional explicit integer candidate sizes.
		#' @param objective Selection objective. Currently \code{"ci_width"}.
		#' @param target Target summary. Currently \code{"ci"}.
		#' @param volatility_window Rolling window size used to measure local
		#'   volatility across candidate sizes.
		#' @param subsampling_type Optional empirical-resampling scheme.
		#' @param scaling Scaling sequence for centered subsampling pivots.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param min_finite_fraction Minimum finite-subsample fraction required for
		#'   a candidate size to be eligible.
		#'
		#' @details This implements the same minimum-volatility idea used in PTE's
		#'   m-selection workflow, applied to the PRW block/subsample-size choice:
		#'   scan admissible intermediate sizes and choose a stable region of the
		#'   target statistic rather than assuming one exponent is uniformly optimal.
		#'
		#' @return An \code{EDISubsamplingBSelection} list with the selected
		#'   \code{b}, mapped exponent, candidate table, status, and reason.
		select_optimal_b_subsampling_impl = function(B = 251, alpha = 0.05, b_pow_of_n_grid = seq(0.5, 0.9, by = 0.05), b_grid = NULL, objective = "ci_width", target = "ci", volatility_window = 3L, subsampling_type = NULL, scaling = "sqrt_n", show_progress = TRUE, min_finite_fraction = 0.8){
			if (should_run_asserts()) {
				assertCount(B, positive = TRUE)
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertChoice(target, c("ci"))
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = subsampling_type)
			if (is.null(b_grid)) {
				b_grid_raw = floor(unit_info$n_units ^ b_pow_of_n_grid)
				grid_df = unique(data.frame(b = as.integer(b_grid_raw), pow = as.numeric(b_pow_of_n_grid), stringsAsFactors = FALSE))
				b_grid = grid_df$b
				b_pow_mapped = grid_df$pow
			} else {
				b_grid = as.integer(b_grid)
				b_pow_mapped = rep(NA_real_, length(b_grid))
			}
			evaluator = function(b_size) private$evaluate_subsampling_size(b = b_size, B = B, alpha = alpha, subsampling_type = subsampling_type, scaling = scaling, show_progress = FALSE)
			private$select_optimal_resample_size(
				size_grid = b_grid,
				evaluator = evaluator,
				objective = objective,
				volatility_window = volatility_window,
				min_finite_fraction = min_finite_fraction,
				size_name = "b",
				exponent_grid = b_pow_mapped,
				selection_class = "EDISubsamplingBSelection"
			)
		},
		#' @description Computes PRW subsampling sensitivity over candidate sizes.
		#'
		#' @param B Number of subsamples per candidate size. Default 251.
		#' @param alpha Significance level for interval-width objectives.
		#' @param b_pow_of_n_grid Candidate exponent grid used when
		#'   \code{b_grid = NULL}. Defaults to \code{seq(0.5, 0.9, by = 0.05)}.
		#' @param b_grid Optional explicit integer candidate sizes.
		#' @param objective Selection objective. Currently \code{"ci_width"}.
		#' @param target Target summary. Currently \code{"ci"}.
		#' @param volatility_window Rolling window size used to measure local
		#'   volatility across candidate sizes.
		#' @param subsampling_type Optional empirical-resampling scheme.
		#' @param scaling Scaling sequence for centered subsampling pivots.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param min_finite_fraction Minimum finite-subsample fraction required for
		#'   a candidate size to be eligible. Defaults to 0 for sensitivity scans.
		#'
		#' @return An \code{EDISubsamplingSensitivity} list containing the candidate
		#'   grid table without selecting a final \code{b}.
		compute_subsampling_sensitivity_impl = function(B = 251, alpha = 0.05, b_pow_of_n_grid = seq(0.5, 0.9, by = 0.05), b_grid = NULL, objective = "ci_width", target = "ci", volatility_window = 3L, subsampling_type = NULL, scaling = "sqrt_n", show_progress = TRUE, min_finite_fraction = 0){
			res = self$select_optimal_b_subsampling(
				B = B,
				alpha = alpha,
				b_pow_of_n_grid = b_pow_of_n_grid,
				b_grid = b_grid,
				objective = objective,
				target = target,
				volatility_window = volatility_window,
				subsampling_type = subsampling_type,
				scaling = scaling,
				show_progress = show_progress,
				min_finite_fraction = min_finite_fraction
			)
			res$status = NULL
			res$reason = NULL
			res$b_optimal = NULL
			res$b_pow_of_n_optimal = NULL
			class(res) = c("EDISubsamplingSensitivity", "list")
			res
		}
		,
		subsampling_cache_key = function(B, b, subsampling_type = NULL, scaling = "sqrt_n", center = "full_estimate"){
			paste(as.integer(B), as.integer(b), subsampling_type %||% "NULL", private$resampling_scaling_key(scaling), center, sep = "|")
		},
		subsampling_centered_pivot = function(B, b, unit_info, show_progress = TRUE, subsampling_type = NULL, scaling = "sqrt_n"){
			cache_key = private$subsampling_cache_key(B = B, b = b, subsampling_type = subsampling_type, scaling = scaling, center = "full_estimate")
			cached = private$get_cached_centered_resampling_pivot("subsampling", cache_key)
			if (!is.null(cached)) return(cached)
			est = tryCatch(as.numeric(self$compute_estimate(estimate_only = TRUE))[1L], error = function(e) NA_real_)
			if (!is.finite(est)) {
				return(list(ok = FALSE, reason = "subsampling_original_estimate_unavailable"))
			}
			sub = self$approximate_subsampling_distribution_beta_hat_T(
				B = B,
				b = b,
				show_progress = show_progress,
				subsampling_type = subsampling_type,
				scaling = scaling
			)
			finite = sub[is.finite(sub)]
			pivot = list(
				ok = TRUE,
				reason = NA_character_,
				est = est,
				finite = finite,
				centered_scaled = private$resampling_scaling_factor(b, scaling) * (finite - est),
				b = as.integer(b),
				n_units = unit_info$n_units,
				scaling = scaling
			)
			private$set_cached_centered_resampling_pivot("subsampling", cache_key, pivot)
		},
		subsampling_sample_indices = function(b, subsampling_type = NULL){
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = subsampling_type)
			selected = private$sample_exchangeable_unit_ids(
				unit_info = unit_info,
				size = as.integer(b),
				replace = FALSE,
				stratified = !identical(subsampling_type, "resample_blocks")
			)
			private$build_resampling_draw_from_units(unit_info, selected, size_label = "b", preserve_order = TRUE)
		},
		load_subsampling_draw_into_worker = function(worker_state, draw, ...){
			private$load_bootstrap_sample_into_worker(worker_state, draw)
			invisible(worker_state)
		},
		compute_subsampling_worker_estimate = function(worker_state){
			private$estimate_bootstrap_worker(worker_state)
		},
		evaluate_subsampling_size = function(b, B, alpha, subsampling_type = NULL, scaling = "sqrt_n", show_progress = FALSE){
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = subsampling_type)
			b = private$resolve_resampling_size(b, unit_info$n_units, "b")
			est = tryCatch(as.numeric(self$compute_estimate(estimate_only = TRUE))[1L], error = function(e) NA_real_)
			if (!is.finite(est)) {
				return(list(status = "nonestimable", dominant_failure_reason = "subsampling_original_estimate_unavailable", finite_fraction = 0, n_finite = 0L))
			}
			sub = self$approximate_subsampling_distribution_beta_hat_T(B = B, b = b, show_progress = show_progress, subsampling_type = subsampling_type, scaling = scaling)
			finite = sub[is.finite(sub)]
			centered_scaled = private$resampling_scaling_factor(b, scaling) * (finite - est)
			ci = private$resampling_ci_from_centered_distribution(centered_scaled, alpha = alpha, est = est, n_units = unit_info$n_units, scaling = scaling)
			pval = private$resampling_centered_pval(centered_scaled, est = est, delta = 0, n_units = unit_info$n_units, scaling = scaling)
			ci_ok = length(finite) > 0L &&
				all(is.finite(ci[1:2])) &&
				!private$bootstrap_confidence_interval_extreme(ci, est = est)
			list(
				status = if (ci_ok) "ok" else "nonestimable",
				estimate = est,
				ci = ci,
				pval = pval,
				n_finite = length(finite),
				finite_fraction = length(finite) / as.integer(B),
				dominant_failure_reason = if (!length(finite)) {
					"subsampling_too_few_finite_estimates"
				} else if (!isTRUE(ci_ok)) {
					"subsampling_extreme_confidence_interval"
				} else {
					NA_character_
				}
			)
		}
	)
)
