#' m-out-of-n Bootstrap Extension
#'
#' Public and private methods for m-out-of-n bootstrap inference. Pattern-1
#' file-split extension, spliced into exactly one class
#' (\code{InferenceNonParamBootstrap}) -- not a reusable mixin.
#'
#' @keywords internal
#' @noRd
InferenceExtMOutOfNBootstrap = list(
	public = list(),
	private = list(
		#' @description Creates the m-out-of-n bootstrap distribution of the
		#'   treatment-effect estimate.
		#'
		#' @param B Number of resamples. Default 501.
		#' @param m Number of exchangeable resampling units drawn with replacement.
		#'   If \code{NULL} (default), use the deterministic intermediate-size rule
		#'   \code{floor(n_units^0.7)}, where \code{n_units} is the number of
		#'   exchangeable units used by the design (observations, clusters, pairs, or
		#'   matched sets). The resolved value must satisfy
		#'   \code{max(5, p_eff + 2) <= m <= floor(n_units / 2)}.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param debug If \code{TRUE}, return distribution diagnostics in addition
		#'   to the resampled estimates.
		#' @param bootstrap_type Optional empirical-resampling scheme. See
		#'   \code{approximate_bootstrap_distribution_beta_hat_T()}.
		#' @param scaling Scaling sequence for centered m-out-of-n pivots. The
		#'   default \code{"sqrt_n"} uses \code{sqrt(m)} for the m-sample
		#'   distribution and converts back to the full-sample scale using
		#'   \code{sqrt(n_units)}.
		#' @param center Centering convention for diagnostics and cache keys.
		#'
		#' @details The default \code{m = NULL} rule is a cheap deterministic
		#'   intermediate sequence: \eqn{m \to \infty} and \eqn{m / n \to 0}, as
		#'   required by the standard m-out-of-n bootstrap asymptotic setup
		#'   (Bickel, Gotze, and van Zwet; Bickel and Sakov). The exponent 0.7 is a
		#'   pragmatic interior point in \eqn{(0, 1)}; it is not a silver-bullet
		#'   optimal choice. Use \code{select_optimal_m_out_of_n_bootstrap()} for
		#'   data-adaptive minimum-volatility selection.
		#'
		#' @return A numeric vector of bootstrap estimates, or when
		#'   \code{debug = TRUE}, a diagnostic list.
		#'
		#' @references Bickel, P. J., Gotze, F., and van Zwet, W. R. (1997).
		#'   Resampling fewer than n observations: gains, losses, and remedies for
		#'   losses. \emph{Statistica Sinica}.
		#'
		#' @references Bickel, P. J. and Sakov, A. (2008). On the choice of m in
		#'   the m out of n bootstrap. \emph{The Annals of Statistics}.
		approximate_m_out_of_n_bootstrap_distribution_beta_hat_T_impl = function(B = 501, m = NULL, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL, scaling = "sqrt_n", center = "full_estimate"){
			private$active_resampling_operation = "m_out_of_n_boot"
			on.exit(private$active_resampling_operation <- NULL, add = TRUE)
			private$check_bootstrap_replicate_deadline("m-out-of-n bootstrap setup")
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("m-out-of-n bootstrap inference")
				private$assert_valid_bootstrap_type(bootstrap_type)
				assertCount(B, positive = TRUE)
				assertFlag(show_progress)
				assertFlag(debug)
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = bootstrap_type)
			if (is.list(m)) {
				selection = do.call(self$select_optimal_m_out_of_n_bootstrap, c(list(
					B = as.integer(m$B %||% min(as.integer(B), 251L)),
					alpha = as.numeric(m$alpha %||% 0.05),
					bootstrap_type = bootstrap_type,
					scaling = scaling,
					show_progress = show_progress
				), m[setdiff(names(m), c("method", "B", "alpha"))]))
				m = selection$m_optimal
				if (!is.finite(m)) stop("m-out-of-n bootstrap m selection failed.", call. = FALSE)
			}
			m = private$resolve_resampling_size(m, unit_info$n_units, "m")
			cache_key = private$m_out_of_n_bootstrap_cache_key(B = B, m = m, bootstrap_type = bootstrap_type, scaling = scaling, center = center)
			if (!isTRUE(debug)) {
				private$ensure_resampling_distribution_cache("m_out_of_n_boot")
				cached = private$get_cached_resampling_distribution("m_out_of_n_boot", cache_key)
				if (!is.null(cached)) return(cached)
			}
			private$check_bootstrap_replicate_deadline("m-out-of-n bootstrap draw generation")
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
				size = m,
				replace = TRUE,
				stratified = !identical(bootstrap_type, "resample_blocks"),
				preserve_order = FALSE,
				size_label = "m"
			)
			res = private$compute_resampling_draw_distribution(
				draws = draws,
				operation = "m_out_of_n_boot",
				show_progress = show_progress,
				debug = debug,
				cache_key = cache_key
			)
			if (isTRUE(debug)) {
				private$check_bootstrap_replicate_deadline("m-out-of-n bootstrap debug summary")
				est = tryCatch(as.numeric(self$compute_estimate(estimate_only = TRUE))[1L], error = private$resampling_error_to_na)
				res$centered_scaled_values = private$resampling_scaling_factor(m, scaling) * (res$values - est)
				res$full_estimate = est
				res$m = m
				res$n_units = unit_info$n_units
				res$scaling = scaling
				res$unit_type = unit_info$unit_type
			}
			res
		},
		#' @description Computes a centered m-out-of-n bootstrap two-sided p-value.
		#'
		#' @param delta Null treatment effect. Default 0.
		#' @param B Number of resamples. Default 501.
		#' @param m Number of exchangeable units drawn with replacement. If
		#'   \code{NULL} (default), use \code{floor(n_units^0.7)} subject to the
		#'   validation bounds documented for
		#'   \code{approximate_m_out_of_n_bootstrap_distribution_beta_hat_T()}.
		#' @param type P-value type. Currently only \code{"centered"} is supported.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param min_number_usable_samples Minimum number of finite resampled
		#'   estimates required after filtering. Default 5.
		#' @param bootstrap_type Optional empirical-resampling scheme.
		#' @param scaling Scaling sequence for centered m-out-of-n pivots.
		#'
		#' @details The default \code{m = NULL} follows the intermediate-sequence
		#'   convention from the m-out-of-n bootstrap literature:
		#'   \eqn{m \to \infty} and \eqn{m / n \to 0}. The deterministic exponent
		#'   0.7 is a first-pass default; for unstable paths prefer the
		#'   minimum-volatility selector.
		#'
		#' @return A numeric two-sided p-value, or \code{NA_real_} if the path is
		#'   non-estimable.
		compute_m_out_of_n_bootstrap_two_sided_pval_impl = function(delta = 0, B = 501, m = NULL, type = "centered", show_progress = TRUE, min_number_usable_samples = 5L, bootstrap_type = NULL, scaling = "sqrt_n"){
			private$check_bootstrap_replicate_deadline("m-out-of-n bootstrap p-value setup")
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertChoice(type, c("centered"))
			}
			if (is.list(m)) {
				selection = do.call(self$select_optimal_m_out_of_n_bootstrap, c(list(
					B = as.integer(m$B %||% min(as.integer(B), 251L)),
					alpha = as.numeric(m$alpha %||% 0.05),
					bootstrap_type = bootstrap_type,
					scaling = scaling,
					show_progress = show_progress
				), m[setdiff(names(m), c("method", "B", "alpha"))]))
				m = selection$m_optimal
				if (!is.finite(m)) {
					if (isTRUE(private$harden)) private$cache_nonestimable_estimate("m_out_of_n_m_selection_failed")
					return(NA_real_)
				}
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = bootstrap_type)
			m = private$resolve_resampling_size(m, unit_info$n_units, "m")
			pivot = private$m_out_of_n_bootstrap_centered_pivot(
				B = B,
				m = m,
				unit_info = unit_info,
				show_progress = show_progress,
				bootstrap_type = bootstrap_type,
				scaling = scaling
			)
			if (!isTRUE(pivot$ok)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate(pivot$reason)
				return(NA_real_)
			}
			if (length(pivot$finite) < as.integer(min_number_usable_samples)) {
				if (isTRUE(private$harden)) private$cache_nonestimable_estimate("m_out_of_n_too_few_finite_estimates")
				return(NA_real_)
			}
			private$resampling_centered_pval(pivot$centered_scaled, est = pivot$est, delta = delta, n_units = pivot$n_units, scaling = scaling)
		},
		#' @description Computes a basic m-out-of-n bootstrap confidence interval.
		#'
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of resamples. Default 501.
		#' @param m Number of exchangeable units drawn with replacement. If
		#'   \code{NULL} (default), use \code{floor(n_units^0.7)} subject to the
		#'   validation bounds documented for
		#'   \code{approximate_m_out_of_n_bootstrap_distribution_beta_hat_T()}.
		#' @param type Confidence-interval type. Currently only \code{"basic"} is
		#'   supported.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param min_number_usable_samples Minimum number of finite resampled
		#'   estimates required after filtering. Default 5.
		#' @param bootstrap_type Optional empirical-resampling scheme.
		#' @param scaling Scaling sequence for centered m-out-of-n pivots.
		#'
		#' @details The \code{NULL} default is grounded in the standard m-out-of-n
		#'   asymptotic condition \eqn{m \to \infty} and \eqn{m / n \to 0}. The
		#'   minimum-volatility selector is available when a fixed deterministic
		#'   exponent is too brittle for a specific estimator/design path.
		#'
		#' @return A length-2 numeric confidence interval, or
		#'   \code{c(NA_real_, NA_real_)} if the path is non-estimable.
		compute_m_out_of_n_bootstrap_confidence_interval_impl = function(alpha = 0.05, B = 501, m = NULL, type = "basic", show_progress = TRUE, min_number_usable_samples = 5L, bootstrap_type = NULL, scaling = "sqrt_n"){
			private$check_bootstrap_replicate_deadline("m-out-of-n bootstrap CI setup")
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertCount(B, positive = TRUE)
				assertCount(min_number_usable_samples, positive = TRUE)
				assertChoice(type, c("basic"))
			}
			if (is.list(m)) {
				selection = do.call(self$select_optimal_m_out_of_n_bootstrap, c(list(
					B = as.integer(m$B %||% min(as.integer(B), 251L)),
					alpha = alpha,
					bootstrap_type = bootstrap_type,
					scaling = scaling,
					show_progress = show_progress
				), m[setdiff(names(m), c("method", "B", "alpha"))]))
				m = selection$m_optimal
				if (!is.finite(m)) {
					if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, "m_out_of_n_m_selection_failed", stage = "estimate"))
					return(c(NA_real_, NA_real_))
				}
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = bootstrap_type)
			m = private$resolve_resampling_size(m, unit_info$n_units, "m")
			pivot = private$m_out_of_n_bootstrap_centered_pivot(
				B = B,
				m = m,
				unit_info = unit_info,
				show_progress = show_progress,
				bootstrap_type = bootstrap_type,
				scaling = scaling
			)
			if (!isTRUE(pivot$ok)) {
				if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, pivot$reason, stage = "estimate"))
				return(c(NA_real_, NA_real_))
			}
			if (length(pivot$finite) < as.integer(min_number_usable_samples)) {
				if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, "m_out_of_n_too_few_finite_estimates", stage = "estimate"))
				return(c(NA_real_, NA_real_))
			}
			ci = private$resampling_ci_from_centered_distribution(pivot$centered_scaled, alpha = alpha, est = pivot$est, n_units = pivot$n_units, scaling = scaling)
			if (!all(is.finite(ci[1:2])) || private$bootstrap_confidence_interval_extreme(ci, est = pivot$est)) {
				if (isTRUE(private$harden)) return(private$missing_bootstrap_ci(alpha, "m_out_of_n_extreme_confidence_interval", stage = "estimate"))
			}
			ci
		},
		#' @description Selects an m-out-of-n bootstrap size by minimum volatility.
		#'
		#' @param B Number of resamples per candidate size. Default 251.
		#' @param alpha Significance level for interval-width objectives.
		#' @param m_pow_of_n_grid Candidate exponent grid used when
		#'   \code{m_grid = NULL}. Defaults to \code{seq(0.5, 0.9, by = 0.05)}.
		#' @param m_grid Optional explicit integer candidate sizes.
		#' @param objective Selection objective. Currently \code{"ci_width"}.
		#' @param target Target summary. Currently \code{"ci"}.
		#' @param volatility_window Rolling window size used to measure local
		#'   volatility across candidate sizes.
		#' @param bootstrap_type Optional empirical-resampling scheme.
		#' @param scaling Scaling sequence for centered m-out-of-n pivots.
		#' @param show_progress A flag indicating whether a progress bar should be
		#'   displayed.
		#' @param min_finite_fraction Minimum finite-resample fraction required for a
		#'   candidate size to be eligible.
		#'
		#' @details This implements the same minimum-volatility idea used in PTE's
		#'   m-selection workflow: scan admissible intermediate sizes and choose a
		#'   stable region of the target statistic rather than assuming one exponent
		#'   is uniformly optimal.
		#'
		#' @return An \code{EDIMOutOfNBootstrapMSelection} list with the selected
		#'   \code{m}, mapped exponent, candidate table, status, and reason.
		select_optimal_m_out_of_n_bootstrap_impl = function(B = 251, alpha = 0.05, m_pow_of_n_grid = seq(0.5, 0.9, by = 0.05), m_grid = NULL, objective = "ci_width", target = "ci", volatility_window = 3L, bootstrap_type = NULL, scaling = "sqrt_n", show_progress = TRUE, min_finite_fraction = 0.8){
			if (should_run_asserts()) {
				assertCount(B, positive = TRUE)
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertChoice(target, c("ci"))
			}
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = bootstrap_type)
			if (is.null(m_grid)) {
				m_grid_raw = floor(unit_info$n_units ^ m_pow_of_n_grid)
				grid_df = unique(data.frame(m = as.integer(m_grid_raw), pow = as.numeric(m_pow_of_n_grid), stringsAsFactors = FALSE))
				m_grid = grid_df$m
				m_pow_mapped = grid_df$pow
			} else {
				m_grid = as.integer(m_grid)
				m_pow_mapped = rep(NA_real_, length(m_grid))
			}
			evaluator = function(m_size) private$evaluate_m_out_of_n_bootstrap_size(m = m_size, B = B, alpha = alpha, bootstrap_type = bootstrap_type, scaling = scaling, show_progress = FALSE)
			private$select_optimal_resample_size(
				size_grid = m_grid,
				evaluator = evaluator,
				objective = objective,
				volatility_window = volatility_window,
				min_finite_fraction = min_finite_fraction,
				size_name = "m",
				exponent_grid = m_pow_mapped,
				selection_class = "EDIMOutOfNBootstrapMSelection"
			)
		}
		,
		m_out_of_n_bootstrap_cache_key = function(B, m, bootstrap_type = NULL, scaling = "sqrt_n", center = "full_estimate"){
			paste(as.integer(B), as.integer(m), bootstrap_type %||% "NULL", private$resampling_scaling_key(scaling), center, sep = "|")
		},
		m_out_of_n_bootstrap_centered_pivot = function(B, m, unit_info, show_progress = TRUE, bootstrap_type = NULL, scaling = "sqrt_n"){
			private$check_bootstrap_replicate_deadline("m-out-of-n bootstrap centered pivot")
			cache_key = private$m_out_of_n_bootstrap_cache_key(B = B, m = m, bootstrap_type = bootstrap_type, scaling = scaling, center = "full_estimate")
			cached = private$get_cached_centered_resampling_pivot("m_out_of_n_boot", cache_key)
			if (!is.null(cached)) return(cached)
			est = tryCatch(as.numeric(self$compute_estimate(estimate_only = TRUE))[1L], error = private$resampling_error_to_na)
			if (!is.finite(est)) {
				return(list(ok = FALSE, reason = "m_out_of_n_original_estimate_unavailable"))
			}
			boot = self$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T(
				B = B,
				m = m,
				show_progress = show_progress,
				bootstrap_type = bootstrap_type,
				scaling = scaling
			)
			finite = boot[is.finite(boot)]
			pivot = list(
				ok = TRUE,
				reason = NA_character_,
				est = est,
				finite = finite,
				centered_scaled = private$resampling_scaling_factor(m, scaling) * (finite - est),
				m = as.integer(m),
				n_units = unit_info$n_units,
				scaling = scaling
			)
			private$set_cached_centered_resampling_pivot("m_out_of_n_boot", cache_key, pivot)
		},
		resampling_scaling_key = function(scaling){
			if (is.character(scaling)) return(paste(scaling, collapse = ","))
			if (is.list(scaling)) return(paste(names(scaling), unlist(scaling), collapse = ","))
			if (is.function(scaling)) return("<function>")
			as.character(scaling)[1L]
		},
		m_out_of_n_bootstrap_sample_indices = function(m, bootstrap_type = NULL){
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = bootstrap_type)
			selected = private$sample_exchangeable_unit_ids(
				unit_info = unit_info,
				size = as.integer(m),
				replace = TRUE,
				stratified = !identical(bootstrap_type, "resample_blocks")
			)
			private$build_resampling_draw_from_units(unit_info, selected, size_label = "m", preserve_order = FALSE)
		},
		load_m_out_of_n_bootstrap_draw_into_worker = function(worker_state, draw, ...){
			private$load_bootstrap_sample_into_worker(worker_state, draw)
			invisible(worker_state)
		},
		evaluate_m_out_of_n_bootstrap_size = function(m, B, alpha, bootstrap_type = NULL, scaling = "sqrt_n", show_progress = FALSE){
			private$check_bootstrap_replicate_deadline("m-out-of-n bootstrap m selection")
			unit_info = private$get_exchangeable_units(unit = "auto", resampling_type = bootstrap_type)
			m = private$resolve_resampling_size(m, unit_info$n_units, "m")
			est = tryCatch(as.numeric(self$compute_estimate(estimate_only = TRUE))[1L], error = private$resampling_error_to_na)
			if (!is.finite(est)) {
				return(list(status = "nonestimable", dominant_failure_reason = "m_out_of_n_original_estimate_unavailable", finite_fraction = 0, n_finite = 0L))
			}
			boot = self$approximate_m_out_of_n_bootstrap_distribution_beta_hat_T(B = B, m = m, show_progress = show_progress, bootstrap_type = bootstrap_type, scaling = scaling)
			finite = boot[is.finite(boot)]
			centered_scaled = private$resampling_scaling_factor(m, scaling) * (finite - est)
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
					"m_out_of_n_too_few_finite_estimates"
				} else if (!isTRUE(ci_ok)) {
					"m_out_of_n_extreme_confidence_interval"
				} else {
					NA_character_
				}
			)
		}
	)
)
