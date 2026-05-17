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
		#' @description Recomputes the treatment estimate under Bayesian-bootstrap
		#'   subject-, block-, cluster-, or matched-set weights.
		#'
		#' This is an abstract hook implemented by concrete inference families that
		#' support weighted re-estimation.
		#'
		#' @param subject_or_block_weights Numeric Bayesian-bootstrap weights at the
		#'   design's exchangeable resampling unit. For ordinary designs these are
		#'   subject-level weights. For blocking, clustering, or matching designs
		#'   these may instead be block-, cluster-, pair-, or matched-set-level
		#'   weights, depending on \code{weighting_unit_type}.
		#' @param estimate_only If \code{TRUE}, compute only the point estimate for
		#'   the weighted replicate.
		#'
		#' @return A numeric treatment-effect estimate for the weighted replicate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			stop(class(self)[1], " must implement weighted bootstrap estimation.")
		},
		#' @description Creates the Bayesian-bootstrap distribution of the treatment
		#'   estimate using Dirichlet weights.
		#'
		#' @param B Number of Bayesian-bootstrap replicates. The default is 501.
		#' @param show_progress A flag indicating whether a progress bar should be displayed.
		#' @param debug If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics including error messages, warning messages,
		#'   counts of each, and summary proportions for iterations with errors,
		#'   warnings, and illegal (non-finite) values. Runs serially. Default
		#'   \code{FALSE}.
		#' @param weighting_unit_type Optional Bayesian-bootstrap weighting-unit
		#'   scheme. Legal public values are:
		#'   \describe{
		#'     \item{\code{NULL}}{Use the design's default weighting-unit logic. For
		#'       ordinary non-blocking designs this is the usual subject-level
		#'       Bayesian bootstrap. For certain blocking designs, \code{NULL} maps
		#'       to the same behavior as \code{"within_blocks"}.}
		#'     \item{\code{"within_blocks"}}{Only legal for blocking-style designs
		#'       that support block-aware weighting:
		#'       \code{DesignFixedBlocking}, \code{DesignFixedOptimalBlocks},
		#'       \code{DesignSeqOneByOneSPBR}, and
		#'       \code{DesignFixedBlockedCluster}. Draws Dirichlet weights on
		#'       observational units within each observed block/stratum. For blocked
		#'       cluster designs this means cluster-within-stratum weights.}
		#'     \item{\code{"resample_blocks"}}{Only legal for the same
		#'       blocking-style designs as \code{"within_blocks"}. Draws Dirichlet
		#'       weights on whole observed blocks/strata rather than on units within
		#'       each block.}
		#'   }
		#'   Any non-\code{NULL} value is rejected for designs outside that blocking
		#'   family.
		#'
		#' @return When \code{debug = FALSE} (default), a numeric vector of length
		#'   \code{B} containing the Bayesian-bootstrap estimates. When
		#'   \code{debug = TRUE}, a list with: \code{values}, \code{errors} (list of
		#'   character vectors, one per iteration), \code{warnings} (list of
		#'   character vectors, one per iteration), \code{num_errors},
		#'   \code{num_warnings}, \code{prop_iterations_with_errors},
		#'   \code{prop_iterations_with_warnings}, and
		#'   \code{prop_illegal_values}.
		approximate_bayesian_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, weighting_unit_type = NULL){
			if (should_run_asserts()) {
				private$assert_design_supports_resampling("Bayesian bootstrap inference")
				private$assert_valid_bootstrap_type(weighting_unit_type)
				assertCount(B, positive = TRUE)
				assertFlag(debug)
			}
			cache_key = private$bayesian_bootstrap_cache_key(B = B, weighting_unit_type = weighting_unit_type)
			if (!isTRUE(debug) && !is.null(private$cached_values$bayes_boot_distr_cache[[cache_key]])) {
				return(private$cached_values$bayes_boot_distr_cache[[cache_key]])
			}
			inf_template = self$duplicate()
			draws = replicate(
				as.integer(B),
				private$bayesian_bootstrap_sample_weights(weighting_unit_type = weighting_unit_type),
				simplify = FALSE
			)
			run_one_iter = function(worker_inf) {
				draw = draws[[1L]]
				worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
				as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
					subject_or_block_weights = draw$subject_or_block_weights,
					estimate_only = TRUE
				))[1L]
			}
			if (isTRUE(debug)) {
				run_debug_iter = function(draw, worker_inf = NULL, worker_state = NULL) {
					iter_warns = character(0)
					iter_errs = character(0)
					iter_val = withCallingHandlers(
						tryCatch({
							if (!is.null(worker_state)) {
								private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
								private$compute_bayesian_bootstrap_worker_estimate(worker_state)
							} else {
								worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
								as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
									subject_or_block_weights = draw$subject_or_block_weights,
									estimate_only = TRUE
								))[1L]
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
						lapply(idxs, function(idx) run_debug_iter(draw = draws[[idx]], worker_state = worker_state))
					}
				} else {
					function(idxs) {
						lapply(idxs, function(idx) {
							worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
							run_debug_iter(draw = draws[[idx]], worker_inf = worker_inf)
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
						draw = draws[[1L]]
						private$load_bayesian_bootstrap_weights_into_worker(worker_state, draw)
						tryCatch(private$compute_bayesian_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
					}
				} else {
					function() {
						worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
						draw = draws[[1L]]
						tryCatch({
							worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
							as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
								subject_or_block_weights = draw$subject_or_block_weights,
								estimate_only = TRUE
							))[1L]
						}, error = function(e) NA_real_)
					}
				}
				system.time(do_warmup_iter())
				t_boot_warmup = system.time(do_warmup_iter())[[3]]
				fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.5
				if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores)) actual_cores = 1L
			}
			boot_distr = if (isTRUE(private$use_reusable_bootstrap_worker())) {
				private$compute_bayesian_bootstrap_distribution_with_reused_workers(
					draws = draws,
					actual_cores = actual_cores,
					show_progress = show_progress
				)
			} else {
				unlist(private$par_lapply(
					seq_along(draws),
					function(idx) {
						worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
						draw = draws[[idx]]
						tryCatch({
							worker_inf$.__enclos_env__$private$current_bayesian_bootstrap_context = draw$context
							as.numeric(worker_inf$compute_estimate_with_bootstrap_weights(
								subject_or_block_weights = draw$subject_or_block_weights,
								estimate_only = TRUE
							))[1L]
						}, error = function(e) NA_real_)
					},
					n_cores = actual_cores,
					show_progress = show_progress,
					export_list = list(inf_template = inf_template, draws = draws)
				))
			}
			boot_distr = as.numeric(boot_distr)
			if (is.null(private$cached_values$bayes_boot_distr_cache)) private$cached_values$bayes_boot_distr_cache = list()
			private$cached_values$bayes_boot_distr_cache[[cache_key]] = boot_distr
			boot_distr
		},
		#' @description Computes a Bayesian-bootstrap-based two-sided p-value for
		#'   the treatment effect.
		#'
		#' @param delta Null hypothesis value. Default 0.
		#' @param B Number of Bayesian-bootstrap replicates. Default 501.
		#' @param type Type of Bayesian-bootstrap p-value. Supported values are
		#'   \code{"percentile"}, \code{"symmetric"}, and \code{"wald"}.
		#' @param na.rm If \code{TRUE}, discard non-finite bootstrap replicates before
		#'   computing the p-value. Otherwise, any non-finite replicate returns
		#'   \code{NA}.
		#' @param show_progress A flag indicating whether a progress bar should be displayed.
		#' @param min_number_usable_samples Minimum number of finite Bayesian-bootstrap
		#'   replicates required after filtering. Default 5.
		#' @param weighting_unit_type Optional Bayesian-bootstrap weighting-unit
		#'   scheme. See
		#'   \code{\link{InferenceBayesianBootstrap$approximate_bayesian_bootstrap_distribution_beta_hat_T}()}.
		#'
		#' @return A numeric two-sided p-value, or \code{NA_real_} if too few usable
		#'   replicates remain or the estimate is non-finite.
		compute_bayesian_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = FALSE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
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
				weighting_unit_type = weighting_unit_type
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
		#' @description Computes a Bayesian-bootstrap confidence interval for the
		#'   treatment effect.
		#'
		#' @param alpha Significance level. Default 0.05.
		#' @param B Number of Bayesian-bootstrap replicates. Default 501.
		#' @param type Type of Bayesian-bootstrap interval. Supported values are
		#'   \code{"percentile"}, \code{"basic"}, and \code{"wald"}.
		#' @param na.rm If \code{TRUE}, discard non-finite bootstrap replicates before
		#'   constructing the interval.
		#' @param show_progress A flag indicating whether a progress bar should be displayed.
		#' @param min_number_usable_samples Minimum number of finite Bayesian-bootstrap
		#'   replicates required after filtering. Default 5.
		#' @param weighting_unit_type Optional Bayesian-bootstrap weighting-unit
		#'   scheme. See
		#'   \code{\link{InferenceBayesianBootstrap$approximate_bayesian_bootstrap_distribution_beta_hat_T}()}.
		#'
		#' @return A length-2 numeric confidence interval. Returns
		#'   \code{c(NA_real_, NA_real_)} when the estimate is non-finite or too few
		#'   usable replicates remain.
		compute_bayesian_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL){
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
				weighting_unit_type = weighting_unit_type
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
		current_bayesian_bootstrap_subject_or_block_weights = NULL,
		bayesian_bootstrap_cache_key = function(B, weighting_unit_type = NULL){
			paste0(as.integer(B), "::", weighting_unit_type %||% "default")
		},
		build_bayesian_bootstrap_context = function(weighting_unit_type = NULL){
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
				if (identical(weighting_unit_type, "resample_blocks")) {
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
		bayesian_bootstrap_sample_weights = function(weighting_unit_type = NULL){
			ctx = private$build_bayesian_bootstrap_context(weighting_unit_type = weighting_unit_type)
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
			worker_priv = worker_state$worker$.__enclos_env__$private
			worker_priv$current_bayesian_bootstrap_subject_or_block_weights = as.numeric(draw$subject_or_block_weights)
			worker_priv$current_bayesian_bootstrap_context = draw$context
		},
		compute_bayesian_bootstrap_worker_estimate = function(worker_state){
			worker_priv = worker_state$worker$.__enclos_env__$private
			theta = as.numeric(worker_state$worker$compute_estimate_with_bootstrap_weights(
				subject_or_block_weights = worker_priv$current_bayesian_bootstrap_subject_or_block_weights,
				estimate_only = TRUE
			))[1L]
			if (is.function(worker_state$worker$is_nonestimable) &&
			    isTRUE(worker_state$worker$is_nonestimable("estimate"))) {
				return(NA_real_)
			}
			theta
		},
		compute_bayesian_bootstrap_distribution_with_reused_workers = function(draws, actual_cores, show_progress = FALSE){
			B = length(draws)
			chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
			chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
			chunks = split(seq_len(B), chunk_id)
			run_chunk = function(idxs) {
				worker_state = private$create_bootstrap_worker_state()
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					draw = draws[[idxs[[k]]]]
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
