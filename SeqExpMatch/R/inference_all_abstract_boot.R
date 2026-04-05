#' Bootstrap-based Inference
#'
#' Abstract class for bootstrap-based inference.
#'
#' @keywords internal
InferenceBoot = R6::R6Class("InferenceBoot",
	lock_objects = FALSE,
	inherit = InferenceRandCI,
	public = list(
		#' @description
		#' Creates the bootstrap distribution of the estimate for the treatment effect
		#'
		#' @param B						Number of bootstrap samples. The default is 501.
		#' @param show_progress			A flag indicating whether a progress bar should be displayed.
		#' @param debug					If \code{TRUE}, return a list with the distribution values and
		#'   per-iteration diagnostics including error messages, warning messages, counts of each,
		#'   and summary proportions for iterations with errors, warnings, and illegal (non-finite)
		#'   values. Runs serially. Default \code{FALSE}.
		#' @return 	When \code{debug = FALSE} (default), a numeric vector of length \code{B}
		#'   containing the bootstrap estimates. When \code{debug = TRUE}, a list with: \code{values},
		#'   \code{errors} (list of character vectors, one per iteration), \code{warnings} (list of
		#'   character vectors, one per iteration), \code{num_errors}, \code{num_warnings},
		#'   \code{prop_iterations_with_errors}, \code{prop_iterations_with_warnings}, and
		#'   \code{prop_illegal_values}.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE){
			private$assert_design_supports_resampling("Bootstrap inference")
			assertCount(B, positive = TRUE); assertFlag(debug)

			# Check cache (skipped in debug mode to always get fresh diagnostic results)
			cache_key = as.character(B)
			if (!isTRUE(debug) && !is.null(private$cached_values$boot_distr_cache[[cache_key]])) {
				return(private$cached_values$boot_distr_cache[[cache_key]])
			}

			if (private$verbose) cat("Computing bootstrap distribution...\n")

			# Duplicate objects for thread safety
			inf_template = self$duplicate()
			des_template = private$des_obj$duplicate()
			has_match_structure_local = private$has_match_structure

			run_one_boot_iter = function(worker_des, worker_inf) {
				worker_des$resample_design()
				worker_inf$.__enclos_env__$private$w = worker_des$.__enclos_env__$private$w
				worker_inf$.__enclos_env__$private$y = worker_des$.__enclos_env__$private$y
				if (has_match_structure_local && !is.null(worker_inf$.__enclos_env__$private$compute_basic_match_data)) {
					worker_inf$.__enclos_env__$private$m = worker_des$.__enclos_env__$private$m
					worker_inf$.__enclos_env__$private$compute_basic_match_data()
				}
				worker_inf$compute_treatment_estimate(estimate_only = TRUE)
			}

			if (isTRUE(debug)) {
				debug_results = vector("list", B)
				for (idx in seq_len(B)) {
					iter_warns = character(0)
					iter_errs = character(0)
					iter_val = withCallingHandlers(
						tryCatch({
							worker_des = des_template$duplicate()
							worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
							run_one_boot_iter(worker_des, worker_inf)
						}, error = function(e) { iter_errs <<- c(iter_errs, conditionMessage(e)); NA_real_ }),
						warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
					)
					debug_results[[idx]] = list(
						val = as.numeric(iter_val)[1L],
						errors = iter_errs,
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

			# Determine cores — warm-up guard: run one iteration serially to estimate per-iteration
			# cost, then only parallelize if computation outweighs overhead per worker.
			# For a persistent fork cluster the per-call overhead is ~10ms (socket round-trip);
			# for mclapply (per-call fork) it is ~500ms.
			# We run the warmup iteration TWICE and use the second timing. The first call often
			# pays cold-start penalties (C++ JIT, OS page-cache misses, R bytecode compilation)
			# that can inflate the estimate 5–15× vs steady-state cost, causing the guard to
			# wrongly choose parallel for small B values like r = 19.
			actual_cores = private$effective_parallel_cores("bootstrap", self$num_cores)
			if (actual_cores > 1L) {
				do_warmup_iter = function() {
					w_des = des_template$duplicate()
					w_inf = inf_template$duplicate(make_fork_cluster = FALSE)
					tryCatch(run_one_boot_iter(w_des, w_inf), error = function(e) NA_real_)
				}
				system.time(do_warmup_iter())  # First call: discarded (cold-start overhead)
				t_boot_warmup = system.time(do_warmup_iter())[[3]]  # Second call: representative cost
				# For fork cluster: round-trip per task ~10ms, but first call also pays ~300ms
				# cluster-creation cost. For mclapply: per-fork cost ~500ms per worker.
				fork_overhead_estimate = if (!is.null(get_global_fork_cluster())) 0.01 else 0.5
				cluster_create_overhead = if (!is.null(get_global_fork_cluster()) && is.null(get_global_fork_cluster())) 0.3 else 0.0
				if (!(t_boot_warmup * B > fork_overhead_estimate * actual_cores + cluster_create_overhead))
					actual_cores = 1L
			}

			# Use private$par_lapply which handles both persistent fork clusters and other strategies.
			boot_distr = unlist(private$par_lapply(1:B, function(idx) {
				worker_des = des_template$duplicate()
				worker_inf = inf_template$duplicate(make_fork_cluster = FALSE)
				tryCatch(run_one_boot_iter(worker_des, worker_inf), error = function(e) NA_real_)
			}, n_cores = actual_cores, show_progress = show_progress,
			export_list = list(
				des_template = des_template,
				inf_template = inf_template,
				has_match_structure_local = has_match_structure_local,
				run_one_boot_iter = run_one_boot_iter
			)))

			if (!is.numeric(boot_distr)) boot_distr = as.numeric(boot_distr)

			if (is.null(private$cached_values$boot_distr_cache)) private$cached_values$boot_distr_cache = list()
			private$cached_values$boot_distr_cache[[cache_key]] = boot_distr
			boot_distr
		},

		#' @description
		#' Computes a bootstrap-based two-sided p-value for the treatment effect.
		#'
		#' @param delta					Null hypothesis value. Default 0.
		#' @param B						Number of bootstrap samples. Default 501.
		#' @param na.rm					Remove non-finite bootstrap replicates. Default FALSE.
		#'
		#' @return 	A bootstrap two-sided p-value.
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE){
			assertNumeric(delta, len = 1)
			assertCount(B, positive = TRUE)
			assertFlag(na.rm)
			boot_distr = self$approximate_bootstrap_distribution_beta_hat_T(B)
			if (isTRUE(na.rm)) boot_distr = boot_distr[is.finite(boot_distr)]
			else if (any(!is.finite(boot_distr))) return(NA_real_)
			if (length(boot_distr) == 0L) return(NA_real_)
			est = as.numeric(self$compute_treatment_estimate())
			if (length(est) == 0L || !is.finite(est[1])) return(NA_real_)
			est = est[1]
			# Shift bootstrap distribution to be centered at delta (null hypothesis)
			boot_null = boot_distr - mean(boot_distr) + delta
			n_bs = length(boot_null)
			min(1, max(2 / n_bs, 2 * min(
				sum(boot_null >= est) / n_bs,
				sum(boot_null <= est) / n_bs
			)))
		},

		#' @description
		#' Computes a bootstrap-based confidence interval.
		#'
		#' @param alpha					The confidence level 1 - \code{alpha}. Default 0.05.
		#' @param B						Number of bootstrap samples. Default 501.
		#' @param type					Bootstrap CI type. Supported values are
		#'   \code{"percentile"}, \code{"basic"}, \code{"studentized"},
		#'   \code{"bootstrap-t"}, \code{"symmetric-percentile-t"},
		#'   \code{"bca"}, \code{"prepivoted"}, \code{"double-bootstrap"},
		#'   \code{"calibrated"}, and \code{"smoothed"}.
		#' @param na.rm                                   Remove non-finite bootstrap replicates.
		#'   Default TRUE. Non-finite replicates are always removed internally.
		#' @param show_progress			Show progress bar.
		#'
		#' @return 	A bootstrap confidence interval.
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = "percentile", na.rm = TRUE, show_progress = TRUE){
			private$assert_design_supports_resampling("Bootstrap inference")
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			type = tolower(type)
			assertChoice(type, c(
				"percentile", "basic", "studentized", "bootstrap-t",
				"symmetric-percentile-t", "bca", "prepivoted",
				"double-bootstrap", "calibrated", "smoothed"
			))

			est = as.numeric(self$compute_treatment_estimate(estimate_only = FALSE))
			if (length(est) == 0L || !is.finite(est[1])) stop("Bootstrap confidence interval returned NA bounds")
			est = est[1]

			boot_stats = if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t", "bca", "prepivoted", "double-bootstrap", "calibrated", "smoothed")) {
				private$approximate_bootstrap_statistics_beta_hat_T(B = B, show_progress = show_progress, na.rm = na.rm, smooth = identical(type, "smoothed"))
			} else {
				list(theta = private$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = show_progress), se = NULL)
			}

			boot_distr = boot_stats$theta
			boot_distr = boot_distr[is.finite(boot_distr)]
			if (length(boot_distr) < max(10L, B / 2)) stop("Bootstrap confidence interval returned NA bounds")

			if (type == "percentile") {
				ci = private$ci_from_boot_distribution(boot_distr, alpha, "percentile")
			} else if (type == "basic") {
				ci = private$ci_from_boot_distribution(boot_distr, alpha, "basic", est = est)
			} else if (type %in% c("studentized", "bootstrap-t")) {
				ci = private$ci_studentized(boot_stats, alpha, est)
			} else if (type == "symmetric-percentile-t") {
				ci = private$ci_symmetric_studentized(boot_stats, alpha, est)
			} else if (type == "bca") {
				ci = private$ci_bca(boot_distr, alpha, est)
			} else if (type %in% c("prepivoted", "double-bootstrap", "calibrated")) {
				ci = private$ci_calibrated_bootstrap(alpha, B, type, est, show_progress = show_progress, na.rm = na.rm)
			} else if (type == "smoothed") {
				ci = private$ci_smoothed_bootstrap(alpha, B, est, show_progress = show_progress, na.rm = na.rm)
			} else {
				stop("Unsupported bootstrap CI type: ", type)
			}

			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		}
	),
	private = list(
		# Cache for bootstrap distributions
		boot_distr_cache = list(),

		bootstrap_sample_indices = function(n){
			sample.int(n, n, replace = TRUE)
		},

		bootstrap_subset_inference = function(indices, smooth = FALSE){
			indices = as.integer(indices)
			if (length(indices) == 0L) return(NULL)

			orig_des = private$des_obj
			orig_des_priv = private$des_obj_priv_int
			sub_des = orig_des$duplicate(verbose = FALSE)
			sub_des_priv = sub_des$.__enclos_env__$private

			subset_field = function(x){
				if (is.null(x)) return(NULL)
				if (is.data.frame(x) || is.matrix(x)) return(x[indices, , drop = FALSE])
				if (is.list(x) && !is.data.frame(x)) return(x[indices])
				if (is.atomic(x) && length(x) >= max(indices)) return(x[indices])
				x
			}

			if (!is.null(orig_des_priv$Xraw)) sub_des_priv$Xraw = subset_field(orig_des_priv$Xraw)
			if (!is.null(orig_des_priv$Ximp)) sub_des_priv$Ximp = subset_field(orig_des_priv$Ximp)
			if (!is.null(orig_des_priv$X)) sub_des_priv$X = subset_field(orig_des_priv$X)
			if (!is.null(orig_des_priv$w)) sub_des_priv$w = as.numeric(orig_des_priv$w[indices])
			if (!is.null(orig_des_priv$y)) sub_des_priv$y = as.numeric(orig_des_priv$y[indices])
			if (!is.null(orig_des_priv$dead)) sub_des_priv$dead = as.numeric(orig_des_priv$dead[indices])
			if (!is.null(orig_des_priv$m)) sub_des_priv$m = as.integer(orig_des_priv$m[indices])
			if (!is.null(orig_des_priv$y_i_t_i)) sub_des_priv$y_i_t_i = orig_des_priv$y_i_t_i[indices]
			sub_des_priv$all_subject_data_cache = list()
			sub_des_priv$p_raw_t = if (!is.null(sub_des_priv$Xraw)) ncol(sub_des_priv$Xraw) else NULL
			sub_des_priv$t = length(indices)
			sub_des_priv$n = length(indices)
			sub_des_priv$fixed_sample = TRUE
			if (smooth && !is.null(sub_des_priv$y) && private$des_obj_priv_int$response_type == "continuous") {
				sd_y = stats::sd(as.numeric(sub_des_priv$y), na.rm = TRUE)
				if (is.finite(sd_y) && sd_y > 0) {
					sub_des_priv$y = as.numeric(sub_des_priv$y) + stats::rnorm(length(sub_des_priv$y), 0, sd_y / sqrt(max(1, length(indices))))
				}
			}

			sub_inf = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			sub_inf_priv = sub_inf$.__enclos_env__$private
			sub_inf_priv$des_obj = sub_des
			sub_inf_priv$des_obj_priv_int = sub_des_priv
			sub_inf_priv$y = sub_des_priv$y
			sub_inf_priv$y_temp = sub_des_priv$y
			sub_inf_priv$w = sub_des_priv$w
			sub_inf_priv$dead = sub_des_priv$dead
			sub_inf_priv$n = length(indices)
			sub_inf_priv$X = NULL
			sub_inf_priv$cached_values = list()
			sub_inf_priv$cached_values$rand_distr_cache = list()
			sub_inf_priv$cached_values$permutations_cache = list()
			sub_inf_priv$cached_values$m_cache = list()
			if (!is.null(sub_des_priv$m)) sub_inf_priv$m = sub_des_priv$m
			sub_inf
		},

		bootstrap_replication_stats = function(indices, smooth = FALSE){
			sub_inf = private$bootstrap_subset_inference(indices, smooth = smooth)
			if (is.null(sub_inf)) return(c(theta = NA_real_, se = NA_real_))
			tryCatch({
				sub_inf$compute_treatment_estimate(estimate_only = FALSE)
				theta = as.numeric(sub_inf$.__enclos_env__$private$cached_values$beta_hat_T)[1]
				se = as.numeric(sub_inf$.__enclos_env__$private$cached_values$s_beta_hat_T)[1]
				if (!is.finite(theta)) theta = NA_real_
				if (!is.finite(se)) se = NA_real_
				c(theta = theta, se = se)
			}, error = function(e) c(theta = NA_real_, se = NA_real_))
		},

		approximate_bootstrap_statistics_beta_hat_T = function(B = 501, show_progress = TRUE, na.rm = TRUE, smooth = FALSE){
			assertCount(B, positive = TRUE)
			n = private$des_obj$get_n()
			stats_mat = matrix(NA_real_, nrow = B, ncol = 2L)
			pb = NULL
			if (isTRUE(show_progress) && B > 1L) {
				pb = utils::txtProgressBar(min = 0, max = B, style = 3)
				on.exit(try(close(pb), silent = TRUE), add = TRUE)
			}
			for (b in seq_len(B)) {
				idx = private$bootstrap_sample_indices(n)
				stats_mat[b, ] = private$bootstrap_replication_stats(idx, smooth = smooth)
				if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
			}
			if (isTRUE(na.rm)) {
				ok = is.finite(stats_mat[, 1L]) & is.finite(stats_mat[, 2L])
				stats_mat = stats_mat[ok, , drop = FALSE]
			}
			if (nrow(stats_mat) == 0L) {
				return(list(theta = numeric(0), se = numeric(0)))
			}
			list(theta = stats_mat[, 1L], se = stats_mat[, 2L])
		},

		approximate_jackknife_distribution_beta_hat_T = function(){
			n = private$des_obj$get_n()
			if (n <= 1L) return(numeric(0))
			jack = rep(NA_real_, n)
			for (i in seq_len(n)) {
				idx = seq_len(n)[-i]
				jack[i] = private$bootstrap_replication_stats(idx, smooth = FALSE)[["theta"]]
			}
			jack
		},

		ci_from_boot_distribution = function(boot_distr, alpha, type, est = NULL){
			type = tolower(type)
			if (length(boot_distr) == 0L) stop("Bootstrap confidence interval returned NA bounds")
			if (is.null(est)) est = as.numeric(self$compute_treatment_estimate())[1]
			if (type == "percentile") {
				stats::quantile(boot_distr, probs = c(alpha / 2, 1 - alpha / 2), names = FALSE, type = 8)
			} else {
				2 * est - stats::quantile(boot_distr, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
			}
		},

		ci_studentized = function(boot_stats, alpha, est){
			se_hat = private$infer_original_se()
			if (!is.finite(se_hat) || se_hat <= 0) stop("Studentized bootstrap CI requires a finite standard error.")
			t_vals = (boot_stats$theta - est) / boot_stats$se
			t_vals = t_vals[is.finite(t_vals)]
			if (length(t_vals) < 10L) stop("Studentized bootstrap CI returned too few finite bootstrap draws.")
			q = stats::quantile(t_vals, probs = c(1 - alpha / 2, alpha / 2), names = FALSE, type = 8)
			c(est - q[1L] * se_hat, est - q[2L] * se_hat)
		},

		ci_symmetric_studentized = function(boot_stats, alpha, est){
			se_hat = private$infer_original_se()
			if (!is.finite(se_hat) || se_hat <= 0) stop("Symmetric percentile-t bootstrap CI requires a finite standard error.")
			t_vals = (boot_stats$theta - est) / boot_stats$se
			t_vals = abs(t_vals[is.finite(t_vals)])
			if (length(t_vals) < 10L) stop("Symmetric percentile-t bootstrap CI returned too few finite bootstrap draws.")
			q = stats::quantile(t_vals, probs = 1 - alpha, names = FALSE, type = 8)
			c(est - q[1L] * se_hat, est + q[1L] * se_hat)
		},

		ci_bca = function(boot_distr, alpha, est){
			jack = private$approximate_jackknife_distribution_beta_hat_T()
			jack = jack[is.finite(jack)]
			if (length(jack) < 2L) stop("BCa interval requires jackknife estimates.")
			z0 = stats::qnorm(mean(boot_distr < est))
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			alpha_vec = c(alpha / 2, 1 - alpha / 2)
			z_alpha = stats::qnorm(alpha_vec)
			adj = stats::pnorm(z0 + (z0 + z_alpha) / pmax(.Machine$double.eps, 1 - a * (z0 + z_alpha)))
			adj = pmin(1, pmax(0, adj))
			stats::quantile(boot_distr, probs = adj, names = FALSE, type = 8)
		},

		ci_calibrated_bootstrap = function(alpha, B, type, est, show_progress = TRUE, na.rm = TRUE){
			n_outer = max(25L, min(as.integer(B), 101L))
			n_inner = max(25L, min(as.integer(B), 51L))
			alpha_grid = unique(pmin(0.49, pmax(.Machine$double.eps, c(alpha / 4, alpha / 2, alpha, min(0.25, alpha * 1.5), min(0.4, alpha * 2)))))
			coverage = rep(NA_real_, length(alpha_grid))
			for (j in seq_along(alpha_grid)) {
				a = alpha_grid[j]
				covered = logical(n_outer)
				for (b in seq_len(n_outer)) {
					idx = private$bootstrap_sample_indices(private$des_obj$get_n())
					outer_inf = private$bootstrap_subset_inference(idx, smooth = identical(type, "smoothed"))
					if (is.null(outer_inf)) next
					inner_boot = outer_inf$approximate_bootstrap_distribution_beta_hat_T(B = n_inner, show_progress = FALSE, debug = FALSE)
					inner_boot = inner_boot[is.finite(inner_boot)]
					if (length(inner_boot) < 10L) next
					ci_inner = private$ci_from_boot_distribution(inner_boot, a, "percentile")
					covered[b] = is.finite(ci_inner[1]) && is.finite(ci_inner[2]) && est >= ci_inner[1] && est <= ci_inner[2]
				}
				coverage[j] = mean(covered, na.rm = TRUE)
			}
			target_coverage = 1 - alpha
			idx_best = which.min(abs(coverage - target_coverage))
			alpha_adj = alpha_grid[idx_best]
			outer_boot = private$approximate_bootstrap_statistics_beta_hat_T(
				B = B,
				show_progress = show_progress,
				na.rm = na.rm,
				smooth = identical(type, "smoothed")
			)$theta
			outer_boot = outer_boot[is.finite(outer_boot)]
			if (length(outer_boot) < 10L) stop("Calibrated bootstrap CI returned too few finite bootstrap draws.")
			ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "percentile")
			if (identical(type, "double-bootstrap")) {
				ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "basic", est = est)
			}
			if (identical(type, "prepivoted")) {
				ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "percentile")
			}
			if (identical(type, "calibrated")) {
				ci = private$ci_from_boot_distribution(outer_boot, alpha_adj, "percentile")
			}
			ci
		},

		ci_smoothed_bootstrap = function(alpha, B, est, show_progress = TRUE, na.rm = TRUE){
			boot_stats = private$approximate_bootstrap_statistics_beta_hat_T(B = B, show_progress = show_progress, na.rm = na.rm, smooth = TRUE)
			boot_distr = boot_stats$theta[is.finite(boot_stats$theta)]
			if (length(boot_distr) < 10L) stop("Smoothed bootstrap CI returned too few finite bootstrap draws.")
			private$ci_from_boot_distribution(boot_distr, alpha, "percentile", est = est)
		},

		infer_original_se = function(){
			fresh = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			fresh$.__enclos_env__$private$cached_values = list()
			tryCatch({
				fresh$compute_treatment_estimate(estimate_only = FALSE)
				as.numeric(fresh$.__enclos_env__$private$cached_values$s_beta_hat_T)[1]
			}, error = function(e) NA_real_)
		}
	)
)
