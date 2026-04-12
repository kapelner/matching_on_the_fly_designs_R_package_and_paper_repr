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
		#' @param bootstrap_type      The type of bootstrap resampling to perform.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			private$assert_design_supports_resampling("Bootstrap inference")
			private$assert_valid_bootstrap_type(bootstrap_type)
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
				worker_des$.__enclos_env__$private$resample_assignment()
				worker_inf$.__enclos_env__$private$cached_values = list()
				worker_inf$.__enclos_env__$private$w = worker_des$.__enclos_env__$private$w
				worker_inf$.__enclos_env__$private$y = worker_des$.__enclos_env__$private$y
				worker_inf$.__enclos_env__$private$y_temp = worker_des$.__enclos_env__$private$y
				worker_inf$.__enclos_env__$private$dead = worker_des$.__enclos_env__$private$dead
				if (has_match_structure_local && !is.null(worker_inf$.__enclos_env__$private$compute_basic_match_data)) {
					worker_inf$.__enclos_env__$private$m = worker_des$.__enclos_env__$private$m
					worker_inf$.__enclos_env__$private$compute_basic_match_data()
				}
				worker_inf$compute_treatment_estimate(estimate_only = TRUE)
			}

			if (isTRUE(debug)) {
				debug_results = vector("list", B)
				if (isTRUE(private$supports_reusable_bootstrap_worker())) {
					worker_state = private$create_bootstrap_worker_state()
					for (idx in seq_len(B)) {
						iter_warns = character(0)
						iter_errs = character(0)
						iter_val = withCallingHandlers(
							tryCatch({
								boot_draw = private$bootstrap_sample_indices(private$n, bootstrap_type)
								private$load_bootstrap_sample_into_worker(worker_state, boot_draw$i_b)
								private$compute_bootstrap_worker_estimate(worker_state)
							}, error = function(e) { iter_errs <<- c(iter_errs, conditionMessage(e)); NA_real_ }),
							warning = function(w) { iter_warns <<- c(iter_warns, conditionMessage(w)); invokeRestart("muffleWarning") }
						)
						debug_results[[idx]] = list(
							val = as.numeric(iter_val)[1L],
							errors = iter_errs,
							warnings = iter_warns
						)
					}
				} else {
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
				}
				values = sapply(debug_results, `[[`, "val")
				errors_list = lapply(debug_results, `[[`, "errors")
				warnings_list = lapply(debug_results, `[[`, "warnings")
				num_errors_vec = lengths(errors_list)
				num_warnings_vec = lengths(warnings_list)
				# Populate the normal cache so subsequent non-debug calls can reuse the values
				if (is.null(private$cached_values$boot_distr_cache)) private$cached_values$boot_distr_cache = list()
				private$cached_values$boot_distr_cache[[cache_key]] = values
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
				do_warmup_iter = if (isTRUE(private$supports_reusable_bootstrap_worker())) {
					function() {
						worker_state = private$create_bootstrap_worker_state()
						boot_draw = private$bootstrap_sample_indices(private$n, bootstrap_type)
						private$load_bootstrap_sample_into_worker(worker_state, boot_draw$i_b)
						tryCatch(private$compute_bootstrap_worker_estimate(worker_state), error = function(e) NA_real_)
					}
				} else {
					function() {
						w_des = des_template$duplicate()
						w_inf = inf_template$duplicate(make_fork_cluster = FALSE)
						tryCatch(run_one_boot_iter(w_des, w_inf), error = function(e) NA_real_)
					}
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

			boot_distr = if (isTRUE(private$supports_reusable_bootstrap_worker())) {
				private$compute_bootstrap_distribution_with_reused_workers(
					B = B,
					actual_cores = actual_cores,
					show_progress = show_progress,
					bootstrap_type = bootstrap_type
				)
			} else {
				# Use private$par_lapply which handles both persistent fork clusters and other strategies.
				unlist(private$par_lapply(1:B, function(idx) {
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
			}

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
		#' @param type					Bootstrap p-value type. Supported values are
		#'   \code{"percentile"} (default), \code{"symmetric"}, \code{"studentized"},
		#'   \code{"bootstrap-t"}, and \code{"bca"}.
		#'   \code{"percentile"}: shifts the bootstrap distribution to be centred at
		#'   \code{delta} and counts the two-tail proportion (Hall 1992).
		#'   \code{"symmetric"}: uses \eqn{|T^* - \bar{T}^*| \ge |t_{\rm obs} - \delta|}
		#'   for a symmetric one-sample test; recommended by Hall & Wilson (1991) when the
		#'   null distribution may be skewed.
		#'   \code{"studentized"} / \code{"bootstrap-t"}: pivots by the per-replicate
		#'   standard error, giving O(n^{-1}) error versus O(n^{-1/2}) for the percentile
		#'   method (Hall 1992; Davidson & MacKinnon 1999).
		#'   \code{"bca"}: bias-corrected and accelerated p-value via closed-form CI
		#'   inversion using the jackknife acceleration and bias-correction constants;
		#'   second-order accurate (Efron 1987; Efron & Tibshirani 1993).
		#' @param na.rm					Remove non-finite bootstrap replicates. Default FALSE.
		#'
		#' @return 	A bootstrap two-sided p-value.
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = FALSE){
			assertNumeric(delta, len = 1)
			assertCount(B, positive = TRUE)
			assertFlag(na.rm)
			type = tolower(private$get_bootstrap_type(type))
			assertChoice(type, c("percentile", "symmetric", "studentized", "bootstrap-t", "bca"))

			need_se = type %in% c("studentized", "bootstrap-t")
			if (need_se) {
				boot_stats = private$approximate_bootstrap_statistics_beta_hat_T(B = B, na.rm = isTRUE(na.rm))
				boot_distr = boot_stats$theta
			} else {
				boot_distr = self$approximate_bootstrap_distribution_beta_hat_T(B)
				boot_stats = NULL
			}

			if (isTRUE(na.rm)) boot_distr = boot_distr[is.finite(boot_distr)]
			else if (any(!is.finite(boot_distr))) return(NA_real_)
			if (length(boot_distr) == 0L) return(NA_real_)

			est = as.numeric(self$compute_treatment_estimate())
			if (length(est) == 0L || !is.finite(est[1])) return(NA_real_)
			est = est[1]

			if (type == "percentile") {
				# Shift bootstrap distribution to be centred at delta (null hypothesis)
				boot_null = boot_distr - mean(boot_distr) + delta
				n_bs = length(boot_null)
				min(1, max(2 / n_bs, 2 * min(
					sum(boot_null >= est) / n_bs,
					sum(boot_null <= est) / n_bs
				)))
			} else if (type == "symmetric") {
				# Hall & Wilson (1991) symmetric test: pool both tails via absolute deviations.
				# p = P(|T* - mean(T*)| >= |t_obs - delta|)
				D_obs = abs(est - delta)
				D_boot = abs(boot_distr - mean(boot_distr))
				n_bs = length(D_boot)
				min(1, max(1 / n_bs, mean(D_boot >= D_obs)))
			} else if (type %in% c("studentized", "bootstrap-t")) {
				# Studentized (bootstrap-t) p-value: pivot by standard error.
				# t_obs = (est - delta) / se_hat
				# t*_b  = (T*_b  - est)  / se*_b    [centred at original estimate, not null]
				# p = P(|t*_b| >= |t_obs|)
				se_hat = tryCatch(private$infer_original_se(), error = function(e) NA_real_)
				if (!is.finite(se_hat) || se_hat <= 0) return(NA_real_)
				t_obs = abs(est - delta) / se_hat
				se_boot = boot_stats$se
				ok = is.finite(boot_distr) & is.finite(se_boot) & se_boot > 0
				t_boot = abs(boot_distr[ok] - est) / se_boot[ok]
				t_boot = t_boot[is.finite(t_boot)]
				if (length(t_boot) == 0L) return(NA_real_)
				min(1, max(1 / length(t_boot), mean(t_boot >= t_obs)))
			} else if (type == "bca") {
				# BCa p-value via closed-form CI inversion (Efron 1987; Efron & Tibshirani 1993).
				private$pval_bca(boot_distr, est, delta)
			}
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
		compute_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE){
			private$assert_design_supports_resampling("Bootstrap inference")
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			type = tolower(private$get_bootstrap_type(type))
			assertChoice(type, c(
				"percentile", "basic", "studentized", "bootstrap-t",
				"symmetric-percentile-t", "bca", "prepivoted",
				"double-bootstrap", "calibrated", "smoothed"
			))

			est = as.numeric(self$compute_treatment_estimate(estimate_only = FALSE))
			if (length(est) == 0L || !is.finite(est[1])) stop("Bootstrap confidence interval returned NA bounds")
			est = est[1]

			# BCa only needs theta (not se), so route it through the parallel bootstrap path.
			# studentized/symmetric-percentile-t need se per replicate → serial statistics path.
			boot_stats = if (type %in% c("studentized", "bootstrap-t", "symmetric-percentile-t", "prepivoted", "double-bootstrap", "calibrated", "smoothed")) {
				private$approximate_bootstrap_statistics_beta_hat_T(B = B, show_progress = show_progress, na.rm = na.rm, smooth = identical(type, "smoothed"))
			} else {
				list(theta = self$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = show_progress), se = NULL)
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

		assert_valid_bootstrap_type = function(bootstrap_type){
			if (is.null(bootstrap_type)) return(invisible(NULL))
			assertChoice(bootstrap_type, c("within_blocks", "resample_blocks"))
			valid_blocking_classes = c("FixedDesignBlocking", "FixedDesignOptimalBlocks", "DesignSeqOneByOneSPBR", "FixedDesignBlockedCluster")
			if (!any(vapply(valid_blocking_classes, function(cls) is(private$des_obj, cls), logical(1)))){
				stop("bootstrap_type can only be set for blocking designs: ", paste(valid_blocking_classes, collapse = ", "))
			}
			invisible(NULL)
		},

		get_bootstrap_type = function(type) {
			if (!is.null(type)) return(type)
			edi_bootstrap_dispatch_policy(class(self), object = self)
		},

		supports_reusable_bootstrap_worker = function(){
			FALSE
		},

		create_bootstrap_worker_state = function(){
			NULL
		},

		create_design_backed_bootstrap_worker_state = function(){
			worker = self$duplicate(verbose = FALSE, make_fork_cluster = FALSE)
			worker$num_cores = 1L
			worker_priv = worker$.__enclos_env__$private
			worker_des = if (!is.null(worker_priv$des_obj)) worker_priv$des_obj$duplicate(verbose = FALSE) else NULL
			worker_des_priv = if (!is.null(worker_des)) worker_des$.__enclos_env__$private else NULL
			source_des_priv = private$des_obj_priv_int
			if (!is.null(worker_des)) {
				worker_priv$des_obj = worker_des
				worker_priv$des_obj_priv_int = worker_des_priv
			}
			worker_priv$X = private$get_X()
			list(
				worker = worker,
				worker_priv = worker_priv,
				worker_des_priv = worker_des_priv,
				base_Xraw = if (!is.null(source_des_priv$Xraw)) source_des_priv$Xraw else NULL,
				base_Ximp = if (!is.null(source_des_priv$Ximp)) source_des_priv$Ximp else NULL,
				base_X = if (!is.null(private$X)) private$X else private$get_X(),
				base_w = if (!is.null(source_des_priv$w)) as.numeric(source_des_priv$w) else NULL,
				base_y = if (!is.null(source_des_priv$y)) source_des_priv$y else NULL,
				base_dead = if (!is.null(source_des_priv$dead)) as.numeric(source_des_priv$dead) else NULL,
				base_m = if (!is.null(source_des_priv$m)) source_des_priv$m else NULL,
				n = private$n
			)
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			stop("Reusable bootstrap workers are not implemented for this class.")
		},

		load_bootstrap_sample_into_design_backed_worker = function(worker_state, indices){
			indices = as.integer(indices)
			use_resampled_covariates =
				length(indices) != worker_state$n ||
				is(worker_state$worker_priv$des_obj, "FixedDesignBlockedCluster")
			w_priv = worker_state$worker_priv
			w_priv$X = if (use_resampled_covariates && !is.null(worker_state$base_X)) {
				worker_state$base_X[indices, , drop = FALSE]
			} else {
				worker_state$base_X
			}
			w_priv$w = if (!is.null(worker_state$base_w)) as.numeric(worker_state$base_w[indices]) else NULL
			w_priv$y = if (!is.null(worker_state$base_y)) as.numeric(worker_state$base_y[indices]) else NULL
			w_priv$dead = if (!is.null(worker_state$base_dead)) as.numeric(worker_state$base_dead[indices]) else NULL
			w_priv$y_temp = w_priv$y
			if (!is.null(worker_state$base_m)) w_priv$m = worker_state$base_m[indices]
			w_priv$n = if (use_resampled_covariates) length(indices) else worker_state$n
			w_priv$cached_values = list()
			w_priv$reduced_design_keep_cache = NULL
			w_priv$fixed_covariate_keep_cache = NULL
			w_priv$cached_mod = NULL

			des_priv = worker_state$worker_des_priv
			if (!is.null(des_priv)) {
				if (use_resampled_covariates) {
					des_priv$X = w_priv$X
					des_priv$t = length(indices)
					des_priv$n = length(indices)
				}
				des_priv$w = w_priv$w
				des_priv$y = w_priv$y
				des_priv$dead = w_priv$dead
				if (!is.null(worker_state$base_m)) des_priv$m = w_priv$m
				des_priv$all_subject_data_cache = list()
				if (use_resampled_covariates) des_priv$lin_centered_covariates = NULL
			}
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			stop("Reusable bootstrap workers are not implemented for this class.")
		},

		compute_bootstrap_worker_estimate_via_compute_treatment_estimate = function(worker_state){
			as.numeric(worker_state$worker$compute_treatment_estimate(estimate_only = TRUE))[1L]
		},

		compute_bootstrap_distribution_with_reused_workers = function(B, actual_cores, show_progress = FALSE, bootstrap_type = NULL){
			chunk_n = max(1L, min(as.integer(actual_cores), as.integer(B)))
			chunk_id = ceiling(seq_len(B) / ceiling(B / chunk_n))
			chunks = split(seq_len(B), chunk_id)

			run_chunk = function(idxs) {
				worker_state = private$create_bootstrap_worker_state()
				out = numeric(length(idxs))
				for (k in seq_along(idxs)) {
					boot_draw = private$bootstrap_sample_indices(private$n, bootstrap_type)
					out[k] = tryCatch({
						private$load_bootstrap_sample_into_worker(worker_state, boot_draw$i_b)
						private$compute_bootstrap_worker_estimate(worker_state)
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

		bootstrap_sample_indices = function(n, bootstrap_type = NULL){
			if (!is.null(private$des_obj)){
				return(private$des_obj_priv_int$draw_bootstrap_indices(bootstrap_type))
			}
			list(i_b = sample.int(n, n, replace = TRUE), m_vec_b = NULL)
		},

		bootstrap_subset_inference = function(boot_draw, smooth = FALSE){
			# boot_draw is list(i_b, m_vec_b) from des_obj$draw_bootstrap_indices()
			# For backward compatibility also accept a bare integer vector.
			if (is.list(boot_draw)){
				indices = as.integer(boot_draw$i_b)
				m_vec_b = boot_draw$m_vec_b
			} else {
				indices = as.integer(boot_draw)
				m_vec_b = NULL
			}
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
			# Use Design-provided m_vec_b (pair-aware) if available; otherwise subset original m
			if (!is.null(orig_des_priv$m)) {
				sub_des_priv$m = if (!is.null(m_vec_b)) as.integer(m_vec_b) else as.integer(orig_des_priv$m[indices])
			}
			if (!is.null(orig_des_priv$y_i_t_i)) sub_des_priv$y_i_t_i = orig_des_priv$y_i_t_i[indices]
			sub_des_priv$all_subject_data_cache = list()
			sub_des_priv$p_raw_t = if (!is.null(sub_des_priv$Xraw)) ncol(sub_des_priv$Xraw) else NULL
			sub_des_priv$t = length(indices)
			sub_des_priv$n = length(indices)
			sub_des_priv$fixed_sample = TRUE
			# Reset bootstrap structure cache — indices changed, pair structure no longer valid
			sub_des_priv$kk_boot_pair_rows   = NULL
			sub_des_priv$kk_boot_i_reservoir = NULL
			sub_des_priv$kk_boot_n_reservoir  = NULL
			# Reset model-specific design caches that hold n-row matrices — must be recomputed
			# on the subset (e.g. Lin centered covariates, which are cached as an n-row Xc).
			sub_des_priv$lin_centered_covariates = NULL
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
			sub_inf_priv$X = if (!is.null(private$X)) subset_field(private$X) else NULL
			sub_inf_priv$y = sub_des_priv$y
			sub_inf_priv$y_temp = sub_des_priv$y
			sub_inf_priv$w = sub_des_priv$w
			sub_inf_priv$dead = sub_des_priv$dead
			sub_inf_priv$n = length(indices)
			sub_inf_priv$cached_values = list()
			sub_inf_priv$cached_values$rand_distr_cache = list()
			sub_inf_priv$cached_values$m_cache = list()
			sub_inf_priv$reduced_design_keep_cache = NULL
			sub_inf_priv$fixed_covariate_keep_cache = NULL
			sub_inf_priv$cached_mod = NULL
			if (!is.null(sub_des_priv$m)) sub_inf_priv$m = sub_des_priv$m
			sub_inf
		},

		bootstrap_replication_stats = function(boot_draw, smooth = FALSE){
			sub_inf = private$bootstrap_subset_inference(boot_draw, smooth = smooth)
			if (is.null(sub_inf)) return(c(theta = NA_real_, se = NA_real_))
			tryCatch({
				# Use the return value of compute_treatment_estimate() for theta.
				# Reading cached_values$beta_hat_T is unreliable: some classes (e.g.
				# GComp, KMDiff) return estimates directly without storing beta_hat_T.
				theta = as.numeric(sub_inf$compute_treatment_estimate(estimate_only = FALSE))[1L]
				se = as.numeric(sub_inf$.__enclos_env__$private$cached_values$s_beta_hat_T)[1L]
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
				stats_mat[b, ] = private$bootstrap_replication_stats(idx, smooth = smooth)  # idx is a boot_draw list
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
			actual_cores = private$effective_parallel_cores("jackknife", self$num_cores)
			jack = unlist(private$par_lapply(seq_len(n), function(i) {
				idx = seq_len(n)[-i]
				sub_inf = private$bootstrap_subset_inference(idx, smooth = FALSE)
				if (is.null(sub_inf)) return(NA_real_)
				tryCatch({
					theta = as.numeric(sub_inf$compute_treatment_estimate(estimate_only = TRUE))[1L]
					if (is.finite(theta)) theta else NA_real_
				}, error = function(e) NA_real_)
			}, n_cores = actual_cores, show_progress = FALSE), use.names = FALSE)
			as.numeric(jack)
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
					idx = private$bootstrap_sample_indices(private$des_obj$get_n())  # returns boot_draw list
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
		},

		# BCa p-value via closed-form CI inversion (Efron 1987; Efron & Tibshirani 1993).
		# Derives the bias-correction (z0) and acceleration (a) constants from the bootstrap
		# distribution and jackknife estimates, then computes the BCa-adjusted CDF value at
		# delta_0.  The two-sided p-value is  2 * min(Phi(adj_z), 1 - Phi(adj_z))  where
		#   z_delta = Phi^{-1}(F_boot(delta_0)),   s = z_delta - z0,
		#   adj_z   = s / (1 + a*s) - z0.
		pval_bca = function(boot_distr, est, delta){
			jack = private$approximate_jackknife_distribution_beta_hat_T()
			jack = jack[is.finite(jack)]
			if (length(jack) < 2L) stop("BCa p-value requires jackknife estimates.")
			z0 = stats::qnorm(mean(boot_distr < est))
			jack_bar = mean(jack)
			num = sum((jack_bar - jack)^3)
			den = 6 * (sum((jack_bar - jack)^2)^(3/2))
			a = if (is.finite(den) && den > 0) num / den else 0
			p_delta = mean(boot_distr < delta)
			p_delta = pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, p_delta))
			z_delta = stats::qnorm(p_delta)
			s = z_delta - z0
			denom = 1 + a * s
			if (!is.finite(denom) || abs(denom) < .Machine$double.eps) return(NA_real_)
			adj_z = s / denom - z0
			if (!is.finite(adj_z)) return(NA_real_)
			min(1, 2 * min(stats::pnorm(adj_z), 1 - stats::pnorm(adj_z)))
		}
	)
)
