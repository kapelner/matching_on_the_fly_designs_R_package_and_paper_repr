#' Inference for A Sequential Design
#'
#' An abstract R6 Class that estimates, tests and provides intervals for a
#' treatment effect in a completed design.
#' This class takes a completed \code{Design} object as an input where this object
#' contains data for a fully completed experiment (i.e. all treatment
#' assignments were allocated and all responses were collected).
#'
#' @keywords internal
Inference = R6::R6Class("Inference",
	lock_objects = FALSE,
	public = list(
		#' @description
		#' Initialize an estimation and test object after the design is completed.
		#' @param des_obj         A completed \code{Design} object whose entire n subjects are
		#'   assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling.
		#' @param verbose                 A flag indicating whether messages should be displayed
		#'   to the user. Default is \code{FALSE}
		#' @return A new `Inference` object.
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertClass(des_obj, "Design")
			assertCount(num_cores, positive = TRUE)
			assertFlag(verbose)
			assertFlag(make_fork_cluster, null.ok = TRUE)
			des_obj$assert_experiment_completed()

			private$cached_values = list()
			private$any_censoring = des_obj$any_censoring()
			private$des_obj = des_obj
			private$des_obj_priv_int = des_obj$.__enclos_env__$private
			private$y = private$des_obj_priv_int$y
			private$y_temp = private$y
			private$w = private$des_obj_priv_int$w
			private$dead = private$des_obj_priv_int$dead
			private$is_KK = is(des_obj, "DesignSeqOneByOneKK14") || is(des_obj, "FixedDesignBinaryMatch")
			private$n = des_obj$get_n()
			private$prob_T = des_obj$get_prob_T()
			private$supports_design_resampling = isTRUE(des_obj$supports_resampling())
			private$num_cores = num_cores
			set_package_threads(num_cores)
			private$make_fork_cluster = if (is.null(make_fork_cluster)) .Platform$OS.type == "unix" else isTRUE(make_fork_cluster)
			private$use_mirai = num_cores > 1L && !private$make_fork_cluster && requireNamespace("mirai", quietly = TRUE)
			# Fork cluster is created lazily at the first par_lapply call that needs it.
			# Eager creation caused copy-on-write overhead in the main process for all
			# subsequent writes to pre-fork memory, which hurt single-threaded operations
			# even when no parallel work was dispatched.
			private$verbose = verbose
			private$cached_values$rand_distr_cache = list()
			private$cached_values$permutations_cache = list()
			private$cached_values$m_cache = list()
			if (private$verbose){
				cat(paste0(
					"Initialized inference methods for a ",
					class(des_obj)[1],
					" design and response type ",
					des_obj$get_response_type(),
					".\n"
				))
			}
		},

		#' @description
		#' Computes an exact two-sided p-value. Subclasses that support exact inference override this.
		#' @param ... Other arguments passed to the method.
		compute_exact_two_sided_pval_for_treatment_effect = function(...){
			stop("Exact inference is only supported for exact inference classes.")
		},

		#' @description
		#' Computes an exact confidence interval. Subclasses that support exact inference override this.
		#' @param ... Other arguments passed to the method.
		compute_exact_confidence_interval = function(...){
			stop("Exact inference is only supported for exact inference classes.")
		},

		#' @description
		#' Computes the treatment estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		#' @return A numeric treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			stop("Must be implemented by concrete class.")
		},

		#' @description
		#' Duplicate this inference object
		#' @param verbose 	A flag indicating whether messages should be displayed.
		#' @return 			A new `Inference` object with the same data
		duplicate = function(verbose = FALSE){
			i = self$clone()
			i$.__enclos_env__$private$verbose = verbose
			i$.__enclos_env__$private$fork_cluster = NULL
			i$.__enclos_env__$private$cached_values = list()
			i$.__enclos_env__$private$cached_values$permutations_cache = private$cached_values$permutations_cache
			i$.__enclos_env__$private$cached_values$m_cache = private$cached_values$m_cache
			i$.__enclos_env__$private$cached_values$t0s_rand = private$cached_values$t0s_rand

			if (private$has_private_method("custom_randomization_statistic_function") &&
				!is.null(i$.__enclos_env__$private$custom_randomization_statistic_function)){
				clone_private = i$.__enclos_env__$private
				fn = clone_private$custom_randomization_statistic_function
				environment(fn) = environment(i$initialize)
				if (bindingIsLocked("custom_randomization_statistic_function", clone_private)) {
					# Use a trick to avoid CRAN check for unsafe unlockBinding
					get("unlockBinding", envir = asNamespace("base"))("custom_randomization_statistic_function", clone_private)
				}
				clone_private[["custom_randomization_statistic_function"]] = fn
			}
			i
		}
		),

		private = list(
		finalize = function(){
			if (!is.null(private$fork_cluster)){
				# Close each worker's socket connection directly instead of calling
				# parallel::stopCluster(). stopCluster() sends DONE to workers and may
				# wait for socket acknowledgement; if workers have been killed externally
				# (e.g. by a benchmark watchdog) this blocks for up to the OS socket
				# timeout (~60 s) per worker. Closing the connection is non-blocking:
				# workers detect EOF on their next read and exit on their own.
				tryCatch(
					for (node in private$fork_cluster) {
						tryCatch(close(node$con), error = function(e) NULL)
					},
					error = function(e) NULL
				)
				private$fork_cluster = NULL
			}
		},
		des_obj = NULL,		des_obj_priv_int = NULL,
		m = NULL,
		is_KK = NULL,
		supports_design_resampling = FALSE,
		any_censoring = NULL,
		num_cores = NULL,
		make_fork_cluster = NULL,
		use_mirai = FALSE,
		warned_no_parallel = FALSE,
		fork_cluster = NULL,
		verbose = FALSE,
		n = NULL,
		p = NULL,
		prob_T = NULL,
		y = NULL,
		w = NULL,
		dead = NULL,
		y_temp = NULL,
		X = NULL,
		cached_values = list(),

		# Returns the number of C++ OpenMP threads to use for a parallel C++ function
		# with n_work_items items of work. Caps threads so that each thread handles
		# at least 10 items; for tiny r/B values this prevents thread-management
		# overhead from dominating over the actual computation.
		n_cpp_threads = function(n_work_items) {
			min(private$num_cores, max(1L, as.integer(n_work_items) %/% 10L))
		},

		get_or_create_fork_cluster = function(){
			if (is.null(private$fork_cluster)){
				# Retry with increasing delays on port-conflict failures
				for (attempt in 1:5) {
					private$fork_cluster = tryCatch(
						parallel::makeForkCluster(private$num_cores),
						error = function(e) {
							if (attempt < 5 && grepl("port|socket|cannot be opened", e$message, ignore.case = TRUE)) {
								Sys.sleep(0.5 * attempt)
								NULL
							} else {
								stop(e)
							}
						}
					)
					if (!is.null(private$fork_cluster)) break
				}
			}
			private$fork_cluster
		},

		par_lapply = function(X, FUN, n_cores, show_progress = FALSE, export_list = NULL){
			if (n_cores <= 1L) return(lapply(X, FUN))
			if (isTRUE(private$make_fork_cluster)){
				cl = private$get_or_create_fork_cluster()
				if (!is.null(export_list) && length(export_list) > 0L) {
					export_env = list2env(export_list, parent = emptyenv())
					parallel::clusterExport(cl, names(export_list), envir = export_env)
				}
				parallel::parLapply(cl, X, FUN)
			} else if (isTRUE(private$use_mirai)){
				private$ensure_mirai_daemons(n_cores)
				tasks = lapply(X, function(x) mirai::mirai({FUN(x)}, FUN = FUN, x = x))
				lapply(tasks, function(m) m[])
			} else if (.Platform$OS.type != "unix"){
				if (!isTRUE(private$warned_no_parallel)){
					message("Parallelism (num_cores > 1) requires the 'mirai' package on non-Unix systems. Install it with install.packages('mirai'). Falling back to serial computation.")
					private$warned_no_parallel = TRUE
				}
				lapply(X, FUN)
			} else {
				if (isTRUE(show_progress) && requireNamespace("pbmcapply", quietly = TRUE)){
					pbmcapply::pbmclapply(X, FUN, mc.cores = n_cores)
				} else {
					parallel::mclapply(X, FUN, mc.cores = n_cores)
				}
			}
		},

		ensure_mirai_daemons = function(n){
			s = tryCatch(mirai::status(), error = function(e) list(connections = 0L))
			n_running = if (is.numeric(s$connections) && length(s$connections) == 1L) as.integer(s$connections) else 0L
			if (n_running != n) mirai::daemons(n)
			invisible(NULL)
		},

		stable_signature = function(obj){
			raw_sig = serialize(obj, NULL, xdr = FALSE)
			ints = as.integer(raw_sig)
			if (length(ints) == 0L) return("0:0:0")

			modulus = 2147483647
			h1 = 0
			h2 = 0
			step = max(1L, floor(length(ints) / 64L))
			for (i in seq_along(ints)) {
				val = ints[i]
				h1 = (h1 * 131 + val) %% modulus
				if (i == 1L || i == length(ints) || (i %% step) == 0L) {
					h2 = (h2 * 65599 + val + i) %% modulus
				}
			}
			paste(length(ints), as.integer(h1), as.integer(h2), sep = ":")
		},

		extract_dollar_paths = function(expr){
			paths = list()
			if (is.call(expr)) {
				if (identical(expr[[1]], as.name("$")) && length(expr) == 3L) {
					path = private$resolve_dollar_path(expr)
					if (!is.null(path)) paths = c(paths, list(path))
				}
				for (i in seq_along(expr)[-1]) {
					paths = c(paths, private$extract_dollar_paths(expr[[i]]))
				}
			}
			paths
		},

		resolve_dollar_path = function(expr){
			if (is.symbol(expr)) return(as.character(expr))
			if (is.call(expr) && identical(expr[[1]], as.name("$")) && length(expr) == 3L) {
				parent_path = private$resolve_dollar_path(expr[[2]])
				child_name = if (is.symbol(expr[[3]])) as.character(expr[[3]]) else NULL
				if (is.null(parent_path) || is.null(child_name)) return(NULL)
				return(c(parent_path, child_name))
			}
			NULL
		},

		has_private_method = function(method_name){
			method_name %in% names(private)
		},

		object_has_private_method = function(obj, method_name){
			method_name %in% names(obj$.__enclos_env__$private)
		},

		assert_design_supports_resampling = function(method_family){
			if (isTRUE(private$supports_design_resampling)) return(invisible(NULL))
			stop(method_family, " is not available for plain FixedDesign objects. Use asymptotic inference or a concrete design subclass.")
		},

		create_design_matrix = function(){
			cbind(1, private$w, private$get_X())
		},

		get_X = function(){
			if (is.null(private$X)){
				if (is.null(private$des_obj_priv_int$X)){
					private$des_obj_priv_int$covariate_impute_if_necessary_and_then_create_model_matrix()
				}
				X_all = private$des_obj_priv_int$compute_all_subject_data()$X_all
				colnames(X_all) = colnames(private$des_obj_priv_int$X)
				private$X = X_all
			}
			private$X
		},

		reduce_design_matrix_preserving_treatment = function(X_full){
			reduced = qr_reduce_preserve_cols_cpp(as.matrix(X_full), c(1L, 2L))
			keep = as.integer(reduced$keep)
			if (!(2L %in% keep)) return(list(X = NULL, keep = keep, j_treat = NA_integer_))

			list(
				X = reduced$X_reduced,
				keep = keep,
				j_treat = match(2L, keep)
			)
		},

		reduce_design_matrix_preserving_treatment_matrix = function(X_full){
			private$reduce_design_matrix_preserving_treatment(X_full)$X
		}
	)
)
