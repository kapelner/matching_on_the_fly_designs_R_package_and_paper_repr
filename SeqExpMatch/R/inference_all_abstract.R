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
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, verbose = FALSE){
			assertClass(des_obj, "Design")

			assertFlag(verbose)
			des_obj$assert_all_responses_recorded()

			private$cached_values = list()
			private$any_censoring = des_obj$any_censoring()
			private$des_obj = des_obj
			private$des_obj_priv_int = des_obj$.__enclos_env__$private
			private$y = private$des_obj_priv_int$y
			private$y_temp = private$y
			private$w = private$des_obj_priv_int$w
			private$dead = private$des_obj_priv_int$dead
			private$is_KK = is(des_obj, "DesignSeqOneByOneKK14")
			private$has_match_structure = private$is_KK || is(des_obj, "FixedDesignBinaryMatch")
			private$n = des_obj$get_n()
			private$prob_T = des_obj$get_prob_T()
			private$supports_design_resampling = isTRUE(des_obj$supports_resampling())
			
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
		#' @param make_fork_cluster 	Whether the duplicate should be allowed to create a fork cluster. Default FALSE.
		#' @return 			A new \code{Inference} object with the same data
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
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

	active = list(
		#' @field num_cores Current number of cores for this inference object.
		#'   Defaults to the global budget unless overridden on the object.
		num_cores = function(value) {
			if (missing(value)) {
				if (!is.null(private$num_cores_override)) return(private$num_cores_override)
				return(get_num_cores())
			}
			checkmate::assertCount(value, positive = TRUE)
			private$num_cores_override = as.integer(value)
			invisible(self)
		}
	),

	private = list(
		finalize = function(){
			# We no longer own the cluster, it is global.
		},
		des_obj = NULL,		des_obj_priv_int = NULL,
		m = NULL,
		is_KK = NULL,
		has_match_structure = NULL,
		supports_design_resampling = FALSE,
		any_censoring = NULL,
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
			num_cores_override = NULL,

		# Returns the number of C++ OpenMP threads to use for a parallel C++ function
		# with n_work_items items of work. Caps threads so that each thread handles
		# at least 10 items; for tiny r/B values this prevents thread-management
		# overhead from dominating over the actual computation.
		n_cpp_threads = function(n_work_items) {
			min(self$num_cores, max(1L, as.integer(n_work_items) %/% 10L))
		},

		parallel_dispatch_policy = function(operation) {
			edi_parallel_dispatch_policy(
				inference_class = class(self)[1],
				response_type = private$des_obj$get_response_type(),
				operation = operation
			)
		},

		effective_parallel_cores = function(operation, requested_cores = self$num_cores) {
			requested_cores = max(1L, as.integer(requested_cores))
			policy = private$parallel_dispatch_policy(operation)
			if (requested_cores > 1L && isTRUE(policy$force_serial)) {
				return(1L)
			}
			requested_cores
		},

			par_lapply = function(X, FUN, n_cores = self$num_cores, budget = 1L, show_progress = FALSE, export_list = NULL){
				if (length(X) == 0L) return(list())
				n_cores = max(1L, min(as.integer(n_cores), length(X)))
				budget = max(1L, as.integer(budget))
				chunk_count = min(length(X), max(1L, 4L * n_cores))
				chunk_size = max(1L, ceiling(length(X) / chunk_count))
				chunks = split(X, ceiling(seq_along(X) / chunk_size))

				# Run a whole chunk under the requested worker budget so we do not
				# pay scheduler/export overhead once per iteration.
				RUN_CHUNK = function(chunk) {
					ns = asNamespace("EDI")
					edi_env = ns$edi_env
					prev_override = edi_env$num_cores_override
					prev_threads = getOption(".edi_last_set_threads")
					if (is.null(prev_threads) || length(prev_threads) != 1L || !is.finite(prev_threads)) {
						prev_threads = 1L
					}
					edi_env$num_cores_override = budget
					ns$set_package_threads(budget)
					on.exit({
						edi_env$num_cores_override = prev_override
						ns$set_package_threads(prev_threads)
					}, add = TRUE)
					lapply(chunk, FUN)
				}

				flatten_chunk_results = function(results) {
					if (length(results) == 0L) return(list())
					unlist(results, recursive = FALSE, use.names = FALSE)
				}

				if (n_cores <= 1L) return(RUN_CHUNK(X))

				global_cl = get_global_fork_cluster()
				global_mirai_cores = get_global_mirai_cores()

				if (!is.null(global_cl)){
					worker_cl = global_cl[seq_len(min(n_cores, length(global_cl)))]
					tryCatch({
						if (!is.null(export_list) && length(export_list) > 0L) {
							export_env = list2env(export_list, parent = emptyenv())
							parallel::clusterExport(worker_cl, names(export_list), envir = export_env)
						}
						flatten_chunk_results(parallel::parLapply(worker_cl, chunks, RUN_CHUNK))
					}, error = function(e) {
					# If the persistent cluster has been killed externally (for example by a
					# timeout watchdog), clear it and fall back to a one-shot Unix fork apply.
					# This prevents stale socket connections from poisoning subsequent calls.
					msg = conditionMessage(e)
						if (.Platform$OS.type == "unix" &&
							grepl("connection|serialize|unserialize|postNode|sendData|recvData", msg, ignore.case = TRUE)) {
							edi_env$global_fork_cluster = NULL
							try(parallel::stopCluster(global_cl), silent = TRUE)
							if (isTRUE(show_progress) && requireNamespace("pbmcapply", quietly = TRUE)) {
								return(flatten_chunk_results(pbmcapply::pbmclapply(chunks, RUN_CHUNK, mc.cores = n_cores)))
							}
							return(flatten_chunk_results(parallel::mclapply(chunks, RUN_CHUNK, mc.cores = n_cores)))
						}
							stop(e)
					})
				} else if (!is.null(global_mirai_cores)){
					requested_mirai_cores = min(n_cores, global_mirai_cores)
					private$ensure_mirai_daemons(requested_mirai_cores)
					on.exit(private$ensure_mirai_daemons(global_mirai_cores), add = TRUE)
					tasks = lapply(chunks, function(chunk) mirai::mirai({RUN_CHUNK(chunk)}, RUN_CHUNK = RUN_CHUNK, chunk = chunk))
					flatten_chunk_results(lapply(tasks, function(m) m[]))
				} else if (.Platform$OS.type != "unix"){
					if (!isTRUE(private$warned_no_parallel)){
						message("Parallelism (num_cores > 1) requires the 'mirai' package on non-Unix systems. Install it with install.packages('mirai'). Falling back to serial computation.")
						private$warned_no_parallel = TRUE
					}
					RUN_CHUNK(X)
				} else {
					if (isTRUE(show_progress) && requireNamespace("pbmcapply", quietly = TRUE)){
						flatten_chunk_results(pbmcapply::pbmclapply(chunks, RUN_CHUNK, mc.cores = n_cores))
					} else {
						flatten_chunk_results(parallel::mclapply(chunks, RUN_CHUNK, mc.cores = n_cores))
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

			get_or_create_fork_cluster = function(){
				cl = get_global_fork_cluster()
				if (is.null(cl)) {
					stop("No global fork cluster is initialized. Call set_num_cores() first.")
				}
				cl
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
