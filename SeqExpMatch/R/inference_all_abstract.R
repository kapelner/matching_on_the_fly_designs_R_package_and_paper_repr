#' Inference for A Sequential Design
#'
#' @description
#' An abstract R6 Class that estimates, tests and provides intervals for a treatment effect in a completed design.
#' This class takes a completed \code{Design} object as an input where this object
#' contains data for a fully completed experiment (i.e. all treatment
#' assignments were allocated and all responses were collected).
#'
#' @keywords internal
Inference = R6::R6Class("Inference",
	public = list(
		# @description
		# Initialize an estimation and test object after the design is completed.
		# @param des_obj		A completed \code{Design} object whose entire n subjects are assigned and response y is recorded within.
		# @param num_cores			The number of CPU cores to use to parallelize the sampling.
		# @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{FALSE}
		# @return A new `Inference` object.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertClass(des_obj, "Design")
			assertCount(num_cores, positive = TRUE)
			assertFlag(verbose)
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

		# @description
		# Duplicate this inference object
		# @param verbose 	A flag indicating whether messages should be displayed.
		# @return 			A new `Inference` object with the same data
		duplicate = function(verbose = FALSE){
			i = self$clone()
			i$.__enclos_env__$private$verbose = verbose
			i$.__enclos_env__$private$cached_values = list()
			i$.__enclos_env__$private$cached_values$permutations_cache = private$cached_values$permutations_cache
			i$.__enclos_env__$private$cached_values$m_cache = private$cached_values$m_cache
			i$.__enclos_env__$private$cached_values$t0s_rand = private$cached_values$t0s_rand
			
			if (private$has_private_method("custom_randomization_statistic_function") && 
				!is.null(i$.__enclos_env__$private$custom_randomization_statistic_function)){
				clone_private = i$.__enclos_env__$private
				fn = clone_private$custom_randomization_statistic_function
				environment(fn) = environment(i$initialize)
				field = "custom_randomization_statistic_function"
				if (bindingIsLocked(field, clone_private)) {
					unlockBinding(field, clone_private)
				}
				assign(field, fn, envir = clone_private)
			}
			i
		}
	),

	private = list(
		des_obj = NULL,
		des_obj_priv_int = NULL,
		m = NULL,
		is_KK = NULL,
		supports_design_resampling = FALSE,
		any_censoring = NULL,
		num_cores = NULL,
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
