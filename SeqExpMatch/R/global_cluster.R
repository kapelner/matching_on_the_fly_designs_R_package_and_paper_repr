# Internal environment for the EDI package to store global state
edi_env = new.env(parent = emptyenv())

#' Set the number of cores for parallelization
#' 
#' This function initializes a persistent parallel cluster (either a fork cluster
#' on Unix-like systems or a mirai cluster on others) to be used by all Design
#' and Inference objects. This avoids the overhead of creating clusters
#' repeatedly.
#'
#' @details
#' \code{set_num_cores()} sets a global upper bound for parallel work. It does
#' not guarantee that every inference routine will use all requested workers.
#' EDI's inference dispatcher applies a blocklist-first heuristic informed by
#' package benchmarks: workloads that have shown consistent multicore slowdowns
#' are forced to run serially, while the remaining workloads are allowed to use
#' their method-specific warmup heuristics and native thread caps.
#'
#' The current forced-serial blocklist covers incidence randomization confidence
#' intervals, bootstrap for non-regression KK Wilcoxon inference, bootstrap for
#' non-KK survival procedures, and bootstrap for incidence procedures. Do not
#' expect a universal "more cores is faster" rule.
#' 
#' @param num_cores Integer number of worker processes to make available.
#' @param force_mirai If \code{TRUE}, forces the use of the \code{mirai} package
#'   even on systems where forking is available.
#' 
#' @return Invisible \code{NULL}.
#' 
#' @export
set_num_cores = function(num_cores, force_mirai = FALSE) {
  checkmate::assertCount(num_cores, positive = TRUE)
  checkmate::assertFlag(force_mirai)
  
  # Clear any existing clusters first
  unset_num_cores()

  # Set package/native thread budgets before workers are created so forked
  # children inherit the intended global thread settings.
  set_package_threads(num_cores)

  if (as.integer(num_cores) <= 1L) {
    return(invisible(NULL))
  }
  
  if (force_mirai || .Platform$OS.type != "unix") {
    if (!requireNamespace("mirai", quietly = TRUE)) {
      stop("The 'mirai' package is required for parallelization on this system or when force_mirai = TRUE. Please install it.")
    }
    # Initialize mirai daemons
    mirai::daemons(num_cores)
    edi_env$global_mirai_num_cores = num_cores
  } else {
    # Unix-like system, use forking
    edi_env$global_fork_cluster = parallel::makeForkCluster(num_cores)
  }
  
  invisible(NULL)
}

#' Unset the number of cores and stop parallel clusters
#' 
#' This function stops any global fork or mirai clusters stored in the package
#' environment and resets the core count to serial execution.
#' 
#' @return Invisible \code{NULL}.
#' 
#' @export
unset_num_cores = function() {
  # Handle fork cluster
  if (!is.null(edi_env$global_fork_cluster)) {
    cl = edi_env$global_fork_cluster
    edi_env$global_fork_cluster = NULL
    tryCatch(parallel::stopCluster(cl), error = function(e) invisible(NULL))
  }
  
  # Handle mirai cluster
  if (!is.null(edi_env$global_mirai_num_cores)) {
    if (requireNamespace("mirai", quietly = TRUE)) {
      tryCatch(mirai::daemons(0), error = function(e) invisible(NULL)) # Stop daemons
    }
    edi_env$global_mirai_num_cores = NULL
  }
  
  # Reset package threads to 1
  set_package_threads(1L)
  
  invisible(NULL)
}

# Internal helper to get the global fork cluster
get_global_fork_cluster = function() {
  edi_env$global_fork_cluster
}

# Internal helper to get the global mirai core count
get_global_mirai_cores = function() {
  edi_env$global_mirai_num_cores
}

# Internal helper to get the current core count budget
get_num_cores = function() {
  if (!is.null(edi_env$num_cores_override)) return(edi_env$num_cores_override)
  cl = get_global_fork_cluster()
  if (!is.null(cl)) return(length(cl))
  mirai_cores = get_global_mirai_cores()
  if (!is.null(mirai_cores)) return(mirai_cores)
  1L
}

# Internal helper for empirical parallel dispatch policy.
# This is intentionally conservative and only forces serial execution for
# inference families that have shown repeat multicore regressions in the package
# benchmark suite. Everything else remains eligible for the usual warmup-based
# parallel heuristics. Do not expect a universal "more cores is faster" rule.
edi_parallel_dispatch_policy = function(inference_class, response_type, operation) {
  inference_class = as.character(inference_class[[1]])
  response_type = as.character(response_type[[1]])
  operation = as.character(operation[[1]])

  is_incid = grepl("^InferenceIncid", inference_class) ||
    grepl("^InferenceIncidence", inference_class)
  is_survival_nonkk = grepl("^InferenceSurvival", inference_class) &&
    !grepl("KK", inference_class, fixed = TRUE)
  is_nonregr_kk_wilcox = identical(inference_class, "InferenceAllKKWilcoxIVWC")

  reason = NULL

  if (identical(operation, "bootstrap")) {
    if (is_incid && identical(response_type, "incidence")) {
      reason = "bootstrap for incidence inference is forced serial by benchmark policy"
    } else if (is_survival_nonkk && identical(response_type, "survival")) {
      reason = "bootstrap for non-KK survival inference is forced serial by benchmark policy"
    } else if (is_nonregr_kk_wilcox) {
      reason = "bootstrap for non-regression KK Wilcoxon inference is forced serial by benchmark policy"
    }
  } else if (identical(operation, "rand_ci")) {
    if (is_incid && identical(response_type, "incidence")) {
      reason = "randomization confidence intervals for incidence inference are forced serial by benchmark policy"
    }
  }

  list(
    force_serial = !is.null(reason),
    reason = reason,
    inference_class = inference_class,
    response_type = response_type,
    operation = operation
  )
}

# Internal helper
set_package_threads = function(num_cores) {
  # Ensure it's an integer
  num_cores = as.integer(num_cores)

  # Check if we are already at the target to avoid slow Sys.setenv and other calls
  if (!is.null(options(".edi_last_set_threads")[[1]]) && options(".edi_last_set_threads")[[1]] == num_cores) {
    return(invisible(NULL))
  }

  # R packages with global thread setters
	  if (requireNamespace("data.table", quietly = TRUE)) {
	    data.table::setDTthreads(num_cores)
	  }
	  if (requireNamespace("fixest", quietly = TRUE)) {
	    suppressWarnings(try(fixest::setFixest_nthreads(num_cores), silent = TRUE))
	  }

  # Environment variables for OpenMP and BLAS/LAPACK
  # This helps prevent thread explosion in child processes
  # that call multi-threaded native libraries.
  # Sys.setenv is relatively slow, so we only do it if needed.
  Sys.setenv(OMP_NUM_THREADS =    num_cores)
  Sys.setenv(MKL_NUM_THREADS =    num_cores)
  Sys.setenv(OPENBLAS_NUM_THREADS =   num_cores)
  Sys.setenv(VECLIB_MAXIMUM_THREADS = num_cores)
  Sys.setenv(NUMEXPR_NUM_THREADS =    num_cores)
  options(".edi_last_set_threads" =   num_cores)
  invisible(NULL)
}
