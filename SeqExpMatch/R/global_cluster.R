# Internal environment for the EDI package to store global state
edi_env = new.env(parent = emptyenv())

#' Set the number of cores for parallelization
#' 
#' This function initializes a persistent parallel cluster (either a fork cluster
#' on Unix-like systems or a mirai cluster on others) to be used by all Design
#' and Inference objects. This avoids the overhead of creating clusters
#' repeatedly.
#' 
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

  if (as.integer(num_cores) <= 1L) {
    set_package_threads(1L)
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
  
  # Set package threads globally
  set_package_threads(num_cores)
  
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
    parallel::stopCluster(edi_env$global_fork_cluster)
    edi_env$global_fork_cluster = NULL
  }
  
  # Handle mirai cluster
  if (!is.null(edi_env$global_mirai_num_cores)) {
    if (requireNamespace("mirai", quietly = TRUE)) {
      mirai::daemons(0) # Stop daemons
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
    try(fixest::setFixest_nthreads(num_cores), silent = TRUE)
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
