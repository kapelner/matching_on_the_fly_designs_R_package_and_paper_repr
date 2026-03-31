# Internal environment for the EDI package to store global state
edi_env = new.env(parent = emptyenv())

#' Create a Global Fork Cluster
#' 
#' This function creates a persistent socket-based fork cluster using the specified
#' number of cores and stores it within the package environment. This cluster can
#' be automatically reused by inference objects, avoiding the ~300ms startup penalty
#' incurred by creating a cluster for each individual object.
#' 
#' @param num_cores The number of CPU cores to use for the cluster.
#' 
#' @details This function is only available on Linux systems. On other operating
#' systems, socket-based fork clusters are either unsupported or significantly 
#' less efficient. For cross-platform parallelization, specify \code{num_cores} 
#' directly in the design or inference object constructors.
#' 
#' @return Invisible \code{NULL}.
#' 
#' @export
create_global_fork_cluster = function(num_cores) {
  if (Sys.info()["sysname"] != "Linux") {
    stop("Global fork clusters are only supported on Linux. For other systems, please specify num_cores in the design and inference initialization functions.")
  }
  
  checkmate::assertCount(num_cores, positive = TRUE)
  
  # Clear existing cluster if it exists
  if (!is.null(edi_env$global_fork_cluster)) {
    stop_global_fork_cluster()
  }
  
  edi_env$global_fork_cluster = parallel::makeForkCluster(num_cores)
  invisible(NULL)
}

#' Stop the Global Fork Cluster
#' 
#' This function stops the global fork cluster stored in the package environment
#' and clears the reference.
#' 
#' @return Invisible \code{NULL}.
#' 
#' @export
stop_global_fork_cluster = function() {
  if (is.null(edi_env$global_fork_cluster)) {
    warning("No global fork cluster exists to clear.")
    return(invisible(NULL))
  }
  
  parallel::stopCluster(edi_env$global_fork_cluster)
  edi_env$global_fork_cluster = NULL
  invisible(NULL)
}

# Internal helper to get the global cluster
get_global_fork_cluster = function() {
  edi_env$global_fork_cluster
}
