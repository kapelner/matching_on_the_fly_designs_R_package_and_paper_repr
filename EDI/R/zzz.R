#' @importFrom pbmcapply pbmclapply
#' @importFrom randomizr block_ra block_and_cluster_ra cluster_ra
NULL

.onLoad = function(libname, pkgname) {
	if (is.null(getOption("datatable.quiet"))) {
		options(datatable.quiet = TRUE)
	}
	
	# Set Java heap size for rJava-based dependencies (like GreedyExperimentalDesign)
	# This must be set before the JVM is initialized.
	if (is.null(getOption("java.parameters"))) {
		options(java.parameters = "-Xmx10g")
	}
}

.onAttach = function(libname, pkgname){
	packageStartupMessage(
			paste("Welcome to EDI v",
					utils::packageVersion("EDI"),
					sep = "")
	)
}
