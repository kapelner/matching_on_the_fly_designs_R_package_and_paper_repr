#' @importFrom pbmcapply pbmclapply
#' @importFrom randomizr block_ra block_and_cluster_ra cluster_ra
NULL

.onLoad = function(libname, pkgname) {
	if (is.null(getOption("datatable.quiet"))) {
	options(datatable.quiet = TRUE)
	}
}

.onAttach = function(libname, pkgname){
	packageStartupMessage(
			paste("Welcome to EDI v",
					utils::packageVersion("EDI"),
					sep = "")
	)
}
