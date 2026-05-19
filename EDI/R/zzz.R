#' @importFrom randomizr block_ra block_and_cluster_ra cluster_ra
NULL

.onLoad = function(libname, pkgname) {
	if (is.null(getOption("datatable.quiet"))) {
		options(datatable.quiet = TRUE)
	}
	
	# Set default for assertion execution
	if (is.null(getOption("edi.run_asserts"))) {
		options(edi.run_asserts = TRUE)
	}
	
}

.onAttach = function(libname, pkgname){
	version <- tryCatch(
		as.character(utils::packageDescription(pkgname, lib.loc = libname, fields = "Version")),
		error = function(e) NA_character_
	)
	if (!is.character(version) || length(version) != 1L || is.na(version) || version == "") {
		version <- "unknown"
	}
	packageStartupMessage(
			paste("Welcome to EDI v", version, sep = "")
	)
}
