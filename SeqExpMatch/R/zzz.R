.onLoad = function(libname, pkgname) {
  if (is.null(getOption("datatable.quiet"))) {
    options(datatable.quiet = TRUE)
  }
}

.onAttach = function(libname, pkgname){
	packageStartupMessage(
			paste("Welcome to SeqExpMatch v", 
					utils::packageVersion("SeqExpMatch"), 
					sep = "")
	)	
}