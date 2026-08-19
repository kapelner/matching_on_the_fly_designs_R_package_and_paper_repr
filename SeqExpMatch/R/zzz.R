.onAttach = function(libname, pkgname){
	packageStartupMessage(
			paste("Welcome to SeqExpMatch v",
					utils::packageVersion("SeqExpMatch"),
					sep = "")
	)
	packageStartupMessage(
			"SeqExpMatch is DEPRECATED and no longer maintained. All designs and inference ",
			"procedures in this package have been superseded by the 'EDI' package ",
			"(https://github.com/kapelner/EDI). Please migrate to EDI -- see ?SeqDesign and ",
			"?SeqDesignInference for specific migration code."
	)
}