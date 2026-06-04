options(keep.source = TRUE)

load_source_with_filenames = function(path) {
	env = new.env(parent = globalenv())
	methods::setPackageName("roxygen_devtest", env)

	deps = desc::desc_get_deps(path)
	pkgs = deps$package[deps$type %in% c("Depends", "Imports") & deps$package != "R"]
	lapply(pkgs, require, character.only = TRUE)

	lapply(roxygen2:::package_files(path), sys.source, envir = env, keep.source = TRUE)
	env
}

original_extract_r6_methods = roxygen2:::extract_r6_methods
patched_extract_r6_methods = function(x) {
	methods = original_extract_r6_methods(x)
	if (!("initialize" %in% methods$name)) {
		class_file = methods$file[1]
	} else {
		class_file = methods$file[match("initialize", methods$name)]
	}
	if (is.na(class_file)) {
		return(methods)
	}
	methods[is.na(methods$file) | methods$file == class_file, , drop = FALSE]
}
unlockBinding("extract_r6_methods", asNamespace("roxygen2"))
assign("extract_r6_methods", patched_extract_r6_methods, envir = asNamespace("roxygen2"))
lockBinding("extract_r6_methods", asNamespace("roxygen2"))

roxygen2::roxygenize("EDI", load_code = load_source_with_filenames)

