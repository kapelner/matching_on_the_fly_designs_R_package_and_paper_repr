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

# `methods$file[1]` (the old heuristic below) picks whichever method happens to
# be first in the merged public/private list. When a class composes in mixin
# methods via `c(SomeMixin$public, list(...))`, the mixin's methods (sourced
# from a DIFFERENT file) sort first, so `methods$file[1]` names the mixin's
# file instead of the class's own file -- every method actually documented in
# the class's own file then gets filtered OUT (roxygen2 reports their doc
# blocks as "Cannot find matching R6 method"), while the undocumented mixin
# methods get kept and flagged as "Undocumented R6 method". `block$file` -- the
# file of the class-level roxygen block right above `R6::R6Class(...)` -- is
# the reliable source of truth for which file a class "belongs to", so stash
# it from `object_defaults.r6class` (which has the block in scope) and prefer
# it in `extract_r6_methods` (which only receives the class generator).
.fast_roxygenize_current_block_file = new.env(parent = emptyenv())
.fast_roxygenize_package_name = desc::desc_get("Package", file = file.path("EDI", "DESCRIPTION"))[[1]]

original_object_defaults_r6class = roxygen2:::object_defaults.r6class
patched_object_defaults_r6class = function(x, block) {
	assign("file", block$file, envir = .fast_roxygenize_current_block_file)
	original_object_defaults_r6class(x, block)
}
# `object_defaults` dispatches via UseMethod(), which resolves `object_defaults.r6class`
# through roxygen2's internal `.__S3MethodsTable__.`, NOT via a plain namespace lookup --
# so the namespace binding below must be patched together with the S3 methods table entry,
# or dispatch keeps calling the original (unpatched) method.
unlockBinding("object_defaults.r6class", asNamespace("roxygen2"))
assign("object_defaults.r6class", patched_object_defaults_r6class, envir = asNamespace("roxygen2"))
lockBinding("object_defaults.r6class", asNamespace("roxygen2"))

s3_methods_table = get(".__S3MethodsTable__.", envir = asNamespace("roxygen2"))
unlockBinding("object_defaults.r6class", s3_methods_table)
assign("object_defaults.r6class", patched_object_defaults_r6class, envir = s3_methods_table)
lockBinding("object_defaults.r6class", s3_methods_table)

original_extract_r6_methods = roxygen2:::extract_r6_methods
patched_extract_r6_methods = function(x) {
	methods = original_extract_r6_methods(x)
	block_file = if (exists("file", envir = .fast_roxygenize_current_block_file, inherits = FALSE)) {
		get("file", envir = .fast_roxygenize_current_block_file, inherits = FALSE)
	} else {
		NA_character_
	}
	# `block$file` is an absolute path, but `methods$file` (from
	# `extract_r6_methods`) only stores the basename -- match on basename,
	# same as roxygen2's own `find_method_for_tag()` does.
	class_file = if (!is.na(block_file) && any(methods$file == basename(block_file), na.rm = TRUE)) {
		basename(block_file)
	} else if (!("initialize" %in% methods$name)) {
		methods$file[1]
	} else {
		methods$file[match("initialize", methods$name)]
	}
	if (is.na(class_file)) {
		return(methods)
	}
	methods[is.na(methods$file) | methods$file == class_file, , drop = FALSE]
}
unlockBinding("extract_r6_methods", asNamespace("roxygen2"))
assign("extract_r6_methods", patched_extract_r6_methods, envir = asNamespace("roxygen2"))
lockBinding("extract_r6_methods", asNamespace("roxygen2"))

# R6 stores source-loaded local superclasses with `package = ""` in this fast
# path. roxygen2 passes that value to rdtools::topic_exists(), whose current
# implementation eventually calls exists("") and aborts. Treat the empty package
# marker as the package currently being documented, matching the normal
# roxygenize/pkgload path for local R6 inheritance.
original_r6_extract_superclasses = roxygen2:::r6_extract_superclasses
patched_r6_extract_superclasses = function(r6data, env, class) {
	super = r6data$super
	cls = unique(super$classes$classname)
	if (length(cls) == 0) {
		return(roxygen2:::rd_r6_super(class))
	}
	pkgs = super$classes$package[match(cls, super$classes$classname)]
	current_package = roxygen2:::roxy_meta_get("current_package")
	if (is.null(current_package) || is.na(current_package) || !nzchar(current_package)) {
		current_package = .fast_roxygenize_package_name
	}
	pkgs[is.na(pkgs) | !nzchar(pkgs)] = current_package
	ht = purrr::map2_lgl(cls, pkgs, roxygen2:::has_topic)
	roxygen2:::rd_r6_super(class, package = pkgs, classname = cls, has_topic = ht)
}
unlockBinding("r6_extract_superclasses", asNamespace("roxygen2"))
assign("r6_extract_superclasses", patched_r6_extract_superclasses, envir = asNamespace("roxygen2"))
lockBinding("r6_extract_superclasses", asNamespace("roxygen2"))

original_r6_extract_inherited_methods = roxygen2:::r6_extract_inherited_methods
patched_r6_extract_inherited_methods = function(r6data) {
	inherited = original_r6_extract_inherited_methods(r6data)
	if (length(inherited$package) == 0) {
		return(inherited)
	}
	current_package = roxygen2:::roxy_meta_get("current_package")
	if (is.null(current_package) || is.na(current_package) || !nzchar(current_package)) {
		current_package = .fast_roxygenize_package_name
	}
	inherited$package[is.na(inherited$package) | !nzchar(inherited$package)] = current_package
	inherited
}
unlockBinding("r6_extract_inherited_methods", asNamespace("roxygen2"))
assign("r6_extract_inherited_methods", patched_r6_extract_inherited_methods, envir = asNamespace("roxygen2"))
lockBinding("r6_extract_inherited_methods", asNamespace("roxygen2"))

# format.rd_r6_inherited() is no longer monkey-patched here: the empty/NA
# package fallback this used to provide (for superclasses loaded via
# load_source_with_filenames() below, which don't go through pkgload) is now
# handled directly in roxygen2's own format.rd_r6_inherited(), which also
# respects the deep_public_method_links option -- patching it here would
# silently freeze out that flag with a stale, hardcoded rendering.

Rcpp::compileAttributes("EDI")
roxygen2::roxygenize("EDI", load_code = load_source_with_filenames)
