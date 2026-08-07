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
	keep_files = class_file
	methods[is.na(methods$file) | methods$file %in% keep_files, , drop = FALSE]
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

# Workaround: this Rcpp version sometimes emits RcppExports.cpp without
# #include <RcppEigen.h>, even though many exported functions return/take
# Eigen types -- breaking the C++ compile ("'Eigen' does not name a type").
# Re-add it if compileAttributes() dropped it.
rcpp_exports_cpp = file.path("EDI", "src", "RcppExports.cpp")
rcpp_exports_lines = readLines(rcpp_exports_cpp)
if (!any(grepl("^#include <RcppEigen\\.h>$", rcpp_exports_lines))) {
	rcpp_include_idx = which(grepl("^#include <Rcpp\\.h>$", rcpp_exports_lines))[1]
	rcpp_exports_lines = append(rcpp_exports_lines, "#include <RcppEigen.h>", after = rcpp_include_idx - 1L)
	writeLines(rcpp_exports_lines, rcpp_exports_cpp)
}

roxygen2::roxygenize("EDI", load_code = load_source_with_filenames)

# KLUDGE, not a real fix: some "mixin" files document a plain public = list(...)
# object (spliced into multiple R6 classes elsewhere via c(SomeMixin$public,
# list(...))), not an R6::R6Class(...) call. roxygen2's R6 method-doc
# extraction (patched above for real R6 classes) doesn't apply to plain
# lists, so every per-method @param from every method in the list gets
# accumulated onto the ONE topic documenting the bare list object -- whose
# \usage{} is just that object's name, with no formals to match against.
# R CMD check's "Rd \usage sections" check then flags every one of those
# @param entries as "documented but not in usage" (correctly -- there IS no
# usage for them to be in). The @param text must stay in the source: it's
# what feeds the correctly-generated per-method \usage on the real R6
# subclasses that splice this mixin in. So instead of touching the source,
# strip the resulting (structurally unfixable) \arguments{} block from just
# these bare-object Rd pages after generation -- \name/\alias/\usage/
# \description/\value/\section/\link cross-references are untouched, so nothing
# that links to these pages breaks.
#
# The real fix is teaching fast_roxygenize.R's R6-method-doc patches (see the
# object_defaults.r6class / extract_r6_methods patches above) to also handle
# plain-list mixins, generating a proper per-method \usage the same way it
# already does for R6 classes. Until then, this keeps `R CMD check` clean.
strip_arguments_block_from_bare_object_rd = function(path) {
	lines = readLines(path)
	usage_start = which(lines == "\\usage{")
	if (length(usage_start) != 1L) return(invisible(FALSE))
	usage_end = which(lines == "}")
	usage_end = usage_end[usage_end > usage_start][1]
	if (is.na(usage_end)) return(invisible(FALSE))
	usage_body = trimws(lines[(usage_start + 1L):(usage_end - 1L)])
	usage_body = usage_body[nzchar(usage_body)]
	# A bare object (not a function call) has no usage text containing "(".
	if (length(usage_body) != 1L || grepl("(", usage_body, fixed = TRUE)) return(invisible(FALSE))
	if (!grepl("^[A-Za-z._][A-Za-z0-9._]*$", usage_body)) return(invisible(FALSE))

	args_start = which(lines == "\\arguments{")
	if (length(args_start) != 1L) return(invisible(FALSE))
	depth = 0L
	args_end = NA_integer_
	for (i in args_start:length(lines)) {
		depth = depth + lengths(regmatches(lines[i], gregexpr("\\{", lines[i]))) -
		                lengths(regmatches(lines[i], gregexpr("\\}", lines[i])))
		if (depth == 0L) { args_end = i; break }
	}
	if (is.na(args_end)) return(invisible(FALSE))
	# Also drop the blank line that roxygen2 puts before \arguments{.
	drop_start = if (args_start > 1L && !nzchar(lines[args_start - 1L])) args_start - 1L else args_start
	writeLines(lines[-(drop_start:args_end)], path)
	invisible(TRUE)
}
for (rd_file in list.files(file.path("EDI", "man"), pattern = "\\.Rd$", full.names = TRUE)) {
	strip_arguments_block_from_bare_object_rd(rd_file)
}
