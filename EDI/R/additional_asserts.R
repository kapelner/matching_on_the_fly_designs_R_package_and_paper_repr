# Internal cross-argument assertions promoted from the public argument
# combination runner. These helpers intentionally remain unexported.

assertBootstrapArgs = function(B, min_number_usable_samples = NULL, type = NULL, allowed_types = NULL, context = "bootstrap"){
	if (!should_run_asserts()) return(invisible(NULL))
	assertCount(B, positive = TRUE)
	if (!is.null(min_number_usable_samples)) {
		assertCount(min_number_usable_samples, positive = TRUE)
		if (as.integer(min_number_usable_samples) > as.integer(B)) {
			stop(
				"min_number_usable_samples must be less than or equal to B for ",
				context,
				", but min_number_usable_samples = ",
				as.integer(min_number_usable_samples),
				" and B = ",
				as.integer(B),
				".",
				call. = FALSE
			)
		}
	}
	if (!is.null(type) && !is.null(allowed_types)) {
		assertChoice(type, allowed_types)
	}
	invisible(NULL)
}

assertFormulaContext = function(model_formula, data, context = "model_formula"){
	if (!should_run_asserts()) return(invisible(NULL))
	assertFormula(model_formula, null.ok = TRUE)
	if (is.null(model_formula) || is.null(data)) return(invisible(NULL))
	all_features = names(data)
	required_vars = setdiff(all.vars(model_formula), ".")
	missing_vars = setdiff(required_vars, all_features)
	if (length(missing_vars) > 0L) {
		stop(
			context,
			" contains variables not present in the available covariates: ",
			paste(missing_vars, collapse = ", "),
			".",
			call. = FALSE
		)
	}
	invisible(NULL)
}

assertStrataClusterArgs = function(strata_cols = NULL, cluster_col = NULL, data = NULL, strata_cols_must_be_factor = FALSE, context = "strata/cluster arguments"){
	if (!should_run_asserts()) return(invisible(NULL))
	if (!is.null(strata_cols)) assertCharacter(strata_cols, min.len = 1, any.missing = FALSE)
	if (!is.null(cluster_col)) assertCharacter(cluster_col, len = 1, any.missing = FALSE)
	if (!is.null(strata_cols) && !is.null(cluster_col) && cluster_col %in% strata_cols) {
		stop("cluster_col must not also appear in strata_cols for ", context, ".", call. = FALSE)
	}
	if (is.null(data)) return(invisible(NULL))
	data_names = names(data)
	needed_cols = c(strata_cols %||% character(), cluster_col %||% character())
	missing_cols = setdiff(needed_cols, data_names)
	if (length(missing_cols) > 0L) {
		stop(
			context,
			" references column(s) not present in the supplied covariates: ",
			paste(missing_cols, collapse = ", "),
			".",
			call. = FALSE
		)
	}
	if (isTRUE(strata_cols_must_be_factor) && length(strata_cols)) {
		is_categorical = vapply(strata_cols, function(col) is.factor(data[[col]]) || is.character(data[[col]]), logical(1))
		non_categorical_cols = strata_cols[!is_categorical]
		if (length(non_categorical_cols) > 0L) {
			stop(
				"strata_cols must refer to factor/categorical columns for ",
				context,
				"; non-categorical column(s): ",
				paste(non_categorical_cols, collapse = ", "),
				".",
				call. = FALSE
			)
		}
	}
	invisible(NULL)
}
