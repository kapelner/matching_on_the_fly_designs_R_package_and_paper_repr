#!/usr/bin/env Rscript

`%||%` = function(x, y) if (is.null(x)) y else x

valid_invariant_statuses = c("ok", "failed", "skipped", "not_applicable")

invariant_result = function(invariant, status, reason) {
	if (!status %in% valid_invariant_statuses) {
		stop("Unknown invariant status: ", status, call. = FALSE)
	}
	data.frame(
		invariant = as.character(invariant),
		status = as.character(status),
		reason = as.character(reason),
		stringsAsFactors = FALSE
	)
}

split_invariant_expr = function(invariants_expr) {
	if (is.null(invariants_expr) || length(invariants_expr) == 0L || is.na(invariants_expr[1L])) return(character())
	invariants = strsplit(as.character(invariants_expr[1L]), ";", fixed = TRUE)[[1L]]
	invariants[nzchar(invariants)]
}

is_documented_nonestimable = function(context) {
	isTRUE(context$documented_nonestimable_behavior) ||
		isTRUE(context$documented_nonestimable_ok) ||
		isTRUE(context$allow_documented_missing)
}

is_missing_numeric_output = function(output) {
	is.numeric(output) && length(output) > 0L && all(is.na(output))
}

invariant_no_unexpected_error = function(output, execution_status, context = list()) {
	if (identical(execution_status, "ok")) {
		return(invariant_result("no_unexpected_error", "ok", "public call completed without unexpected error"))
	}
	invariant_result("no_unexpected_error", "failed", paste("public call status was", execution_status))
}

invariant_numeric_finite_or_documented_missing = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("numeric_finite_or_documented_missing", "skipped", "execution did not complete"))
	if (!is.numeric(output)) return(invariant_result("numeric_finite_or_documented_missing", "failed", "output is not numeric"))
	if (length(output) == 0L) return(invariant_result("numeric_finite_or_documented_missing", "failed", "numeric output is empty"))
	if (all(is.finite(output))) return(invariant_result("numeric_finite_or_documented_missing", "ok", "numeric output is finite"))
	if (is_missing_numeric_output(output) && is_documented_nonestimable(context)) {
		return(invariant_result("numeric_finite_or_documented_missing", "ok", "missing numeric output is documented as non-estimable"))
	}
	invariant_result("numeric_finite_or_documented_missing", "failed", "numeric output contains non-finite or undocumented missing values")
}

invariant_pval_in_0_1 = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("pval_in_0_1", "skipped", "execution did not complete"))
	if (!is.numeric(output) || length(output) < 1L) return(invariant_result("pval_in_0_1", "failed", "p-value output is not numeric"))
	if (is_missing_numeric_output(output) && is_documented_nonestimable(context)) return(invariant_result("pval_in_0_1", "ok", "missing p-value is documented"))
	bad = is.na(output) | output < 0 | output > 1
	if (any(bad)) return(invariant_result("pval_in_0_1", "failed", "p-value is outside [0, 1] or missing"))
	invariant_result("pval_in_0_1", "ok", "p-value is within [0, 1]")
}

invariant_ci_length_2 = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("ci_length_2", "skipped", "execution did not complete"))
	if (!is.numeric(output)) return(invariant_result("ci_length_2", "failed", "confidence interval output is not numeric"))
	if (length(output) != 2L) return(invariant_result("ci_length_2", "failed", "confidence interval does not have length 2"))
	invariant_result("ci_length_2", "ok", "confidence interval has length 2")
}

invariant_ci_ordered = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("ci_ordered", "skipped", "execution did not complete"))
	if (!is.numeric(output) || length(output) != 2L) return(invariant_result("ci_ordered", "failed", "confidence interval is not numeric length 2"))
	if (all(is.na(output)) && is_documented_nonestimable(context)) return(invariant_result("ci_ordered", "ok", "missing confidence interval is documented"))
	if (any(is.na(output))) return(invariant_result("ci_ordered", "failed", "confidence interval contains undocumented missing values"))
	if (output[1L] > output[2L]) return(invariant_result("ci_ordered", "failed", "confidence interval lower endpoint exceeds upper endpoint"))
	invariant_result("ci_ordered", "ok", "confidence interval endpoints are ordered")
}

invariant_estimate_shape_valid = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("estimate_shape_valid", "skipped", "execution did not complete"))
	if (!is.numeric(output) || length(output) < 1L) return(invariant_result("estimate_shape_valid", "failed", "estimate output is not a non-empty numeric vector"))
	if (length(output) > 1L && is.null(names(output))) return(invariant_result("estimate_shape_valid", "ok", "estimate is numeric vector without documented names"))
	invariant_result("estimate_shape_valid", "ok", "estimate output has valid numeric shape")
}

invariant_distribution_length_valid = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("distribution_length_valid", "skipped", "execution did not complete"))
	if (!is.numeric(output)) return(invariant_result("distribution_length_valid", "failed", "distribution output is not numeric"))
	if (length(output) < 1L) return(invariant_result("distribution_length_valid", "failed", "distribution output is empty"))
	expected_length = context$expected_distribution_length %||% context$B %||% NULL
	if (!is.null(expected_length) && !is.na(expected_length) && length(output) != as.integer(expected_length)) {
		return(invariant_result("distribution_length_valid", "failed", "distribution length does not match documented expected length"))
	}
	invariant_result("distribution_length_valid", "ok", "distribution output has valid length")
}

invariant_documented_nonestimable_ok = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("documented_nonestimable_ok", "skipped", "execution did not complete"))
	if (!is_missing_numeric_output(output)) return(invariant_result("documented_nonestimable_ok", "not_applicable", "output is estimable or not represented as all-NA numeric"))
	if (is_documented_nonestimable(context)) return(invariant_result("documented_nonestimable_ok", "ok", "non-estimable output is documented"))
	invariant_result("documented_nonestimable_ok", "failed", "missing numeric output is not documented as non-estimable")
}

invariant_seeded_call_deterministic = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("seeded_call_deterministic", "skipped", "execution did not complete"))
	if (is.null(context$repeated_output)) return(invariant_result("seeded_call_deterministic", "skipped", "no repeated seeded output supplied"))
	if (identical(output, context$repeated_output)) return(invariant_result("seeded_call_deterministic", "ok", "seeded repeated output is identical"))
	invariant_result("seeded_call_deterministic", "failed", "seeded repeated output differs")
}

invariant_output_names_stable_when_documented = function(output, execution_status, context = list()) {
	if (!identical(execution_status, "ok")) return(invariant_result("output_names_stable_when_documented", "skipped", "execution did not complete"))
	expected_names = context$expected_names %||% NULL
	if (is.null(expected_names)) return(invariant_result("output_names_stable_when_documented", "not_applicable", "no documented output names supplied"))
	actual_names = names(output)
	if (identical(actual_names, expected_names)) return(invariant_result("output_names_stable_when_documented", "ok", "output names match documented names"))
	invariant_result("output_names_stable_when_documented", "failed", "output names do not match documented names")
}

public_argument_invariant_functions = list(
	no_unexpected_error = invariant_no_unexpected_error,
	numeric_finite_or_documented_missing = invariant_numeric_finite_or_documented_missing,
	pval_in_0_1 = invariant_pval_in_0_1,
	ci_length_2 = invariant_ci_length_2,
	ci_ordered = invariant_ci_ordered,
	estimate_shape_valid = invariant_estimate_shape_valid,
	distribution_length_valid = invariant_distribution_length_valid,
	documented_nonestimable_ok = invariant_documented_nonestimable_ok,
	seeded_call_deterministic = invariant_seeded_call_deterministic,
	output_names_stable_when_documented = invariant_output_names_stable_when_documented
)

evaluate_public_argument_invariants = function(output, execution_status = "ok", invariants_expr = "", context = list()) {
	invariants = split_invariant_expr(invariants_expr)
	if (!length(invariants)) {
		return(invariant_result("no_registered_invariants", "not_applicable", "no invariants registered for case"))
	}
	rows = lapply(invariants, function(invariant) {
		fn = public_argument_invariant_functions[[invariant]]
		if (is.null(fn)) {
			return(invariant_result(invariant, "skipped", "unknown invariant name"))
		}
		fn(output, execution_status, context)
	})
	out = do.call(rbind, rows)
	row.names(out) = NULL
	out
}

summarize_public_argument_invariants = function(results) {
	if (!is.data.frame(results) || !nrow(results)) return("not_checked")
	failed = results[results$status == "failed", , drop = FALSE]
	if (nrow(failed)) return(paste0("failed:", paste(failed$invariant, collapse = ";")))
	if (any(results$status == "ok")) return("ok")
	if (all(results$status == "skipped")) return("skipped")
	"not_applicable"
}

check_invariants = function(output, status, invariants_expr, context = list()) {
	summarize_public_argument_invariants(
		evaluate_public_argument_invariants(output, status, invariants_expr, context)
	)
}
