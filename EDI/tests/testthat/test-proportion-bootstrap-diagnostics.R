build_airquality_proportion_dataset <- function(max_n_dataset = 150L){
	aq = stats::na.omit(airquality)
	idx = sample.int(nrow(aq), max_n_dataset, replace = TRUE)
	aqs = aq[idx, , drop = FALSE]

	X = stats::model.matrix(Wind ~ 0 + ., data = aqs)
	X = apply(X, 2, scale)
	X = X / ncol(X)
	X = as.data.frame(X, check.names = FALSE)
	X = X[, !vapply(X, function(col) any(is.na(col)), logical(1)), drop = FALSE]

	y_cont = aqs$Wind
	y_prop = pmin(1, pmax(0, (y_cont - min(y_cont) + 1e-6) / max(y_cont - min(y_cont) + 2e-6)))

	list(X = X, y = y_prop)
}

build_seq_proportion_design <- function(design_name, X, y, beta_T = 0, sd_noise = 0.1){
	des = switch(
		design_name,
		Bernoulli = DesignSeqOneByOneBernoulli$new(response_type = "proportion", n = nrow(X), verbose = FALSE),
		iBCRD = DesignSeqOneByOneiBCRD$new(response_type = "proportion", n = nrow(X), verbose = FALSE),
		Efron = DesignSeqOneByOneEfron$new(response_type = "proportion", n = nrow(X), verbose = FALSE),
		stop("Unsupported design_name: ", design_name)
	)

	for (i in seq_len(nrow(X))) {
		w_t = des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
		eps = stats::rnorm(1, 0, sd_noise)
		bt = ifelse(w_t == 1, beta_T, 0)
		y_t = pmin(1, pmax(0, y[i] + bt + eps))
		des$add_one_subject_response(i, y_t, 1)
	}

	des
}

collect_bootstrap_diagnostics <- function(inf_obj, B = 51L){
	dbg = tryCatch(
		inf_obj$approximate_bootstrap_distribution_beta_hat_T(B = B, show_progress = FALSE, debug = TRUE),
		error = function(e) e
	)
	ci = tryCatch(
		inf_obj$compute_bootstrap_confidence_interval(B = B, na.rm = TRUE, show_progress = FALSE),
		error = function(e) e
	)
	pv = tryCatch(
		inf_obj$compute_bootstrap_two_sided_pval(B = B, na.rm = TRUE),
		error = function(e) e
	)

	if (inherits(dbg, "error")) {
		return(data.frame(
			debug_error = TRUE,
			debug_error_message = conditionMessage(dbg),
			prop_illegal = NA_real_,
			prop_err = NA_real_,
			prop_warn = NA_real_,
			n_finite = NA_integer_,
			n_unique_finite = NA_integer_,
			boot_ci_failed = inherits(ci, "error"),
			boot_ci_degenerate = FALSE,
			boot_ci_lower = NA_real_,
			boot_ci_upper = NA_real_,
			boot_p_failed = inherits(pv, "error") || !is.finite(as.numeric(pv)[1L]),
			boot_p = if (inherits(pv, "error")) NA_real_ else as.numeric(pv)[1L],
			stringsAsFactors = FALSE
		))
	}

	values = as.numeric(dbg$values)
	finite_vals = values[is.finite(values)]
	ci_vals = if (inherits(ci, "error")) c(NA_real_, NA_real_) else as.numeric(ci)
	pv_val = if (inherits(pv, "error")) NA_real_ else as.numeric(pv)[1L]

	data.frame(
		debug_error = FALSE,
		debug_error_message = "",
		prop_illegal = as.numeric(dbg$prop_illegal_values),
		prop_err = as.numeric(dbg$prop_iterations_with_errors),
		prop_warn = as.numeric(dbg$prop_iterations_with_warnings),
		n_finite = length(finite_vals),
		n_unique_finite = length(unique(signif(finite_vals, 8))),
		boot_ci_failed = inherits(ci, "error") || any(!is.finite(ci_vals)),
		boot_ci_degenerate = all(is.finite(ci_vals)) && isTRUE(all.equal(ci_vals[1L], ci_vals[2L])),
		boot_ci_lower = ci_vals[1L],
		boot_ci_upper = ci_vals[2L],
		boot_p_failed = inherits(pv, "error") || !is.finite(pv_val),
		boot_p = pv_val,
		stringsAsFactors = FALSE
	)
}

test_that("proportion bootstrap diagnostics explain any comprehensive-style failures", {
	set.seed(1)

	class_names = c(
		"InferencePropUniGCompMeanDiff",
		"InferencePropMultiGCompMeanDiff",
		"InferencePropUniZeroOneInflatedBetaRegr",
		"InferencePropMultiZeroOneInflatedBetaRegr"
	)
	design_names = c("Bernoulli", "iBCRD", "Efron")
	B = 51L

	dat = build_airquality_proportion_dataset(max_n_dataset = 150L)
	rows = list()

	for (design_name in design_names) {
		des = build_seq_proportion_design(design_name, dat$X, dat$y, beta_T = 0, sd_noise = 0.1)
		for (class_name in class_names) {
			inf_obj = get(class_name, inherits = TRUE)$new(des, verbose = FALSE)
			diag_row = collect_bootstrap_diagnostics(inf_obj, B = B)
			diag_row$inference_class = class_name
			diag_row$design = design_name
			rows[[length(rows) + 1L]] = diag_row
		}
	}

	diagnostics = do.call(rbind, rows)
	diagnostics = diagnostics[, c(
		"inference_class", "design",
		"debug_error", "debug_error_message",
		"prop_illegal", "prop_err", "prop_warn",
		"n_finite", "n_unique_finite",
		"boot_ci_failed", "boot_ci_degenerate", "boot_ci_lower", "boot_ci_upper",
		"boot_p_failed", "boot_p"
	)]

	message("Comprehensive-style bootstrap diagnostics (B = ", B, "):")
	print(diagnostics)

	expect_equal(nrow(diagnostics), length(class_names) * length(design_names))
	expect_false(any(diagnostics$debug_error))

	failed = diagnostics[diagnostics$boot_ci_failed | diagnostics$boot_ci_degenerate | diagnostics$boot_p_failed, , drop = FALSE]
	if (nrow(failed) > 0L) {
		expect_true(all(
			failed$prop_illegal > 0 |
			failed$prop_err > 0 |
			failed$n_finite < B |
			failed$n_unique_finite <= 1
		))
	}
})
