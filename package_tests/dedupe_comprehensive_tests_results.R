args = commandArgs(trailingOnly = TRUE)

if (!requireNamespace("data.table", quietly = TRUE)) {
	stop("The data.table package is required.")
}

files = if (length(args) == 0L) {
	Sys.glob("package_tests/comprehensive_tests_results_nc_*.csv")
} else {
	args
}

if (length(files) == 0L) {
	stop("No comprehensive test result files were found.")
}

key_cols = c(
	"rep",
	"beta_T",
	"dataset",
	"response_type",
	"design",
	"inference_class",
	"function_run"
)

for (file in files) {
	if (!file.exists(file)) {
		message("Skipping missing file: ", file)
		next
	}

	dt = data.table::fread(file)
	if (nrow(dt) == 0L) {
		message("Skipping empty file: ", file)
		next
	}

	missing_cols = setdiff(key_cols, names(dt))
	if (length(missing_cols) > 0L) {
		stop("File ", file, " is missing required columns: ", paste(missing_cols, collapse = ", "))
	}

	order_cols = intersect(c("run_row_id"), names(dt))
	if (length(order_cols) > 0L) {
		data.table::setorderv(dt, order_cols)
	}

	before_n = nrow(dt)
	dt = dt[, .SD[.N], by = key_cols]
	if ("run_row_id" %in% names(dt)) {
		data.table::setorderv(dt, "run_row_id")
	}

	data.table::fwrite(dt, file)
	message(
		"Deduped ", file,
		": kept ", nrow(dt),
		" of ", before_n,
		" rows."
	)
}
