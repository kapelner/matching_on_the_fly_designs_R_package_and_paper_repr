#libreoffice --calc /home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/package_tests/calc_external_links/comprehensive_tests_results_nc_1_links.fods
rm(list = ls())
set.seed(1)
devtools::load_all("EDI")
pacman::p_load(doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
max_n_dataset = 148 #needs to be divisible by 4 for some blocking designs
source("package_tests/_dataset_load.R")
# options(error = recover)
# options(warn=2)

args = commandArgs(trailingOnly = TRUE)
Nrep      = if (length(args) >= 1) as.integer(args[1]) else 40L
NUM_CORES = if (length(args) >= 2) as.integer(args[2]) else 2L
ALL_RESPONSE_TYPES = c("continuous", "incidence", "proportion", "count", "survival", "ordinal")
RESPONSE_TYPE_FILTER = if (length(args) >= 3) as.character(args[3]) else NA_character_
if (!is.na(RESPONSE_TYPE_FILTER) && !(RESPONSE_TYPE_FILTER %in% ALL_RESPONSE_TYPES)) {
	stop(
		"Unsupported response_type filter: ",
		RESPONSE_TYPE_FILTER,
		". Supported values are: ",
		paste(ALL_RESPONSE_TYPES, collapse = ", ")
	)
}
set_num_cores(NUM_CORES)
prob_censoring = 0.15
r = 351
pval_epsilon = 0.007
test_compute_confidence_interval_rand = TRUE
beta_T_values = c(0, 1)
SD_NOISE = 0.1
pending_rep_header = NULL
pending_beta_header = NULL
pending_dataset_header = NULL
pending_response_header = NULL
pending_design_header = NULL
pending_banner = NULL

results_file = if (is.na(RESPONSE_TYPE_FILTER)) {
	paste0("package_tests/comprehensive_tests_results_nc_", NUM_CORES, ".csv")
} else {
	paste0("package_tests/comprehensive_tests_results_nc_", NUM_CORES, "_", RESPONSE_TYPE_FILTER, ".csv")
}
existing_results_dt = if (file.exists(results_file)) data.table::fread(results_file) else data.table::data.table()
run_row_id = if ("run_row_id" %in% colnames(existing_results_dt) && nrow(existing_results_dt) > 0L) {
	as.integer(max(existing_results_dt$run_row_id, na.rm = TRUE))
} else {
	0L
}

serialize_beta_T = function(value){
	if (is.na(value)) return("NA")
	if (is.numeric(value)) return(format(as.numeric(value), scientific = TRUE, digits = 17))
	as.character(value)
}

round_result_field = function(value){
	if (is.null(value) || is.na(value)) return(NA_character_)
	if (value == "") return("")
	num = suppressWarnings(as.numeric(value))
	if (!is.finite(num)) return(value)
	sprintf("%.3f", num)
}

round_duration_field = function(value){
	if (is.null(value) || is.na(value)) return(NA_real_)
	num = suppressWarnings(as.numeric(value))
	if (!is.finite(num)) return(num)
	round(num, 3)
}

add_assignment_only_cluster_id = function(X_design, strata_cols = character(0), cluster_size = 2L){
	X_out = as.data.frame(X_design)
	cluster_col = ".assignment_only_cluster_id"
	cluster_ids = integer(nrow(X_out))
	next_cluster_id = 1L

	if (!length(strata_cols)) {
		cluster_ids = ((seq_len(nrow(X_out)) - 1L) %/% cluster_size) + 1L
	} else {
		strata_df = X_out[, strata_cols, drop = FALSE]
		strata_key = do.call(paste, c(lapply(strata_df, as.character), sep = "\r"))
		for (key in unique(strata_key)) {
			idx = which(strata_key == key)
			cluster_ids[idx] = ((seq_along(idx) - 1L) %/% cluster_size) + next_cluster_id
			next_cluster_id = max(cluster_ids[idx]) + 1L
		}
	}

	X_out[[cluster_col]] = factor(cluster_ids)
	list(X = X_out, cluster_col = cluster_col)
}

build_result_key = function(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val){
	paste(
		as.integer(rep_val),
		serialize_beta_T(beta_val),
		as.character(dataset_val),
		as.character(response_val),
		as.character(design_val),
		as.character(inference_val),
		as.character(function_run_val),
		sep = "||"
	)
}

completed_rows_cache = new.env(parent = emptyenv())

mark_row_completed = function(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val){
	key = build_result_key(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val)
	completed_rows_cache[[key]] <- TRUE
}

last_results_mtime = NULL

reload_completed_rows = function() {
	if (file.exists(results_file) && file.info(results_file)$size > 0) {
		current_mtime = file.info(results_file)$mtime
		if (!is.null(last_results_mtime) && current_mtime == last_results_mtime) {
			return(invisible(NULL))
		}
		
		# Reset the cache
		completed_rows_cache <<- new.env(parent = emptyenv())
		
		# Define columns needed for building the key and filtering status
		needed_cols = c("rep", "beta_T", "dataset", "response_type", "design", "inference_class", "function_run", "status")
		
		# Try to read the file, handling potential locks or partial writes with retries
		max_retries = 10
		retry_count = 0
		dt = NULL
		
		while (retry_count < max_retries) {
			dt = tryCatch({
				# Check which columns actually exist first
				header = names(data.table::fread(results_file, nrows = 0))
				cols_to_read = intersect(needed_cols, header)
				data.table::fread(results_file, select = cols_to_read)
			}, error = function(e) {
				NULL
			})
			
			if (!is.null(dt)) break
			
			retry_count = retry_count + 1
			if (retry_count < max_retries) {
				Sys.sleep(0.5)
			}
		}
		
		if (is.null(dt)) {
			# If still failing after retries, return early (cache remains empty)
			return(invisible(NULL))
		}
		
		last_results_mtime <<- current_mtime
		
		if (nrow(dt) > 0L) {
			rows_to_cache = if ("status" %in% colnames(dt)) {
				dt[status == "ok"]
			} else {
				dt
			}
			for (row_idx in seq_len(nrow(rows_to_cache))) {
				row = rows_to_cache[row_idx]
				mark_row_completed(
					row$rep,
					row$beta_T,
					row$dataset,
					row$response_type,
					row$design,
					row$inference_class,
					row$function_run
				)
			}
		}
	}
}

is_row_completed = function(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val){
	# Reload from disk each time we check
	reload_completed_rows()
	
	key = build_result_key(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val)
	!is.null(completed_rows_cache[[key]])
}

if (nrow(existing_results_dt) > 0L) {
	rows_to_cache = if ("status" %in% colnames(existing_results_dt)) {
		existing_results_dt[status == "ok"]
	} else {
		existing_results_dt
	}
	for (row_idx in seq_len(nrow(rows_to_cache))) {
		row = rows_to_cache[row_idx]
		mark_row_completed(
			row$rep,
			row$beta_T,
			row$dataset,
			row$response_type,
			row$design,
			row$inference_class,
			row$function_run
		)
	}
}
results_dt = data.table(
	rep = integer(),
	beta_T = numeric(),
	dataset = character(),
	response_type = character(),
	design = character(),
	inference_class = character(),
	function_run = character(),
	timestamp = character(),
	duration_time_sec = numeric(),
	result_1 = character(),
	result_2 = character(),
	beta_T_in_confidence_interval = logical(),
	error_message = character(),
	run_row_id = integer(),
	r = integer(),
	pval_epsilon = numeric(),
	prob_censoring = numeric(),
	sd_noise = numeric(),
	num_cores = integer(),
	dataset_n_rows = integer(),
	dataset_n_cols = integer(),
	result = character(),
	status = character()
)
write_results_if_needed = function(force = FALSE){
	if ((force || nrow(results_dt) > 0L) && nrow(results_dt) > 0L){
		append_mode = file.exists(results_file) && file.info(results_file)$size > 0
		data.table::fwrite(
			results_dt,
			results_file,
			append = append_mode,
			col.names = !append_mode,
			na = "NA"
		)
		results_dt <<- results_dt[0]
	}
}

log_progress = function(msg){
	message(msg)
	flush.console()
}

record_result = function(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, inference_class, function_run, result, status, duration_time_sec, error_message = NA_character_){
	result_vec = if (is.null(result)) {
		NA_character_
	} else if (length(result) == 0) {
		character(0)
	} else {
		as.character(result)
	}

	result_str = if (is.null(result)) {
		NA_character_
	} else if (length(result) == 0) {
		""
	} else if (length(result) == 1) {
		as.character(result)
	} else {
		paste(as.character(result), collapse = " ")
	}
	result_1 = if (length(result_vec) >= 1) round_result_field(result_vec[1]) else NA_character_
	is_debug_distribution = grepl("_debug$", function_run)
	result_2 = if ((grepl("confidence_interval", function_run, fixed = TRUE) || is_debug_distribution) && length(result_vec) >= 2) round_result_field(result_vec[2]) else NA_character_
	beta_T_in_confidence_interval = NA
	if (grepl("confidence_interval", function_run, fixed = TRUE) && length(result) >= 2 && all(is.finite(result[1:2]))){
		ci_lo = min(result[1:2])
		ci_hi = max(result[1:2])
		beta_T_in_confidence_interval = (beta_T >= ci_lo && beta_T <= ci_hi)
	}
	run_row_id <<- run_row_id + 1L
	results_dt <<- data.table::rbindlist(list(
		results_dt,
		data.table(
			rep = as.integer(rep_curr),
			beta_T = beta_T,
			dataset = dataset_name,
			response_type = response_type,
			design = design_type,
			inference_class = inference_class,
			function_run = function_run,
			timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
			duration_time_sec = round_duration_field(duration_time_sec),
			result_1 = result_1,
			result_2 = result_2,
			beta_T_in_confidence_interval = beta_T_in_confidence_interval,
			error_message = ifelse(is.null(error_message), NA_character_, as.character(error_message)),
			run_row_id = run_row_id,
			r = as.integer(r),
			pval_epsilon = pval_epsilon,
			prob_censoring = prob_censoring,
			sd_noise = SD_NOISE,
			num_cores = as.integer(NUM_CORES),
			dataset_n_rows = as.integer(dataset_n_rows),
			dataset_n_cols = as.integer(dataset_n_cols),
			result = result_str,
			status = status
		)
	), use.names = TRUE)
	if (identical(status, "ok")) {
		mark_row_completed(rep_curr, beta_T, dataset_name, response_type, design_type, inference_class, function_run)
	}
	write_results_if_needed(force = TRUE)
}

run_inference_checks = function(seq_des_inf, response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols){
	skip_slow = FALSE
	skip_bootstrap = is(seq_des_inf, "InferenceAbstractKKGEE") || is(seq_des_inf, "InferenceAbstractKKGLMM") || is(seq_des_inf, "InferenceContinMultGLS") || is(seq_des_inf, "InferenceAbstractKKClaytonCopulaIVWC") || is(seq_des_inf, "InferenceAbstractKKClaytonCopulaCombinedLikelihood") || is(seq_des_inf, "InferenceAbstractKKWeibullFrailtyIVWC") || is(seq_des_inf, "InferenceAbstractKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "InferenceAllKKWilcoxIVWC") || is(seq_des_inf, "InferenceAbstractKKWilcoxRegrIVWC") || is(seq_des_inf, "InferenceSurvivalUnivKKRankRegrIVWC") || is(seq_des_inf, "InferenceIncidExactZhangAbstract") || is(seq_des_inf, "InferenceAllSimpleWilcox") || is(seq_des_inf, "InferenceOrdinalUnivKKGEE") || is(seq_des_inf, "InferenceOrdinalUnivKKGLMM") || is(seq_des_inf, "InferenceOrdinalMultiKKGLMM") || is(seq_des_inf, "InferenceOrdinalUnivKKGLMMProbit") || is(seq_des_inf, "InferenceOrdinalMultiKKGLMMProbit") || is(seq_des_inf, "InferenceOrdinalPairedSignTest") || is(seq_des_inf, "InferenceOrdinalUnivKKCondPropOddsCombinedRegr") || is(seq_des_inf, "InferenceOrdinalUnivKKCondContRatioRegr") || is(seq_des_inf, "InferenceOrdinalUnivKKCondAdjCatLogitRegr") || is(seq_des_inf, "InferenceOrdinalUniGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalMultiGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalMultiCLLRegr") || is(seq_des_inf, "InferenceOrdinalUniOrderedProbitRegr") || is(seq_des_inf, "InferenceOrdinalMultiOrderedProbitRegr") || is(seq_des_inf, "InferenceOrdinalUniCauchitRegr") || is(seq_des_inf, "InferenceOrdinalMultiCauchitRegr") || is(seq_des_inf, "InferenceOrdinalMultiContRatioRegr") || is(seq_des_inf, "InferenceOrdinalMultiKKGEE") || is(seq_des_inf, "InferenceOrdinalMultiKKCondContRatioRegr") || is(seq_des_inf, "InferenceOrdinalMultiKKCondAdjCatLogitRegr")
	skip_rand      = is(seq_des_inf, "InferenceAbstractKKGEE") || is(seq_des_inf, "InferenceAbstractKKGLMM") || is(seq_des_inf, "InferenceIncidExactZhangAbstract") || is(seq_des_inf, "InferencePropUniGCompMeanDiff") || is(seq_des_inf, "InferencePropMultiGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalUnivKKGEE") || is(seq_des_inf, "InferenceOrdinalUnivKKGLMM") || is(seq_des_inf, "InferenceOrdinalMultiKKGLMM") || is(seq_des_inf, "InferenceOrdinalUnivKKGLMMProbit") || is(seq_des_inf, "InferenceOrdinalMultiKKGLMMProbit") || is(seq_des_inf, "InferenceOrdinalPairedSignTest") || is(seq_des_inf, "InferenceOrdinalUnivKKCondPropOddsCombinedRegr") || is(seq_des_inf, "InferenceOrdinalUnivKKCondContRatioRegr") || is(seq_des_inf, "InferenceOrdinalUnivKKCondAdjCatLogitRegr") || is(seq_des_inf, "InferenceOrdinalUniGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalMultiGCompMeanDiff") || is(seq_des_inf, "InferenceOrdinalMultiCLLRegr") || is(seq_des_inf, "InferenceOrdinalUniOrderedProbitRegr") || is(seq_des_inf, "InferenceOrdinalMultiOrderedProbitRegr") || is(seq_des_inf, "InferenceOrdinalUniCauchitRegr") || is(seq_des_inf, "InferenceOrdinalMultiCauchitRegr") || is(seq_des_inf, "InferenceOrdinalMultiContRatioRegr") || is(seq_des_inf, "InferenceOrdinalMultiKKGEE") || is(seq_des_inf, "InferenceOrdinalMultiKKCondContRatioRegr") || is(seq_des_inf, "InferenceOrdinalMultiKKCondAdjCatLogitRegr")
	skip_mle_pval  = is(seq_des_inf, "InferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood")
	skip_rand_pval = is(seq_des_inf, "InferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "InferenceContinMultGLS") || is(seq_des_inf, "InferencePropUniGCompMeanDiff") || is(seq_des_inf, "InferencePropMultiGCompMeanDiff")
	skip_ci_rand   = is(seq_des_inf, "InferenceContinMultKKQuantileRegrIVWC") || is(seq_des_inf, "InferencePropMultiKKQuantileRegrIVWC") || is(seq_des_inf, "InferenceContinMultKKQuantileRegrCombinedLikelihood") || is(seq_des_inf, "InferencePropMultiKKQuantileRegrCombinedLikelihood") || is(seq_des_inf, "InferencePropUniGCompMeanDiff") || is(seq_des_inf, "InferencePropMultiGCompMeanDiff") || (response_type != "continuous" && (is(seq_des_inf, "InferenceAllSimpleMeanDiff") || is(seq_des_inf, "InferenceAllKKCompoundMeanDiff")))
	skip_ci_rand_custom = FALSE
	
	skip_ci = beta_T == 1 && (
		is(seq_des_inf, "InferenceIncidMultiLogRegr") ||
		is(seq_des_inf, "InferencePropUniBetaRegr") ||
		is(seq_des_inf, "InferencePropMultiBetaRegr") ||
		is(seq_des_inf, "InferenceSurvivalUniCoxPHRegr") ||
		is(seq_des_inf, "InferenceSurvivalMultiCoxPHRegr") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKLWACoxIVWC") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKStratCoxIVWC") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKClaytonCopulaIVWC") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKLWACoxCombinedLikelihood") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKStratCoxCombinedLikelihood") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKWeibullFrailtyIVWC") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood") ||
		is(seq_des_inf, "InferenceAllKKWilcoxRegrMultiIVWC") ||
		is(seq_des_inf, "InferenceSurvivalMultiKKRankRegrIVWC")
	)
	snap_small_numeric_to_zero = function(x, tol = sqrt(.Machine$double.eps)){
		if (is.null(x)) return(x)
		if (is.list(x)) return(lapply(x, snap_small_numeric_to_zero, tol = tol))
		if (is.atomic(x) && is.numeric(x)){
			x[is.finite(x) & abs(x) < tol] = 0
			return(x)
		}
		x
	}

	has_invalid_numeric = function(x){
		if (is.null(x)) return(FALSE)
		if (is.list(x)) return(any(vapply(x, has_invalid_numeric, logical(1))))
		if (is.atomic(x) && is.numeric(x)) return(any(!is.finite(x) | is.na(x) | is.nan(x)))
		FALSE
	}

	is_zero_zero_confidence_interval = function(label, result){
		if (!grepl("confidence_interval", label, fixed = TRUE)) return(FALSE)
		if (!(is.atomic(result) && is.numeric(result))) return(FALSE)
		if (length(result) < 2) return(FALSE)
		isTRUE(all(result[1:2] == 0))
	}

	is_allowed_missing_output = function(label, result){
		if (!has_invalid_numeric(result)) return(FALSE)
		identical(response_type, "ordinal") &&
			identical(label, "compute_asymp_two_sided_pval_for_treatment_effect")
	}

safe_call = function(label, expr){
	if (is_row_completed(
		rep_curr,
		beta_T,
		dataset_name,
		response_type,
		design_type,
		class(seq_des_inf)[1],
		label
	)) {
		return(invisible(NULL))
	}
	
	if (!is.null(pending_rep_header)) { message(pending_rep_header); pending_rep_header <<- NULL }
	if (!is.null(pending_beta_header)) { message(pending_beta_header); pending_beta_header <<- NULL }
	if (!is.null(pending_dataset_header)) { message(pending_dataset_header); pending_dataset_header <<- NULL }
	if (!is.null(pending_response_header)) { message(pending_response_header); pending_response_header <<- NULL }
	if (!is.null(pending_design_header)) { message(pending_design_header); pending_design_header <<- NULL }
	if (!is.null(pending_banner)){
		message(pending_banner)
		pending_banner <<- NULL
	}

	message("          Calling ", label, "()")
	start_elapsed = unname(proc.time()[["elapsed"]])
	tryCatch({
		result <- expr
			if (is_allowed_missing_output(label, result)) {
				message("Recording missing output for ", label, " as ok (ordinal asymptotic p-value not estimable).")
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, NA_character_, status = "ok", duration_time_sec = duration_time_sec)
				return(invisible(NULL))
			}
			if (has_invalid_numeric(result)) {
				msg = paste0("Invalid output detected (NA/NaN/Inf) in ", label)
				message("Skipping ", label, " (non-fatal): ", msg)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = msg)
				return(invisible(NULL))
			}
			if (is_zero_zero_confidence_interval(label, result)) {
				msg = paste0("Degenerate confidence interval [0, 0] detected in ", label)
				message("Skipping ", label, " (non-fatal): ", msg)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = msg)
				return(invisible(NULL))
			}
			result = snap_small_numeric_to_zero(result)
			cat("            ", paste(format(result, digits = 3), collapse = " "), "\n")
			duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
			cat(sprintf("              (Duration: %.3gs)\n", duration_time_sec))
			record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, result, status = "ok", duration_time_sec = duration_time_sec)
			result
		}, error = function(e){
			msg = if (length(e$message) == 0L) "" else e$message
			is_non_fatal = grepl("not implemented", msg, fixed = TRUE) ||
			                 grepl("must implement", msg, fixed = TRUE) ||
			                 grepl("Exact inference is only supported for exact inference classes.", msg, fixed = TRUE) ||
			                 grepl("singular matrix in 'backsolve'", msg, fixed = TRUE) ||
			                 grepl("G-computation RD: could not compute a finite delta-method standard error.", msg, fixed = TRUE) ||
			                 grepl("G-computation RR: could not compute a finite delta-method standard error.", msg, fixed = TRUE) ||
			                 grepl("G-computation RR: could not compute a finite delta-method confidence interval.", msg, fixed = TRUE) ||
			                 grepl("KK g-computation RR: could not compute a finite delta-method confidence interval.", msg, fixed = TRUE) ||
			                 grepl("G-computation mean difference: could not compute a finite delta-method standard error.", msg, fixed = TRUE) ||
			                 grepl("Zero/one-inflated beta requires y in [0, 1]", msg, fixed = TRUE) ||
			                 grepl("Zhang incidence inference is only supported", msg, fixed = TRUE) ||

					 grepl("This type of inference is only available for incidence", msg, fixed = TRUE) ||
						 grepl("not enough discordant pairs", msg, ignore.case = TRUE) ||
						 grepl("Degenerate confidence interval", msg, fixed = TRUE) ||
						 grepl("inconsistent estimator units", msg, ignore.case = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned NA bounds", msg, fixed = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned non-finite bounds", msg, fixed = TRUE) ||
						 					 						 grepl("Weibull regression failed to converge", msg, fixed = TRUE) ||
				 grepl("Negative binomial regression failed to converge", msg, fixed = TRUE) ||
						 					 						 grepl("Invalid output detected", msg, fixed = TRUE) ||
						 					 						 grepl("missing value where TRUE/FALSE needed", msg, fixed = TRUE) ||
						 					 						 ((grepl("NA/NaN/Inf", msg, fixed = TRUE) || grepl("non-finite standard error", msg, fixed = TRUE) || grepl("could not compute a finite standard error", msg, fixed = TRUE)) &&
						 					 
						 													 (is(seq_des_inf, "InferenceAbstractKKClogitIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKPoissonCPoissonIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKStratCoxIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKWeibullFrailtyIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKWilcoxRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKSurvivalRankRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKGEE") ||
						 													  is(seq_des_inf, "InferenceAbstractKKGLMM") ||
						 													  is(seq_des_inf, "InferenceAbstractKKRobustRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKRobustRegrCombinedLikelihood") ||
						 													  is(seq_des_inf, "InferenceAbstractKKQuantileRegrIVWC") ||
						 													  is(seq_des_inf, "InferenceAbstractKKQuantileRegrCombinedLikelihood") ||
						 													  													  is(seq_des_inf, "InferencePropUniFractionalLogit") ||
						 													  													  is(seq_des_inf, "InferenceIncidUnivRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidMultiRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidUnivGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidMultiGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidUnivGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "InferenceIncidMultiGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "InferenceIncidUnivKKGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidMultiKKGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "InferenceIncidUnivKKGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "InferenceIncidMultiKKGCompRiskRatio") ||
						 													  is(seq_des_inf, "InferenceIncidUnivModifiedPoisson") ||
						 													  is(seq_des_inf, "InferenceIncidMultiModifiedPoisson") ||
						 													  is(seq_des_inf, "InferencePropUniGCompMeanDiff") ||
						 													  is(seq_des_inf, "InferencePropMultiGCompMeanDiff") ||
						 													  is(seq_des_inf, "InferencePropUniZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "InferencePropMultiZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "InferenceContinMultGLS")))
						 													  										if (isTRUE(is_non_fatal)){
				message("Skipping ", label, " (non-fatal): ", e$message)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = e$message)
			} else {
				stop(e$message)
			}
		})
	}

	if (is(seq_des_inf, "InferenceOrdinalJonckheereTerpstraTest")){
		safe_call("compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		safe_call("compute_treatment_estimate", seq_des_inf$compute_treatment_estimate())
		return(invisible(NULL))
	}

	supports_exact_inference = is(seq_des_inf, "InferenceExact") || is(seq_des_inf, "InferenceIncidExactZhang")
	if (response_type == "incidence" && supports_exact_inference){
		safe_call("compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		safe_call("compute_exact_confidence_interval", seq_des_inf$compute_exact_confidence_interval(args_for_type = list(Zhang = list(combination_method = "Fisher", pval_epsilon = pval_epsilon))))
	}

	safe_call("compute_treatment_estimate", seq_des_inf$compute_treatment_estimate())
	if (!skip_mle_pval){
		safe_call("compute_asymp_two_sided_pval_for_treatment_effect", seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect())
	}
	if (!skip_ci){
		safe_call("compute_asymp_confidence_interval", seq_des_inf$compute_asymp_confidence_interval(0.05))
	}
	safe_call_debug = function(label, expr) {
		if (is_row_completed(rep_curr, beta_T, dataset_name, response_type, design_type, class(seq_des_inf)[1], label)) {
			return(invisible(NULL))
		}
		if (!is.null(pending_rep_header)) { message(pending_rep_header); pending_rep_header <<- NULL }
		if (!is.null(pending_beta_header)) { message(pending_beta_header); pending_beta_header <<- NULL }
		if (!is.null(pending_dataset_header)) { message(pending_dataset_header); pending_dataset_header <<- NULL }
		if (!is.null(pending_response_header)) { message(pending_response_header); pending_response_header <<- NULL }
		if (!is.null(pending_design_header)) { message(pending_design_header); pending_design_header <<- NULL }
		if (!is.null(pending_banner)) { message(pending_banner); pending_banner <<- NULL }
		message("          Calling ", label, "()")
		start_elapsed = unname(proc.time()[["elapsed"]])
		debug_result = tryCatch(expr, error = function(e) {
			dur = unname(proc.time()[["elapsed"]]) - start_elapsed
			record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type,
						  class(seq_des_inf)[1], label, NA_character_, status = "error",
						  duration_time_sec = dur, error_message = e$message)
			NULL
		})
		if (is.null(debug_result)) return(invisible(NULL))
		duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
		if (identical(label, "approximate_bootstrap_distribution_beta_hat_T_debug")) {
			stats_vec = c(
				debug_result$prop_illegal_values,
				debug_result$prop_iterations_with_errors,
				debug_result$prop_iterations_with_warnings
			)
			cat(sprintf("            prop_illegal=%.3f  prop_err=%.3f  prop_warn=%.3f\n",
					stats_vec[1], stats_vec[2], stats_vec[3]))
		} else {
			stats_vec = c(
				debug_result$prop_iterations_with_errors,
				debug_result$prop_iterations_with_warnings,
				debug_result$prop_illegal_values
			)
			cat(sprintf("            prop_err=%.3f  prop_warn=%.3f  prop_illegal=%.3f\n",
					stats_vec[1], stats_vec[2], stats_vec[3]))
		}
		cat(sprintf("              (Duration: %.3gs)\n", duration_time_sec))
		record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type,
					  class(seq_des_inf)[1], label, stats_vec, status = "ok",
					  duration_time_sec = duration_time_sec)
	}

	if (!skip_slow && !skip_bootstrap){
		safe_call_debug("approximate_bootstrap_distribution_beta_hat_T_debug",
						seq_des_inf$approximate_bootstrap_distribution_beta_hat_T(B = r, debug = TRUE))
	}
	if (!skip_slow && !skip_ci && !skip_bootstrap){
		safe_call("compute_bootstrap_confidence_interval", seq_des_inf$compute_bootstrap_confidence_interval(B = r, na.rm = TRUE))
	}
	if (!skip_slow && !skip_bootstrap){
		safe_call("compute_bootstrap_two_sided_pval", seq_des_inf$compute_bootstrap_two_sided_pval(B = r, na.rm = TRUE))
	}
	if (!skip_slow && !skip_rand && !skip_rand_pval && response_type %in% c("continuous", "survival", "proportion")){
		safe_call_debug("approximate_randomization_distribution_beta_hat_T_debug",
						seq_des_inf$approximate_randomization_distribution_beta_hat_T(r = r, debug = TRUE))
		safe_call("compute_two_sided_pval_for_treatment_effect_rand", seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(r = r, show_progress = FALSE))
		transform_for_rand = switch(
			response_type,
			continuous = "none",
			proportion = "logit",
			count = "log",
			survival = "log",
			"none"
		)
		delta_for_rand = 0.5
		safe_call("compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)",
				seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(r = r, delta = delta_for_rand, transform_responses = transform_for_rand, show_progress = FALSE))
	}

	if (!skip_slow && !skip_rand && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous", "proportion", "count")){
		safe_call("compute_confidence_interval_rand", seq_des_inf$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE))
	}
	if (response_type != "incidence"){
		seq_des_inf$set_custom_randomization_statistic_function(function(){
			yTs = private$des_obj_priv_int$y[private$des_obj_priv_int$w == 1]
			yCs = private$des_obj_priv_int$y[private$des_obj_priv_int$w == 0]
			(mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
		})
		if (!skip_slow && !skip_rand_pval){
			safe_call("compute_two_sided_pval_for_treatment_effect_rand(custom)", seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(r = r, show_progress = FALSE))
		}
		if (!skip_slow && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous")){
			if (!skip_ci_rand_custom){
				safe_call("compute_confidence_interval_rand(custom)", seq_des_inf$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE))
			} else {
				message("    Skipping compute_confidence_interval_rand(custom) (too slow)")
			}
		}
		seq_des_inf$set_custom_randomization_statistic_function(NULL)
	}
}

run_tests_for_response = function(response_type, design_type, dataset_name){
	inference_banner = function(inf_name){
		pending_banner <<- sprintf("\n\n  == Inference: %s design_type = %s dataset = %s response_type = %s beta_T = [%s] num_cores = [%d] rep = [%d/%d]\n", inf_name, design_type, dataset_name, response_type, format(beta_T), NUM_CORES, rep_curr, Nrep)
	}

	apply_treatment_effect_and_noise = function(y_t, w_t, response_type){
		eps = rnorm(1, 0, SD_NOISE)
		bt = ifelse(w_t == 1, beta_T, 0)
		if (response_type == "continuous") return(y_t + bt + eps)
		if (response_type == "incidence") {
			p_base = if (is.finite(y_t) && y_t >= 0 && y_t <= 1) y_t else stats::plogis(y_t)
			p_base = pmin(0.95, pmax(0.05, p_base))
			p_t = plogis(qlogis(p_base) + bt + eps)
			return(as.numeric(stats::rbinom(1, size = 1, prob = p_t)))
		}
		if (response_type == "proportion"){
			return(pmin(1, pmax(0, y_t + bt + eps)))
		}
		if (response_type == "count"){
			lambda_t = pmax(.Machine$double.eps, y_t * exp(bt + eps))
			return(as.numeric(stats::rpois(1, lambda = lambda_t)))
		}
		if (response_type == "survival") return(pmax(.Machine$double.eps, y_t * exp(bt + eps)))
		if (response_type == "ordinal"){
			# For ordinal, we'll just use a simple shift and re-cut if needed, 
			# but here y_t is already the ordinal level.
			# Let's just do a simple shift and round.
			# We use a larger noise multiplier (5x) to ensure some subjects jump categories,
			# otherwise with SD_NOISE=0.1, matched pairs are almost always concordant.
			return(as.integer(max(1, round(y_t + bt + 5 * eps))))
		}
		stop("Unsupported response_type: ", response_type)
	}

	D = datasets_and_response_models[[dataset_name]]
	X_design = as.data.frame(D$X)
	if (identical(design_type, "KK21stepwise") && ncol(X_design) > 20L){
		message(
			"    Truncating KK21stepwise test covariates from ",
			ncol(X_design),
			" to 20 to keep stepwise-weight tests bounded in runtime."
		)
		X_design = X_design[, seq_len(20L), drop = FALSE]
	}

	if (nrow(X_design) %% 4L != 0L){
		n_keep = nrow(X_design) - (nrow(X_design) %% 4L)
		message(
			"    Truncating test rows from ",
			nrow(X_design),
			" to ",
			n_keep,
			" so the design matrix row count is divisible by 4."
		)
		X_design = X_design[seq_len(n_keep), , drop = FALSE]
	}

	n = nrow(X_design)
	dataset_n_rows = nrow(X_design)
	dataset_n_cols = ncol(X_design)
	dead = rep(1, n)
	y = D$y_original[[response_type]]
	t_f = quantile(y, .95)

	if (response_type == "survival"){
		for (i in 1 : n){
			if (runif(1) < prob_censoring || y[i] >= t_f){
				y[i] = runif(1, 0, y[i])
				dead[i] = 0
			}
		}
	}

	cluster_design_setup = NULL
	if (identical(design_type, "FixedCluster")) {
		cluster_design_setup = add_assignment_only_cluster_id(X_design)
	} else if (identical(design_type, "FixedBlockedCluster")) {
		cluster_design_setup = add_assignment_only_cluster_id(X_design, strata_cols = names(X_design)[2:min(2, ncol(X_design))])
	}
	X_design_for_design = if (is.null(cluster_design_setup)) X_design else cluster_design_setup$X

	# For sequential designs that use stratification, we MUST discretize continuous strata
	# because the DesignSeqOneByOne base class explicitly disallows numeric strata.
	X_design_sequential_strata = X_design
	strata_cols_to_use = names(X_design)[1:min(2, ncol(X_design))]
	if (design_type %in% c("SPBR", "PocockSimon", "RandomBlockSize")) {
		for (col in strata_cols_to_use) {
			if (is.numeric(X_design_sequential_strata[[col]])) {
				med = stats::median(X_design_sequential_strata[[col]], na.rm = TRUE)
				X_design_sequential_strata[[col]] = factor(ifelse(X_design_sequential_strata[[col]] <= med, "low", "high"))
			}
		}
	}

	des_obj = tryCatch(switch(design_type,
		KK21 =         DesignSeqOneByOneKK21$new(        response_type = response_type, n = n),
		KK21stepwise = DesignSeqOneByOneKK21stepwise$new(response_type = response_type, n = n),
		KK14 =         DesignSeqOneByOneKK14$new(        response_type = response_type, n = n),
		Bernoulli =    DesignSeqOneByOneBernoulli$new(   response_type = response_type, n = n),
		Efron =        DesignSeqOneByOneEfron$new(       response_type = response_type, n = n),
		Atkinson =     DesignSeqOneByOneAtkinson$new(    response_type = response_type, n = n),
		iBCRD =        DesignSeqOneByOneiBCRD$new(       response_type = response_type, n = n),
		Urn =          DesignSeqOneByOneUrn$new(         response_type = response_type, n = n),
		RandomBlockSize = DesignSeqOneByOneRandomBlockSize$new( strata_cols = strata_cols_to_use, response_type = response_type, n = n),
		SPBR =         DesignSeqOneByOneSPBR$new(        strata_cols = strata_cols_to_use, block_size = 4, response_type = response_type, n = n),
		PocockSimon =  DesignSeqOneByOnePocockSimon$new( strata_cols = strata_cols_to_use, response_type = response_type, n = n),
		FixedBernoulli = FixedDesignBernoulli$new( response_type = response_type, n = n),
		FixediBCRD =     FixedDesigniBCRD$new(     response_type = response_type, n = n),
		FixedBlocking =  FixedDesignBlocking$new(  strata_cols = strata_cols_to_use, response_type = response_type, n = n),
		FixedCluster =   FixedDesignCluster$new(   cluster_col = cluster_design_setup$cluster_col, response_type = response_type, n = n),
		FixedBlockedCluster = FixedDesignBlockedCluster$new( strata_cols = names(X_design)[2:min(2, ncol(X_design))], cluster_col = cluster_design_setup$cluster_col, response_type = response_type, n = n),
		FixedBinaryMatch = FixedDesignBinaryMatch$new( response_type = response_type, n = n),
		FixedGreedy =    FixedDesignGreedy$new(    response_type = response_type, n = n),
		FixedRerandomization = FixedDesignRerandomization$new( response_type = response_type, n = n),
		FixedMatchingGreedy = FixedDesignMatchingGreedyPairSwitching$new( response_type = response_type, n = n),
		FixedDOptimal =  FixedDesignDOptimal$new(  response_type = response_type, n = n),
		FixedAOptimal =  FixedDesignAOptimal$new(  response_type = response_type, n = n),
		stop("Unsupported design_type: ", design_type)
	), error = function(e){ message("    Skipping design (creation error): ", e$message); NULL })
	if (is.null(des_obj)) return(invisible(NULL))

	if (inherits(des_obj, "DesignSeqOneByOne")){
		seq_ok = tryCatch({
			for (t in 1 : n){
				w_t = des_obj$add_one_subject_to_experiment_and_assign(X_design_sequential_strata[t, , drop = FALSE])
				y_t = apply_treatment_effect_and_noise(y[t], w_t, response_type)
				des_obj$add_one_subject_response(t, y_t, dead[t])
			}
			TRUE
		}, error = function(e){ message("    Skipping design (seq error): ", e$message); FALSE })
		if (!seq_ok) return(invisible(NULL))
	} else {
		# It is a FixedDesign but not a DesignSeqOneByOne
		des_obj$add_all_subjects_to_experiment(X_design_for_design)
		randomize_ok = tryCatch({ des_obj$assign_w_to_all_subjects(); TRUE }, error = function(e){ message("    Skipping design: ", e$message); FALSE })
		if (!randomize_ok) return(invisible(NULL))
		w = des_obj$get_w()
		for (t in 1 : n){
			y_t = apply_treatment_effect_and_noise(y[t], w[t], response_type)
			des_obj$add_one_subject_response(t, y_t, dead[t])
		}
	}

	is_kk_design = design_type %in% c("KK21", "KK21stepwise", "KK14")
	if (response_type == "continuous"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		if (is_kk_design){
			if (design_type == "KK14"){
				inference_banner("InferenceBaiAdjustedTKK14")
				run_inference_checks(InferenceBaiAdjustedTKK14$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			}
			if (design_type %in% c("KK21", "KK21stepwise")){
				inference_banner("InferenceBaiAdjustedTKK21")
				run_inference_checks(InferenceBaiAdjustedTKK21$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			}
			inference_banner("InferenceAllKKCompoundMeanDiff")
			run_inference_checks(InferenceAllKKCompoundMeanDiff$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceAllKKWilcoxIVWC")
			run_inference_checks(InferenceAllKKWilcoxIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrUnivIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrMultiIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultOLSKKCombinedLikelihood")
			run_inference_checks(InferenceContinMultOLSKKCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultOLSKKIVWC")
			run_inference_checks(InferenceContinMultOLSKKIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultiKKLinIVWC")
			run_inference_checks(InferenceContinMultiKKLinIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultiKKLinCombinedLikelihood")
			run_inference_checks(InferenceContinMultiKKLinCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultGLS")
			run_inference_checks(InferenceContinMultGLS$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinUnivKKGLMM")
			run_inference_checks(InferenceContinUnivKKGLMM$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultiKKGLMM")
			run_inference_checks(InferenceContinMultiKKGLMM$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinUnivKKRobustRegrIVWC")
			run_inference_checks(InferenceContinUnivKKRobustRegrIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultiKKRobustRegrIVWC")
			run_inference_checks(InferenceContinMultiKKRobustRegrIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinUnivKKRobustRegrCombinedLikelihood")
			run_inference_checks(InferenceContinUnivKKRobustRegrCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultiKKRobustRegrCombinedLikelihood")
			run_inference_checks(InferenceContinMultiKKRobustRegrCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultKKQuantileRegrIVWC")
			run_inference_checks(InferenceContinMultKKQuantileRegrIVWC$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultKKQuantileRegrCombinedLikelihood")
			run_inference_checks(InferenceContinMultKKQuantileRegrCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		} else {
			inference_banner("InferenceContinUnivRobustRegr")
			run_inference_checks(InferenceContinUnivRobustRegr$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultiRobustRegr")
			run_inference_checks(InferenceContinMultiRobustRegr$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinUnivQuantileRegr")
			run_inference_checks(InferenceContinUnivQuantileRegr$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultiQuantileRegr")
			run_inference_checks(InferenceContinMultiQuantileRegr$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultLin")
			run_inference_checks(InferenceContinMultLin$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("InferenceContinMultOLS")
			run_inference_checks(InferenceContinMultOLS$new(des_obj), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		}
	}

	if (response_type == "incidence"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("InferenceAllKKCompoundMeanDiff")
			run_inference_checks(InferenceAllKKCompoundMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKClogitCombinedLikelihood")
			run_inference_checks(InferenceIncidUnivKKClogitCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiKKClogitCombinedLikelihood")
			run_inference_checks(InferenceIncidMultiKKClogitCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKClogitIVWC")
			run_inference_checks(InferenceIncidUnivKKClogitIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiKKClogitIVWC")
			run_inference_checks(InferenceIncidMultiKKClogitIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKGEE")
			run_inference_checks(InferenceIncidUnivKKGEE$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiKKGEE")
			run_inference_checks(InferenceIncidMultiKKGEE$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKNewcombeRiskDiff")
			run_inference_checks(InferenceIncidUnivKKNewcombeRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKGCompRiskDiff")
			run_inference_checks(InferenceIncidUnivKKGCompRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiKKGCompRiskDiff")
			run_inference_checks(InferenceIncidMultiKKGCompRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKGCompRiskRatio")
			run_inference_checks(InferenceIncidUnivKKGCompRiskRatio$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiKKGCompRiskRatio")
			run_inference_checks(InferenceIncidMultiKKGCompRiskRatio$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKModifiedPoisson")
			run_inference_checks(InferenceIncidUnivKKModifiedPoisson$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiKKModifiedPoisson")
			run_inference_checks(InferenceIncidMultiKKModifiedPoisson$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivKKGLMM")
			run_inference_checks(InferenceIncidUnivKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiKKGLMM")
			run_inference_checks(InferenceIncidMultiKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("InferenceIncidUnivLogRegr")
		run_inference_checks(InferenceIncidUnivLogRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceIncidMultiLogRegr")
		run_inference_checks(InferenceIncidMultiLogRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (!is_kk_design){
			inference_banner("InferenceIncidUnivMiettinenNurminenRiskDiff")
			run_inference_checks(InferenceIncidUnivMiettinenNurminenRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivNewcombeRiskDiff")
			run_inference_checks(InferenceIncidUnivNewcombeRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivRiskDiff")
			run_inference_checks(InferenceIncidUnivRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiRiskDiff")
			run_inference_checks(InferenceIncidMultiRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivGCompRiskDiff")
			run_inference_checks(InferenceIncidUnivGCompRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiGCompRiskDiff")
			run_inference_checks(InferenceIncidMultiGCompRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivGCompRiskRatio")
			run_inference_checks(InferenceIncidUnivGCompRiskRatio$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiGCompRiskRatio")
			run_inference_checks(InferenceIncidMultiGCompRiskRatio$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivModifiedPoisson")
			run_inference_checks(InferenceIncidUnivModifiedPoisson$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiModifiedPoisson")
			run_inference_checks(InferenceIncidMultiModifiedPoisson$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivLogBinomial")
			run_inference_checks(InferenceIncidUnivLogBinomial$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiLogBinomial")
			run_inference_checks(InferenceIncidMultiLogBinomial$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidUnivBinomialIdentityRiskDiff")
			run_inference_checks(InferenceIncidUnivBinomialIdentityRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceIncidMultiBinomialIdentityRiskDiff")
			run_inference_checks(InferenceIncidMultiBinomialIdentityRiskDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
	}

	if (response_type == "proportion"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("InferenceAllKKCompoundMeanDiff")
			run_inference_checks(InferenceAllKKCompoundMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxIVWC")
			run_inference_checks(InferenceAllKKWilcoxIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrUnivIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrMultiIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropUnivKKGEE")
			run_inference_checks(InferencePropUnivKKGEE$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropMultiKKGEE")
			run_inference_checks(InferencePropMultiKKGEE$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropUnivKKGLMM")
			run_inference_checks(InferencePropUnivKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropMultiKKGLMM")
			run_inference_checks(InferencePropMultiKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropMultiKKQuantileRegrIVWC")
			run_inference_checks(InferencePropMultiKKQuantileRegrIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropMultiKKQuantileRegrCombinedLikelihood")
			run_inference_checks(InferencePropMultiKKQuantileRegrCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("InferencePropUniBetaRegr")
		run_inference_checks(InferencePropUniBetaRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferencePropMultiBetaRegr")
		run_inference_checks(InferencePropMultiBetaRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (!is_kk_design){
			inference_banner("InferencePropUniZeroOneInflatedBetaRegr")
			run_inference_checks(InferencePropUniZeroOneInflatedBetaRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropMultiZeroOneInflatedBetaRegr")
			run_inference_checks(InferencePropMultiZeroOneInflatedBetaRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropUniGCompMeanDiff")
			run_inference_checks(InferencePropUniGCompMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropMultiGCompMeanDiff")
			run_inference_checks(InferencePropMultiGCompMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropUniFractionalLogit")
			run_inference_checks(InferencePropUniFractionalLogit$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferencePropMultiFractionalLogit")
			run_inference_checks(InferencePropMultiFractionalLogit$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
	}

	if (response_type == "count"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("InferenceAllKKCompoundMeanDiff")
			run_inference_checks(InferenceAllKKCompoundMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxIVWC")
			run_inference_checks(InferenceAllKKWilcoxIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrUnivIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrMultiIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonUnivKKGEE")
			run_inference_checks(InferenceCountPoissonUnivKKGEE$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonMultiKKGEE")
			run_inference_checks(InferenceCountPoissonMultiKKGEE$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonUnivKKCPoissonCombinedLikelihood")
			run_inference_checks(InferenceCountPoissonUnivKKCPoissonCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonMultiKKCPoissonCombinedLikelihood")
			run_inference_checks(InferenceCountPoissonMultiKKCPoissonCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivKKHurdlePoissonCombinedLikelihood")
			run_inference_checks(InferenceCountUnivKKHurdlePoissonCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiKKHurdlePoissonCombinedLikelihood")
			run_inference_checks(InferenceCountMultiKKHurdlePoissonCombinedLikelihood$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivKKHurdlePoissonIVWC")
			run_inference_checks(InferenceCountUnivKKHurdlePoissonIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiKKHurdlePoissonIVWC")
			run_inference_checks(InferenceCountMultiKKHurdlePoissonIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonUnivKKCPoissonIVWC")
			run_inference_checks(InferenceCountPoissonUnivKKCPoissonIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonMultiKKCPoissonIVWC")
			run_inference_checks(InferenceCountPoissonMultiKKCPoissonIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonUnivKKGLMM")
			run_inference_checks(InferenceCountPoissonUnivKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountPoissonMultiKKGLMM")
			run_inference_checks(InferenceCountPoissonMultiKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		if (!is_kk_design){
			inference_banner("InferenceCountUnivPoissonRegr")
			run_inference_checks(InferenceCountUnivPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiPoissonRegr")
			run_inference_checks(InferenceCountMultiPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivRobustPoissonRegr")
			run_inference_checks(InferenceCountUnivRobustPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiRobustPoissonRegr")
			run_inference_checks(InferenceCountMultiRobustPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivQuasiPoissonRegr")
			run_inference_checks(InferenceCountUnivQuasiPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiQuasiPoissonRegr")
			run_inference_checks(InferenceCountMultiQuasiPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivZeroInflatedPoissonRegr")
			run_inference_checks(InferenceCountUnivZeroInflatedPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiZeroInflatedPoissonRegr")
			run_inference_checks(InferenceCountMultiZeroInflatedPoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivZeroInflatedNegBinRegr")
			run_inference_checks(InferenceCountUnivZeroInflatedNegBinRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiZeroInflatedNegBinRegr")
			run_inference_checks(InferenceCountMultiZeroInflatedNegBinRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivHurdlePoissonRegr")
			run_inference_checks(InferenceCountUnivHurdlePoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiHurdlePoissonRegr")
			run_inference_checks(InferenceCountMultiHurdlePoissonRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountUnivHurdleNegBinRegr")
			run_inference_checks(InferenceCountUnivHurdleNegBinRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceCountMultiHurdleNegBinRegr")
			run_inference_checks(InferenceCountMultiHurdleNegBinRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("InferenceCountUnivNegBinRegr")
		run_inference_checks(InferenceCountUnivNegBinRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceCountMultiNegBinRegr")
		run_inference_checks(InferenceCountMultiNegBinRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}

	if (response_type == "survival"){
		if (is_kk_design){
			is_censoring_skip_error = function(msg){
				grepl("only available for uncensored", msg, fixed = TRUE) ||
				grepl("does not currently support censored", msg, fixed = TRUE) ||
				grepl("does not support censored", msg, fixed = TRUE)
			}
			for (kk_surv_class in list(InferenceAllKKCompoundMeanDiff, InferenceAllKKWilcoxIVWC, InferenceAllKKWilcoxRegrUnivIVWC, InferenceAllKKWilcoxRegrMultiIVWC)){
				class_name = kk_surv_class$classname
				inference_banner(class_name)
				err_msg = tryCatch({
					run_inference_checks(kk_surv_class$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
					NULL
				}, error = function(e) if (length(e$message) == 0L) "" else e$message)
				if (!is.null(err_msg)){
					if (is_censoring_skip_error(err_msg)) message("  Skipping ", class_name, " (censored data): ", err_msg)
					else stop(err_msg)
				}
			}
			for (kk_surv_class in list(
				InferenceSurvivalUnivKKGEE,
				InferenceSurvivalMultiKKGammaGEE,
				InferenceSurvivalUnivKKGLMM,
				InferenceSurvivalMultiKKGammaGLMM,
				InferenceSurvivalUnivKKClaytonCopulaIVWC,
				InferenceSurvivalMultiKKClaytonCopulaIVWC,
				InferenceSurvivalUnivKKLWACoxIVWC,
				InferenceSurvivalMultiKKLWACoxIVWC,
				InferenceSurvivalUnivKKStratCoxIVWC,
				InferenceSurvivalMultiKKStratCoxIVWC,
				InferenceSurvivalUnivKKRankRegrIVWC,
				InferenceSurvivalMultiKKRankRegrIVWC,
				InferenceSurvivalUnivKKWeibullFrailtyIVWC,
				InferenceSurvivalMultiKKWeibullFrailtyIVWC,
				InferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood,
				InferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood,
				InferenceSurvivalUnivKKLWACoxCombinedLikelihood,
				InferenceSurvivalMultiKKLWACoxCombinedLikelihood,
				InferenceSurvivalUnivKKStratCoxCombinedLikelihood,
				InferenceSurvivalMultiKKStratCoxCombinedLikelihood,
				InferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood,
				InferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood
			)){
				class_name = kk_surv_class$classname
				inference_banner(class_name)
				err_msg = tryCatch({
					run_inference_checks(kk_surv_class$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
					NULL
				}, error = function(e) if (length(e$message) == 0L) "" else e$message)
				if (!is.null(err_msg)){
					if (grepl("only available for uncensored", err_msg, fixed = TRUE)) message("  Skipping ", class_name, " (censored data): ", err_msg)
					else stop(err_msg)
				}
			}
		}
		inference_banner("InferenceAllSimpleWilcox")
		err_msg_sw = tryCatch({
			run_inference_checks(InferenceAllSimpleWilcox$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			NULL
		}, error = function(e) if (length(e$message) == 0L) "" else e$message)
		if (!is.null(err_msg_sw)){
			if (grepl("does not support censored", err_msg_sw, fixed = TRUE)) message("  Skipping InferenceAllSimpleWilcox (censored data): ", err_msg_sw)
			else stop(err_msg_sw)
		}
		inference_banner("InferenceSurvivalGehanWilcox")
		run_inference_checks(InferenceSurvivalGehanWilcox$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalLogRank")
		run_inference_checks(InferenceSurvivalLogRank$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalRestrictedMeanDiff")
		run_inference_checks(InferenceSurvivalRestrictedMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalKMDiff")
		run_inference_checks(InferenceSurvivalKMDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalUniWeibullRegr")
		run_inference_checks(InferenceSurvivalUniWeibullRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalMultiWeibullRegr")
		run_inference_checks(InferenceSurvivalMultiWeibullRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalUniDepCensTransformRegr")
		run_inference_checks(InferenceSurvivalUniDepCensTransformRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalMultiDepCensTransformRegr")
		run_inference_checks(InferenceSurvivalMultiDepCensTransformRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalUniCoxPHRegr")
		run_inference_checks(InferenceSurvivalUniCoxPHRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalMultiCoxPHRegr")
		run_inference_checks(InferenceSurvivalMultiCoxPHRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalUniStratCoxPHRegr")
		run_inference_checks(InferenceSurvivalUniStratCoxPHRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceSurvivalMultiStratCoxPHRegr")
		run_inference_checks(InferenceSurvivalMultiStratCoxPHRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}

	if (response_type == "ordinal"){
		inference_banner("InferenceAllSimpleMeanDiff")
		run_inference_checks(InferenceAllSimpleMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceAllSimpleWilcox")
		run_inference_checks(InferenceAllSimpleWilcox$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("InferenceAllKKCompoundMeanDiff")
			run_inference_checks(InferenceAllKKCompoundMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxIVWC")
			run_inference_checks(InferenceAllKKWilcoxIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrUnivIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(InferenceAllKKWilcoxRegrMultiIVWC$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalUnivKKGEE")
			run_inference_checks(InferenceOrdinalUnivKKGEE$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalUnivKKGLMM")
			run_inference_checks(InferenceOrdinalUnivKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalMultiKKGLMM")
			run_inference_checks(InferenceOrdinalMultiKKGLMM$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalUnivKKGLMMProbit")
			run_inference_checks(InferenceOrdinalUnivKKGLMMProbit$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalMultiKKGLMMProbit")
			run_inference_checks(InferenceOrdinalMultiKKGLMMProbit$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalUnivKKCondPropOddsRegr")
			run_inference_checks(InferenceOrdinalUnivKKCondPropOddsRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalUnivKKCondPropOddsCombinedRegr")
			run_inference_checks(InferenceOrdinalUnivKKCondPropOddsCombinedRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalUnivKKCondContRatioRegr")
			run_inference_checks(InferenceOrdinalUnivKKCondContRatioRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalMultiKKCondContRatioRegr")
			run_inference_checks(InferenceOrdinalMultiKKCondContRatioRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalUnivKKCondAdjCatLogitRegr")
			run_inference_checks(InferenceOrdinalUnivKKCondAdjCatLogitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalMultiKKCondAdjCatLogitRegr")
			run_inference_checks(InferenceOrdinalMultiKKCondAdjCatLogitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("InferenceOrdinalPairedSignTest")
			run_inference_checks(InferenceOrdinalPairedSignTest$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("InferenceOrdinalUniAdjCatLogitRegr")
		run_inference_checks(InferenceOrdinalUniAdjCatLogitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiAdjCatLogitRegr")
		run_inference_checks(InferenceOrdinalMultiAdjCatLogitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniCumulProbitRegr")
		run_inference_checks(InferenceOrdinalUniCumulProbitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiCumulProbitRegr")
		run_inference_checks(InferenceOrdinalMultiCumulProbitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniOrderedProbitRegr")
		run_inference_checks(InferenceOrdinalUniOrderedProbitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiOrderedProbitRegr")
		run_inference_checks(InferenceOrdinalMultiOrderedProbitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniStereotypeLogitRegr")
		run_inference_checks(InferenceOrdinalUniStereotypeLogitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiStereotypeLogitRegr")
		run_inference_checks(InferenceOrdinalMultiStereotypeLogitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniStereotypeProbitRegr")
		run_inference_checks(InferenceOrdinalUniStereotypeProbitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiStereotypeProbitRegr")
		run_inference_checks(InferenceOrdinalMultiStereotypeProbitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniPropOddsRegr")
		run_inference_checks(InferenceOrdinalUniPropOddsRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniPartialProportionalOddsRegr")
		run_inference_checks(InferenceOrdinalUniPartialProportionalOddsRegr$new(des_obj, verbose = FALSE), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiPartialProportionalOddsRegr")
		run_inference_checks(InferenceOrdinalMultiPartialProportionalOddsRegr$new(des_obj, verbose = FALSE), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniCLLRegr")
		run_inference_checks(InferenceOrdinalUniCLLRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiCLLRegr")
		run_inference_checks(InferenceOrdinalMultiCLLRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniGCompMeanDiff")
		run_inference_checks(InferenceOrdinalUniGCompMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiGCompMeanDiff")
		run_inference_checks(InferenceOrdinalMultiGCompMeanDiff$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalJonckheereTerpstraTest")
		run_inference_checks(InferenceOrdinalJonckheereTerpstraTest$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalUniCauchitRegr")
		run_inference_checks(InferenceOrdinalUniCauchitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiCauchitRegr")
		run_inference_checks(InferenceOrdinalMultiCauchitRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalContRatioRegr")
		run_inference_checks(InferenceOrdinalContRatioRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalMultiContRatioRegr")
		run_inference_checks(InferenceOrdinalMultiContRatioRegr$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("InferenceOrdinalRidit")
		run_inference_checks(InferenceOrdinalRidit$new(des_obj), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}
}


for (rep_curr in 1:Nrep) {
	pending_rep_header <<- paste0("\n\n====== REPETITION ", rep_curr, " OF ", Nrep, " ======\n")
	for (beta_T_iter_curr in seq_along(beta_T_values)){
		beta_T = beta_T_values[beta_T_iter_curr]
		pending_beta_header <<- paste0("  === beta_T = [", beta_T, "] ===")
		for (dataset_name in names(datasets_and_response_models)){
			pending_dataset_header <<- paste0("    === dataset: ", dataset_name, " ===")
			for (response_type in ALL_RESPONSE_TYPES) {
				if (!is.na(RESPONSE_TYPE_FILTER) && !identical(response_type, RESPONSE_TYPE_FILTER)) {
					next
				}
				if (!(response_type %in% names(datasets_and_response_models[[dataset_name]]$y_original))) {
					next
				}
				pending_response_header <<- paste0("      === response_type: ", response_type, " ===")
				for (design_type in c("Bernoulli", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise", "SPBR", "PocockSimon", "Urn", "RandomBlockSize", "FixedBernoulli", "FixediBCRD", "FixedBlocking", "FixedCluster", "FixedBlockedCluster", "FixedBinaryMatch", "FixedGreedy", "FixedRerandomization", "FixedMatchingGreedy", "FixedDOptimal", "FixedAOptimal")) {
					pending_design_header <<- paste0("        === design: ", design_type, " ===")
					run_tests_for_response(response_type, design_type = design_type, dataset_name = dataset_name)
				}
			}
		}
	}
}
message("\n\n----------------------All tests complete!")
write_results_if_needed(force = TRUE)

unset_num_cores()
rm(list=ls())
