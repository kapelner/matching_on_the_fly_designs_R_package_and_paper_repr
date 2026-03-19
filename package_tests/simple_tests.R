#/c/Program\ Files/R/R-devel/bin/R.exe CMD INSTALL -l ~/AppData/Local/R/win-library/4.5/ SeqExpMatch/
rm(list = ls())
set.seed(1)
pacman::p_load(SeqExpMatch, doParallel, PTE, datasets, qgam, mlbench, AppliedPredictiveModeling, dplyr, ggplot2, gridExtra, profvis, data.table, profvis)
max_n_dataset = 150
source("package_tests/_dataset_load.R")
# options(error = recover)
# options(warn=2)

args = commandArgs(trailingOnly = TRUE)
Nrep      = if (length(args) >= 1) as.integer(args[1]) else 40L
NUM_CORES = if (length(args) >= 2) as.integer(args[2]) else 2L
prob_censoring = 0.15
nsim_exact_test = 351
pval_epsilon = 0.007
test_compute_confidence_interval_rand = TRUE
beta_T_values = c(0, 1)
SD_NOISE = 0.1

results_file = paste0("package_tests/simple_tests_results_nc_", NUM_CORES, ".csv")
existing_results_dt = if (file.exists(results_file)) data.table::fread(results_file) else data.table::data.table()
run_row_id = if ("run_row_id" %in% colnames(existing_results_dt) && nrow(existing_results_dt) > 0L) {
	as.integer(max(existing_results_dt$run_row_id, na.rm = TRUE))
} else {
	0L
}

serialize_beta_T = function(value){
	if (is.na(value)) return("NA")
	if (is.numeric(value)) return(format(value, scientific = TRUE, digits = 17))
	as.character(value)
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

is_row_completed = function(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val){
	key = build_result_key(rep_val, beta_val, dataset_val, response_val, design_val, inference_val, function_run_val)
	!is.null(completed_rows_cache[[key]])
}

if (nrow(existing_results_dt) > 0L) {
	for (row_idx in seq_len(nrow(existing_results_dt))) {
		row = existing_results_dt[row_idx]
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
	duration_time_sec = numeric(),
	inference_class = character(),
	dataset = character(),
	design = character(),
	response_type = character(),
	function_run = character(),
	result_1 = character(),
	result_2 = character(),
	beta_T_in_confidence_interval = logical(),
	error_message = character(),
	rep = integer(),
	run_row_id = integer(),
	beta_T = numeric(),
	nsim_exact_test = integer(),
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
			col.names = !append_mode
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
	result_1 = if (length(result_vec) >= 1) result_vec[1] else NA_character_
	result_2 = if (grepl("confidence_interval", function_run, fixed = TRUE) && length(result_vec) >= 2) result_vec[2] else NA_character_
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
			run_row_id = run_row_id,
			duration_time_sec = duration_time_sec,
			beta_T = beta_T,
			nsim_exact_test = as.integer(nsim_exact_test),
			pval_epsilon = pval_epsilon,
			prob_censoring = prob_censoring,
			sd_noise = SD_NOISE,
			num_cores = as.integer(NUM_CORES),
			inference_class = inference_class,
			dataset = dataset_name,
			dataset_n_rows = as.integer(dataset_n_rows),
			dataset_n_cols = as.integer(dataset_n_cols),
			design = design_type,
			response_type = response_type,
			function_run = function_run,
			result_1 = result_1,
			result_2 = result_2,
			beta_T_in_confidence_interval = beta_T_in_confidence_interval,
			result = result_str,
			status = status,
			error_message = ifelse(is.null(error_message), NA_character_, as.character(error_message))
		)
	), use.names = TRUE)
	mark_row_completed(rep_curr, beta_T, dataset_name, response_type, design_type, inference_class, function_run)
	write_results_if_needed(force = TRUE)
}

run_inference_checks = function(seq_des_inf, response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols){
	skip_slow = is(seq_des_inf, "SeqDesignInferencePropMultiBetaRegr") || is(seq_des_inf, "SeqDesignInferencePropUniZeroOneInflatedBetaRegr") || is(seq_des_inf, "SeqDesignInferencePropMultiZeroOneInflatedBetaRegr") || is(seq_des_inf, "SeqDesignInferenceSurvivalUniDepCensTransformRegr") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiDepCensTransformRegr") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiWeibullRegr") || is(seq_des_inf, "SeqDesignInferenceCountMultiNegBinRegr") || is(seq_des_inf, "SeqDesignInferenceCountUnivZeroInflatedPoissonRegr") || is(seq_des_inf, "SeqDesignInferenceCountMultiZeroInflatedPoissonRegr") || is(seq_des_inf, "SeqDesignInferenceCountUnivZeroInflatedNegBinRegr") || is(seq_des_inf, "SeqDesignInferenceCountMultiZeroInflatedNegBinRegr") || is(seq_des_inf, "SeqDesignInferenceCountUnivHurdlePoissonRegr") || is(seq_des_inf, "SeqDesignInferenceCountMultiHurdlePoissonRegr") || is(seq_des_inf, "SeqDesignInferenceCountUnivHurdleNegBinRegr") || is(seq_des_inf, "SeqDesignInferenceCountMultiHurdleNegBinRegr") || is(seq_des_inf, "SeqDesignInferenceCountUnivKKHurdlePoissonCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceCountMultiKKHurdlePoissonCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceCountUnivKKHurdlePoissonIVWC") || is(seq_des_inf, "SeqDesignInferenceCountMultiKKHurdlePoissonIVWC") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiCoxPHRegr") || is(seq_des_inf, "SeqDesignInferenceSurvivalUniCoxPHRegr") || is(seq_des_inf, "SeqDesignInferenceSurvivalUnivKKClaytonCopulaIVWC") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKClaytonCopulaIVWC") || is(seq_des_inf, "SeqDesignInferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKLWACoxIVWC") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKStratCoxIVWC") || is(seq_des_inf, "SeqDesignInferenceCountPoissonMultiKKCPoissonIVWC") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKWeibullFrailtyIVWC") || is(seq_des_inf, "SeqDesignInferenceAllKKWilcoxRegrMultiIVWC") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKRankRegrIVWC") || is(seq_des_inf, "SeqDesignInferenceCountPoissonMultiKKCPoissonCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKLWACoxCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKStratCoxCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceContinMultLin") || is(seq_des_inf, "SeqDesignInferenceContinUnivQuantileRegr") || is(seq_des_inf, "SeqDesignInferenceContinMultiQuantileRegr") || is(seq_des_inf, "SeqDesignInferenceContinMultiKKLinIVWC") || is(seq_des_inf, "SeqDesignInferenceContinMultiKKLinCombinedLikelihood")
	skip_bootstrap = is(seq_des_inf, "SeqDesignInferenceAbstractKKGEE") || is(seq_des_inf, "SeqDesignInferenceAbstractKKGLMM") || is(seq_des_inf, "SeqDesignInferenceContinMultGLS") || is(seq_des_inf, "SeqDesignInferenceAbstractKKClaytonCopulaIVWC") || is(seq_des_inf, "SeqDesignInferenceAbstractKKClaytonCopulaCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceAbstractKKWeibullFrailtyIVWC") || is(seq_des_inf, "SeqDesignInferenceAbstractKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceAllKKWilcoxIVWC") || is(seq_des_inf, "SeqDesignInferenceAbstractKKWilcoxRegrIVWC") || is(seq_des_inf, "SeqDesignInferenceSurvivalUnivKKRankRegrIVWC") || is(seq_des_inf, "SeqDesignInferenceIncidExactZhangAbstract") || is(seq_des_inf, "SeqDesignInferenceAllSimpleWilcox") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKGEE") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKGLMM") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKGLMM") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKGLMMProbit") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKGLMMProbit") || is(seq_des_inf, "SeqDesignInferenceOrdinalPairedSignTest") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKCondContRatioRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKCondAdjCatLogitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUniGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiCLLRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUniOrderedProbitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiOrderedProbitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUniCauchitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiCauchitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiContRatioRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKGEE") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKCondContRatioRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKCondAdjCatLogitRegr")
	skip_rand      = is(seq_des_inf, "SeqDesignInferenceAbstractKKGEE") || is(seq_des_inf, "SeqDesignInferenceAbstractKKGLMM") || is(seq_des_inf, "SeqDesignInferenceIncidExactZhangAbstract") || is(seq_des_inf, "SeqDesignInferencePropUniGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferencePropMultiGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKGEE") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKGLMM") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKGLMM") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKGLMMProbit") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKGLMMProbit") || is(seq_des_inf, "SeqDesignInferenceOrdinalPairedSignTest") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKCondContRatioRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUnivKKCondAdjCatLogitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUniGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiCLLRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUniOrderedProbitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiOrderedProbitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalUniCauchitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiCauchitRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiContRatioRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKGEE") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKCondContRatioRegr") || is(seq_des_inf, "SeqDesignInferenceOrdinalMultiKKCondAdjCatLogitRegr")
	skip_mle_pval  = is(seq_des_inf, "SeqDesignInferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood")
	skip_rand_pval = is(seq_des_inf, "SeqDesignInferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferenceContinMultGLS") || is(seq_des_inf, "SeqDesignInferencePropUniGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferencePropMultiGCompMeanDiff")
	skip_ci_rand   = is(seq_des_inf, "SeqDesignInferenceContinMultKKQuantileRegrIVWC") || is(seq_des_inf, "SeqDesignInferencePropMultiKKQuantileRegrIVWC") || is(seq_des_inf, "SeqDesignInferenceContinMultKKQuantileRegrCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferencePropMultiKKQuantileRegrCombinedLikelihood") || is(seq_des_inf, "SeqDesignInferencePropUniGCompMeanDiff") || is(seq_des_inf, "SeqDesignInferencePropMultiGCompMeanDiff")
	skip_ci_rand_custom = is(seq_des_inf, "SeqDesignInferenceContinMultiKKGLMM") || is(seq_des_inf, "SeqDesignInferenceContinUnivKKGLMM")
	
	skip_ci = beta_T == 1 && (
		is(seq_des_inf, "SeqDesignInferenceIncidMultiLogRegr") ||
		is(seq_des_inf, "SeqDesignInferencePropUniBetaRegr") ||
		is(seq_des_inf, "SeqDesignInferencePropMultiBetaRegr") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalUniCoxPHRegr") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiCoxPHRegr") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKLWACoxIVWC") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKStratCoxIVWC") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKClaytonCopulaIVWC") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKLWACoxCombinedLikelihood") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKStratCoxCombinedLikelihood") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKWeibullFrailtyIVWC") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood") ||
		is(seq_des_inf, "SeqDesignInferenceAllKKWilcoxRegrMultiIVWC") ||
		is(seq_des_inf, "SeqDesignInferenceSurvivalMultiKKRankRegrIVWC")
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
		message("Skipping ", class(seq_des_inf)[1], " / ", label, " (already recorded for rep ", rep_curr, ")")
		return(invisible(NULL))
	}
	start_elapsed = unname(proc.time()[["elapsed"]])
	tryCatch({
		result <- expr
			if (has_invalid_numeric(result)) stop("Invalid output detected (NA/NaN/Inf) in ", label)
			if (is_zero_zero_confidence_interval(label, result)) stop("Degenerate confidence interval [0, 0] detected in ", label)
			result = snap_small_numeric_to_zero(result)
			cat("  ", result, "\n")
			duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
			record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, result, status = "ok", duration_time_sec = duration_time_sec)
			result
		}, error = function(e){
			msg = if (length(e$message) == 0L) "" else e$message
			is_non_fatal = grepl("not implemented", msg, fixed = TRUE) ||
						 grepl("not enough discordant pairs", msg, ignore.case = TRUE) ||
						 grepl("Degenerate confidence interval", msg, fixed = TRUE) ||
						 grepl("inconsistent estimator units", msg, ignore.case = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned NA bounds", msg, fixed = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned non-finite bounds", msg, fixed = TRUE) ||
						 					 						 grepl("Weibull regression failed to converge", msg, fixed = TRUE) ||
						 					 						 grepl("Invalid output detected", msg, fixed = TRUE) ||
						 					 						 grepl("missing value where TRUE/FALSE needed", msg, fixed = TRUE) ||
						 					 						 ((grepl("NA/NaN/Inf", msg, fixed = TRUE) || grepl("non-finite standard error", msg, fixed = TRUE) || grepl("could not compute a finite standard error", msg, fixed = TRUE)) &&
						 					 
						 													 (is(seq_des_inf, "SeqDesignInferenceAbstractKKClogitIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKPoissonCPoissonIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKStratCoxIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKWeibullFrailtyIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKWilcoxRegrIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKSurvivalRankRegrIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKGEE") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKGLMM") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKRobustRegrIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKRobustRegrCombinedLikelihood") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKQuantileRegrIVWC") ||
						 													  is(seq_des_inf, "SeqDesignInferenceAbstractKKQuantileRegrCombinedLikelihood") ||
						 													  													  is(seq_des_inf, "SeqDesignInferencePropUniFractionalLogit") ||
						 													  													  is(seq_des_inf, "SeqDesignInferenceIncidUnivRiskDiff") ||
						 													  													  is(seq_des_inf, "SeqDesignInferenceIncidMultiRiskDiff") ||
						 													  													  is(seq_des_inf, "SeqDesignInferenceIncidUnivGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "SeqDesignInferenceIncidMultiGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "SeqDesignInferenceIncidUnivGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "SeqDesignInferenceIncidMultiGCompRiskRatio") ||
						 													  is(seq_des_inf, "SeqDesignInferenceIncidUnivModifiedPoisson") ||
						 													  is(seq_des_inf, "SeqDesignInferenceIncidMultiModifiedPoisson") ||
						 													  is(seq_des_inf, "SeqDesignInferencePropUniGCompMeanDiff") ||
						 													  is(seq_des_inf, "SeqDesignInferencePropMultiGCompMeanDiff") ||
						 													  is(seq_des_inf, "SeqDesignInferencePropUniZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "SeqDesignInferencePropMultiZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "SeqDesignInferenceContinMultGLS")))
						 													  										if (isTRUE(is_non_fatal)){
				message("Skipping ", label, " (non-fatal): ", e$message)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = e$message)
			} else {
				stop(e$message)
			}
		})
	}

	if (is(seq_des_inf, "SeqDesignInferenceIncidExactZhangAbstract")){
		message("    Calling compute_exact_two_sided_pval_for_treatment_effect()")
		safe_call("compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		message("    Calling compute_exact_confidence_interval()")
		safe_call("compute_exact_confidence_interval", seq_des_inf$compute_exact_confidence_interval(pval_epsilon = pval_epsilon))
		return(invisible(NULL))
	}

	if (is(seq_des_inf, "SeqDesignInferenceOrdinalJonckheereTerpstraTest")){
		message("    Calling compute_exact_two_sided_pval_for_treatment_effect()")
		safe_call("compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		message("    Calling compute_treatment_estimate()")
		safe_call("compute_treatment_estimate", seq_des_inf$compute_treatment_estimate())
		return(invisible(NULL))
	}

	message("    Calling compute_exact_two_sided_pval_for_treatment_effect()")
	safe_call("compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
	message("    Calling compute_exact_confidence_interval()")
	safe_call("compute_exact_confidence_interval", seq_des_inf$compute_exact_confidence_interval(pval_epsilon = pval_epsilon))

	message("    Calling compute_treatment_estimate()")
	safe_call("compute_treatment_estimate", seq_des_inf$compute_treatment_estimate())
	if (!skip_mle_pval){
		message("    Calling compute_asymp_two_sided_pval_for_treatment_effect()")
		safe_call("compute_asymp_two_sided_pval_for_treatment_effect", seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect())
	}
	if (!skip_ci){
		message("    Calling compute_asymp_confidence_interval()")
		safe_call("compute_asymp_confidence_interval", seq_des_inf$compute_asymp_confidence_interval(0.05))
	}
	if (!skip_slow && !skip_ci && !skip_bootstrap){
		message("    Calling compute_bootstrap_confidence_interval()")
		safe_call("compute_bootstrap_confidence_interval", seq_des_inf$compute_bootstrap_confidence_interval(B = nsim_exact_test, na.rm = TRUE))
	}
	if (!skip_slow && !skip_bootstrap){
		message("    Calling compute_bootstrap_two_sided_pval()")
		safe_call("compute_bootstrap_two_sided_pval", seq_des_inf$compute_bootstrap_two_sided_pval(B = nsim_exact_test, na.rm = TRUE))
	}
	if (!skip_slow && !skip_rand && !skip_rand_pval && response_type %in% c("continuous", "survival", "proportion")){
		message("    Calling compute_two_sided_pval_for_treatment_effect_rand()")
		safe_call("compute_two_sided_pval_for_treatment_effect_rand", seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = nsim_exact_test, show_progress = FALSE))
		message("    Calling compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)")
		use_log   = response_type == "survival"
		use_logit = response_type == "proportion"
		transform_for_rand = if (use_log) "log" else if (use_logit) "logit" else "none"
		delta_for_rand = 0.5
		safe_call("compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)",
				seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = nsim_exact_test, delta = delta_for_rand, transform_responses = transform_for_rand, show_progress = FALSE))
	}

	if (!skip_slow && !skip_rand && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous", "proportion", "count")){
		message("    Calling compute_confidence_interval_rand()")
		safe_call("compute_confidence_interval_rand", seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, show_progress = FALSE))
	}
	if (response_type != "incidence"){
		message("    Calling set_custom_randomization_statistic_function()")
		seq_des_inf$set_custom_randomization_statistic_function(function(){
			yTs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 1]
			yCs = private$seq_des_obj_priv_int$y[private$seq_des_obj_priv_int$w == 0]
			(mean(yTs) - mean(yCs)) / sqrt(var(yTs) / length(yTs) + var(yCs) / length(yCs))
		})
		if (!skip_slow && !skip_rand_pval){
			message("    Calling compute_two_sided_pval_for_treatment_effect_rand(custom)")
			safe_call("compute_two_sided_pval_for_treatment_effect_rand(custom)", seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(nsim_exact_test = nsim_exact_test, show_progress = FALSE))
		}
		if (!skip_slow && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous")){
			message("    Calling compute_confidence_interval_rand(custom)")
			if (!skip_ci_rand_custom){
				safe_call("compute_confidence_interval_rand(custom)", seq_des_inf$compute_confidence_interval_rand(nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, show_progress = FALSE))
			} else {
				message("    Skipping compute_confidence_interval_rand(custom) (too slow)")
			}
		}
		seq_des_inf$set_custom_randomization_statistic_function(NULL)
	}
}

run_tests_for_response = function(response_type, design_type, dataset_name){
	inference_banner = function(inf_name){
		message(paste0("\n\n  == Inference: ", inf_name, " design_type = ", design_type, " dataset = ", dataset_name, " response_type = ", response_type, " beta_T = [", beta_T, "] num_cores = [", NUM_CORES, "] rep = [", rep_curr, "/", Nrep, "]\n"))
	}

	apply_treatment_effect_and_noise = function(y_t, w_t, response_type){
		eps = rnorm(1, 0, SD_NOISE)
		bt = ifelse(w_t == 1, beta_T, 0)
		if (response_type == "continuous") return(y_t + bt + eps)
		if (response_type == "incidence"){
			p_base = ifelse(y_t > 0, 0.75, 0.25)
			p_t = plogis(qlogis(p_base) + bt + eps)
			return(as.numeric(stats::rbinom(1, size = 1, prob = p_t)))
		}
		if (response_type == "proportion"){
			y_clamped = pmax(.Machine$double.eps, pmin(1 - .Machine$double.eps, y_t))
			return(plogis(qlogis(y_clamped) + bt + eps))
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
			return(as.integer(max(1, round(y_t + bt + eps))))
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

	seq_des_obj = switch(design_type,
		KK21 =         SeqDesignKK21$new(        response_type = response_type, n = n),
		KK21stepwise = SeqDesignKK21stepwise$new(response_type = response_type, n = n),
		KK14 =         SeqDesignKK14$new(        response_type = response_type, n = n),
		CRD =          SeqDesignCRD$new(         response_type = response_type, n = n),
		Efron =        SeqDesignEfron$new(       response_type = response_type, n = n),
		Atkinson =     SeqDesignAtkinson$new(    response_type = response_type, n = n),
		iBCRD =        SeqDesigniBCRD$new(       response_type = response_type, n = n),
		stop("Unsupported design_type: ", design_type)
	)
	for (t in 1 : n){
		w_t = seq_des_obj$add_subject_to_experiment_and_assign(X_design[t, , drop = FALSE])
		y_t = apply_treatment_effect_and_noise(y[t], w_t, response_type)
		seq_des_obj$add_subject_response(t, y_t, dead[t])
	}

	is_kk_design = design_type %in% c("KK21", "KK21stepwise", "KK14")
	if (response_type == "continuous"){
		inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
		run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		inference_banner("SeqDesignInferenceAllSimpleWilcox")
		run_inference_checks(SeqDesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		if (is_kk_design){
			if (design_type == "KK14"){
				inference_banner("SeqDesignInferenceBaiAdjustedTKK14")
				run_inference_checks(SeqDesignInferenceBaiAdjustedTKK14$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			}
			if (design_type %in% c("KK21", "KK21stepwise")){
				inference_banner("SeqDesignInferenceBaiAdjustedTKK21")
				run_inference_checks(SeqDesignInferenceBaiAdjustedTKK21$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			}
			inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultOLSKKCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceContinMultOLSKKCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultOLSKKIVWC")
			run_inference_checks(SeqDesignInferenceContinMultOLSKKIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultiKKLinIVWC")
			run_inference_checks(SeqDesignInferenceContinMultiKKLinIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultiKKLinCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceContinMultiKKLinCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultGLS")
			run_inference_checks(SeqDesignInferenceContinMultGLS$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinUnivKKGLMM")
			run_inference_checks(SeqDesignInferenceContinUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultiKKGLMM")
			run_inference_checks(SeqDesignInferenceContinMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinUnivKKRobustRegrIVWC")
			run_inference_checks(SeqDesignInferenceContinUnivKKRobustRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultiKKRobustRegrIVWC")
			run_inference_checks(SeqDesignInferenceContinMultiKKRobustRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinUnivKKRobustRegrCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceContinUnivKKRobustRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultiKKRobustRegrCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceContinMultiKKRobustRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultKKQuantileRegrIVWC")
			run_inference_checks(SeqDesignInferenceContinMultKKQuantileRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultKKQuantileRegrCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceContinMultKKQuantileRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		} else {
			inference_banner("SeqDesignInferenceContinUnivRobustRegr")
			run_inference_checks(SeqDesignInferenceContinUnivRobustRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultiRobustRegr")
			run_inference_checks(SeqDesignInferenceContinMultiRobustRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinUnivQuantileRegr")
			run_inference_checks(SeqDesignInferenceContinUnivQuantileRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultiQuantileRegr")
			run_inference_checks(SeqDesignInferenceContinMultiQuantileRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultLin")
			run_inference_checks(SeqDesignInferenceContinMultLin$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("SeqDesignInferenceContinMultOLS")
			run_inference_checks(SeqDesignInferenceContinMultOLS$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		}
	}

	if (response_type == "incidence"){
		inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
		run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (design_type == "CRD"){
			inference_banner("SeqDesignInferenceIncidExactZhang")
			run_inference_checks(SeqDesignInferenceIncidExactZhang$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		if (is_kk_design){
			inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidKKExactZhang")
			run_inference_checks(SeqDesignInferenceIncidKKExactZhang$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKClogitCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceIncidUnivKKClogitCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiKKClogitCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceIncidMultiKKClogitCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKClogitIVWC")
			run_inference_checks(SeqDesignInferenceIncidUnivKKClogitIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiKKClogitIVWC")
			run_inference_checks(SeqDesignInferenceIncidMultiKKClogitIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKGEE")
			run_inference_checks(SeqDesignInferenceIncidUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiKKGEE")
			run_inference_checks(SeqDesignInferenceIncidMultiKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKNewcombeRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidUnivKKNewcombeRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKGCompRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidUnivKKGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiKKGCompRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidMultiKKGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKGCompRiskRatio")
			run_inference_checks(SeqDesignInferenceIncidUnivKKGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiKKGCompRiskRatio")
			run_inference_checks(SeqDesignInferenceIncidMultiKKGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKModifiedPoisson")
			run_inference_checks(SeqDesignInferenceIncidUnivKKModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiKKModifiedPoisson")
			run_inference_checks(SeqDesignInferenceIncidMultiKKModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivKKGLMM")
			run_inference_checks(SeqDesignInferenceIncidUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiKKGLMM")
			run_inference_checks(SeqDesignInferenceIncidMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("SeqDesignInferenceIncidUnivLogRegr")
		run_inference_checks(SeqDesignInferenceIncidUnivLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceIncidMultiLogRegr")
		run_inference_checks(SeqDesignInferenceIncidMultiLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (!is_kk_design){
			inference_banner("SeqDesignInferenceIncidUnivMiettinenNurminenRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidUnivMiettinenNurminenRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivNewcombeRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidUnivNewcombeRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidUnivRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidMultiRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivGCompRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidUnivGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiGCompRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidMultiGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivGCompRiskRatio")
			run_inference_checks(SeqDesignInferenceIncidUnivGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiGCompRiskRatio")
			run_inference_checks(SeqDesignInferenceIncidMultiGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivModifiedPoisson")
			run_inference_checks(SeqDesignInferenceIncidUnivModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiModifiedPoisson")
			run_inference_checks(SeqDesignInferenceIncidMultiModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivLogBinomial")
			run_inference_checks(SeqDesignInferenceIncidUnivLogBinomial$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiLogBinomial")
			run_inference_checks(SeqDesignInferenceIncidMultiLogBinomial$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidUnivBinomialIdentityRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidUnivBinomialIdentityRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceIncidMultiBinomialIdentityRiskDiff")
			run_inference_checks(SeqDesignInferenceIncidMultiBinomialIdentityRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
	}

	if (response_type == "proportion"){
		inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
		run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceAllSimpleWilcox")
		run_inference_checks(SeqDesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropUnivKKGEE")
			run_inference_checks(SeqDesignInferencePropUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropMultiKKGEE")
			run_inference_checks(SeqDesignInferencePropMultiKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropUnivKKGLMM")
			run_inference_checks(SeqDesignInferencePropUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropMultiKKGLMM")
			run_inference_checks(SeqDesignInferencePropMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropMultiKKQuantileRegrIVWC")
			run_inference_checks(SeqDesignInferencePropMultiKKQuantileRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropMultiKKQuantileRegrCombinedLikelihood")
			run_inference_checks(SeqDesignInferencePropMultiKKQuantileRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("SeqDesignInferencePropUniBetaRegr")
		run_inference_checks(SeqDesignInferencePropUniBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferencePropMultiBetaRegr")
		run_inference_checks(SeqDesignInferencePropMultiBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (!is_kk_design){
			inference_banner("SeqDesignInferencePropUniZeroOneInflatedBetaRegr")
			run_inference_checks(SeqDesignInferencePropUniZeroOneInflatedBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropMultiZeroOneInflatedBetaRegr")
			run_inference_checks(SeqDesignInferencePropMultiZeroOneInflatedBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropUniGCompMeanDiff")
			run_inference_checks(SeqDesignInferencePropUniGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropMultiGCompMeanDiff")
			run_inference_checks(SeqDesignInferencePropMultiGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropUniFractionalLogit")
			run_inference_checks(SeqDesignInferencePropUniFractionalLogit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferencePropMultiFractionalLogit")
			run_inference_checks(SeqDesignInferencePropMultiFractionalLogit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
	}

	if (response_type == "count"){
		inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
		run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceAllSimpleWilcox")
		run_inference_checks(SeqDesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonUnivKKGEE")
			run_inference_checks(SeqDesignInferenceCountPoissonUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonMultiKKGEE")
			run_inference_checks(SeqDesignInferenceCountPoissonMultiKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonUnivKKCPoissonCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceCountPoissonUnivKKCPoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonMultiKKCPoissonCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceCountPoissonMultiKKCPoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivKKHurdlePoissonCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceCountUnivKKHurdlePoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiKKHurdlePoissonCombinedLikelihood")
			run_inference_checks(SeqDesignInferenceCountMultiKKHurdlePoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivKKHurdlePoissonIVWC")
			run_inference_checks(SeqDesignInferenceCountUnivKKHurdlePoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiKKHurdlePoissonIVWC")
			run_inference_checks(SeqDesignInferenceCountMultiKKHurdlePoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonUnivKKCPoissonIVWC")
			run_inference_checks(SeqDesignInferenceCountPoissonUnivKKCPoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonMultiKKCPoissonIVWC")
			run_inference_checks(SeqDesignInferenceCountPoissonMultiKKCPoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonUnivKKGLMM")
			run_inference_checks(SeqDesignInferenceCountPoissonUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountPoissonMultiKKGLMM")
			run_inference_checks(SeqDesignInferenceCountPoissonMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		if (!is_kk_design){
			inference_banner("SeqDesignInferenceCountUnivPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountUnivPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountMultiPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivRobustPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountUnivRobustPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiRobustPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountMultiRobustPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivQuasiPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountUnivQuasiPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiQuasiPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountMultiQuasiPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivZeroInflatedPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountUnivZeroInflatedPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiZeroInflatedPoissonRegr")
			run_inference_checks(SeqDesignInferenceCountMultiZeroInflatedPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivZeroInflatedNegBinRegr")
			run_inference_checks(SeqDesignInferenceCountUnivZeroInflatedNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiZeroInflatedNegBinRegr")
			run_inference_checks(SeqDesignInferenceCountMultiZeroInflatedNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivHurdlePoissonRegr")
			run_inference_checks(SeqDesignInferenceCountUnivHurdlePoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiHurdlePoissonRegr")
			run_inference_checks(SeqDesignInferenceCountMultiHurdlePoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountUnivHurdleNegBinRegr")
			run_inference_checks(SeqDesignInferenceCountUnivHurdleNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceCountMultiHurdleNegBinRegr")
			run_inference_checks(SeqDesignInferenceCountMultiHurdleNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("SeqDesignInferenceCountUnivNegBinRegr")
		run_inference_checks(SeqDesignInferenceCountUnivNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceCountMultiNegBinRegr")
		run_inference_checks(SeqDesignInferenceCountMultiNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}

	if (response_type == "survival"){
		if (is_kk_design){
			is_censoring_skip_error = function(msg){
				grepl("only available for uncensored", msg, fixed = TRUE) ||
				grepl("does not currently support censored", msg, fixed = TRUE) ||
				grepl("does not support censored", msg, fixed = TRUE)
			}
			for (kk_surv_class in list(SeqDesignInferenceAllKKCompoundMeanDiff, SeqDesignInferenceAllKKWilcoxIVWC, SeqDesignInferenceAllKKWilcoxRegrUnivIVWC, SeqDesignInferenceAllKKWilcoxRegrMultiIVWC)){
				class_name = kk_surv_class$classname
				inference_banner(class_name)
				err_msg = tryCatch({
					run_inference_checks(kk_surv_class$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
					NULL
				}, error = function(e) if (length(e$message) == 0L) "" else e$message)
				if (!is.null(err_msg)){
					if (is_censoring_skip_error(err_msg)) message("  Skipping ", class_name, " (censored data): ", err_msg)
					else stop(err_msg)
				}
			}
			for (kk_surv_class in list(
				SeqDesignInferenceSurvivalUnivKKGEE,
				SeqDesignInferenceSurvivalMultiKKGammaGEE,
				SeqDesignInferenceSurvivalUnivKKGLMM,
				SeqDesignInferenceSurvivalMultiKKGammaGLMM,
				SeqDesignInferenceSurvivalUnivKKClaytonCopulaIVWC,
				SeqDesignInferenceSurvivalMultiKKClaytonCopulaIVWC,
				SeqDesignInferenceSurvivalUnivKKLWACoxIVWC,
				SeqDesignInferenceSurvivalMultiKKLWACoxIVWC,
				SeqDesignInferenceSurvivalUnivKKStratCoxIVWC,
				SeqDesignInferenceSurvivalMultiKKStratCoxIVWC,
				SeqDesignInferenceSurvivalUnivKKRankRegrIVWC,
				SeqDesignInferenceSurvivalMultiKKRankRegrIVWC,
				SeqDesignInferenceSurvivalUnivKKWeibullFrailtyIVWC,
				SeqDesignInferenceSurvivalMultiKKWeibullFrailtyIVWC,
				SeqDesignInferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood,
				SeqDesignInferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood,
				SeqDesignInferenceSurvivalUnivKKLWACoxCombinedLikelihood,
				SeqDesignInferenceSurvivalMultiKKLWACoxCombinedLikelihood,
				SeqDesignInferenceSurvivalUnivKKStratCoxCombinedLikelihood,
				SeqDesignInferenceSurvivalMultiKKStratCoxCombinedLikelihood,
				SeqDesignInferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood,
				SeqDesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood
			)){
				class_name = kk_surv_class$classname
				inference_banner(class_name)
				err_msg = tryCatch({
					run_inference_checks(kk_surv_class$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
					NULL
				}, error = function(e) if (length(e$message) == 0L) "" else e$message)
				if (!is.null(err_msg)){
					if (grepl("only available for uncensored", err_msg, fixed = TRUE)) message("  Skipping ", class_name, " (censored data): ", err_msg)
					else stop(err_msg)
				}
			}
		}
		inference_banner("SeqDesignInferenceAllSimpleWilcox")
		err_msg_sw = tryCatch({
			run_inference_checks(SeqDesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			NULL
		}, error = function(e) if (length(e$message) == 0L) "" else e$message)
		if (!is.null(err_msg_sw)){
			if (grepl("does not support censored", err_msg_sw, fixed = TRUE)) message("  Skipping SeqDesignInferenceAllSimpleWilcox (censored data): ", err_msg_sw)
			else stop(err_msg_sw)
		}
		inference_banner("SeqDesignInferenceSurvivalGehanWilcox")
		run_inference_checks(SeqDesignInferenceSurvivalGehanWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalLogRank")
		run_inference_checks(SeqDesignInferenceSurvivalLogRank$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalRestrictedMeanDiff")
		run_inference_checks(SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalKMDiff")
		run_inference_checks(SeqDesignInferenceSurvivalKMDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalUniWeibullRegr")
		run_inference_checks(SeqDesignInferenceSurvivalUniWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalMultiWeibullRegr")
		run_inference_checks(SeqDesignInferenceSurvivalMultiWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalUniDepCensTransformRegr")
		run_inference_checks(SeqDesignInferenceSurvivalUniDepCensTransformRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalMultiDepCensTransformRegr")
		run_inference_checks(SeqDesignInferenceSurvivalMultiDepCensTransformRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalUniCoxPHRegr")
		run_inference_checks(SeqDesignInferenceSurvivalUniCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalMultiCoxPHRegr")
		run_inference_checks(SeqDesignInferenceSurvivalMultiCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalUniStratCoxPHRegr")
		run_inference_checks(SeqDesignInferenceSurvivalUniStratCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceSurvivalMultiStratCoxPHRegr")
		run_inference_checks(SeqDesignInferenceSurvivalMultiStratCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}

	if (response_type == "ordinal"){
		inference_banner("SeqDesignInferenceAllSimpleMeanDiff")
		run_inference_checks(SeqDesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceAllSimpleWilcox")
		run_inference_checks(SeqDesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("SeqDesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(SeqDesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalUnivKKGEE")
			run_inference_checks(SeqDesignInferenceOrdinalUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalUnivKKGLMM")
			run_inference_checks(SeqDesignInferenceOrdinalUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalMultiKKGLMM")
			run_inference_checks(SeqDesignInferenceOrdinalMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalUnivKKGLMMProbit")
			run_inference_checks(SeqDesignInferenceOrdinalUnivKKGLMMProbit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalMultiKKGLMMProbit")
			run_inference_checks(SeqDesignInferenceOrdinalMultiKKGLMMProbit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalUnivKKCondPropOddsRegr")
			run_inference_checks(SeqDesignInferenceOrdinalUnivKKCondPropOddsRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr")
			run_inference_checks(SeqDesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalUnivKKCondContRatioRegr")
			run_inference_checks(SeqDesignInferenceOrdinalUnivKKCondContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalMultiKKCondContRatioRegr")
			run_inference_checks(SeqDesignInferenceOrdinalMultiKKCondContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalUnivKKCondAdjCatLogitRegr")
			run_inference_checks(SeqDesignInferenceOrdinalUnivKKCondAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalMultiKKCondAdjCatLogitRegr")
			run_inference_checks(SeqDesignInferenceOrdinalMultiKKCondAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("SeqDesignInferenceOrdinalPairedSignTest")
			run_inference_checks(SeqDesignInferenceOrdinalPairedSignTest$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("SeqDesignInferenceOrdinalUniAdjCatLogitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiAdjCatLogitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniCumulProbitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniCumulProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiCumulProbitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiCumulProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniOrderedProbitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniOrderedProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiOrderedProbitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiOrderedProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniStereotypeLogitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniStereotypeLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiStereotypeLogitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiStereotypeLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniStereotypeProbitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniStereotypeProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiStereotypeProbitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiStereotypeProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniPropOddsRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniPropOddsRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniPartialProportionalOddsRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniPartialProportionalOddsRegr$new(seq_des_obj, verbose = FALSE, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiPartialProportionalOddsRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiPartialProportionalOddsRegr$new(seq_des_obj, verbose = FALSE, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalPartialProportionalOdds")
		run_inference_checks(SeqDesignInferenceOrdinalPartialProportionalOdds$new(seq_des_obj, verbose = FALSE, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniCLLRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniCLLRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiCLLRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiCLLRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniGCompMeanDiff")
		run_inference_checks(SeqDesignInferenceOrdinalUniGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiGCompMeanDiff")
		run_inference_checks(SeqDesignInferenceOrdinalMultiGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalJonckheereTerpstraTest")
		run_inference_checks(SeqDesignInferenceOrdinalJonckheereTerpstraTest$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalUniCauchitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalUniCauchitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiCauchitRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiCauchitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalContRatioRegr")
		run_inference_checks(SeqDesignInferenceOrdinalContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalMultiContRatioRegr")
		run_inference_checks(SeqDesignInferenceOrdinalMultiContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("SeqDesignInferenceOrdinalRidit")
		run_inference_checks(SeqDesignInferenceOrdinalRidit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}
}


for (rep_curr in 1:Nrep) {
	log_progress(paste0("\n\n====== RUNNING REPETITION ", rep_curr, " OF ", Nrep, " ======\n\n"))
	for (beta_T_iter_curr in seq_along(beta_T_values)){
		beta_T = beta_T_values[beta_T_iter_curr]
		log_progress(paste0("\n\n=== Running tests with beta_T = [", beta_T, "] ===\n\n"))
		for (dataset_name in names(datasets_and_response_models)){
			log_progress(paste0("\n\n=== Running tests for dataset: ", dataset_name, " ===\n\n"))
			for (response_type in c("continuous", "incidence", "proportion", "count", "survival", "ordinal")) {
				if (!(response_type %in% names(datasets_and_response_models[[dataset_name]]$y_original))) {
					log_progress(paste0("\n\n  === Skipping response_type: ", response_type, " for dataset: ", dataset_name, " (not found in y_original) ===\n\n"))
					next
				}
				log_progress(paste0("\n\n  === Running tests for response_type: ", response_type, " ===\n\n"))
				for (design_type in c("CRD", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise")) {
					log_progress(paste0("\n\n    === Running tests for design: ", design_type, " ==="))
					run_tests_for_response(response_type, design_type = design_type, dataset_name = dataset_name)
					log_progress(paste0("\n\n  === Finished tests for design_type: ", design_type, " ===\n\n"))
				}
				log_progress(paste0("\n\n  === Finished tests for response_type: ", response_type, " ===\n\n"))
			}
			log_progress(paste0("\n\n  === Finished tests for dataset: ", dataset_name, " ===\n\n"))
		}
	}
}
message("\n\n----------------------All tests complete!")
write_results_if_needed(force = TRUE)

rm(list=ls())
