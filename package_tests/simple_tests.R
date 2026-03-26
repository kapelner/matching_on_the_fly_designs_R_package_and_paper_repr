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
r = 351
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
			r = as.integer(r),
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
	skip_slow = is(seq_des_inf, "DesignInferencePropMultiBetaRegr") || is(seq_des_inf, "DesignInferencePropUniZeroOneInflatedBetaRegr") || is(seq_des_inf, "DesignInferencePropMultiZeroOneInflatedBetaRegr") || is(seq_des_inf, "DesignInferenceSurvivalUniDepCensTransformRegr") || is(seq_des_inf, "DesignInferenceSurvivalMultiDepCensTransformRegr") || is(seq_des_inf, "DesignInferenceSurvivalMultiWeibullRegr") || is(seq_des_inf, "DesignInferenceCountMultiNegBinRegr") || is(seq_des_inf, "DesignInferenceCountUnivZeroInflatedPoissonRegr") || is(seq_des_inf, "DesignInferenceCountMultiZeroInflatedPoissonRegr") || is(seq_des_inf, "DesignInferenceCountUnivZeroInflatedNegBinRegr") || is(seq_des_inf, "DesignInferenceCountMultiZeroInflatedNegBinRegr") || is(seq_des_inf, "DesignInferenceCountUnivHurdlePoissonRegr") || is(seq_des_inf, "DesignInferenceCountMultiHurdlePoissonRegr") || is(seq_des_inf, "DesignInferenceCountUnivHurdleNegBinRegr") || is(seq_des_inf, "DesignInferenceCountMultiHurdleNegBinRegr") || is(seq_des_inf, "DesignInferenceCountUnivKKHurdlePoissonCombinedLikelihood") || is(seq_des_inf, "DesignInferenceCountMultiKKHurdlePoissonCombinedLikelihood") || is(seq_des_inf, "DesignInferenceCountUnivKKHurdlePoissonIVWC") || is(seq_des_inf, "DesignInferenceCountMultiKKHurdlePoissonIVWC") || is(seq_des_inf, "DesignInferenceSurvivalMultiCoxPHRegr") || is(seq_des_inf, "DesignInferenceSurvivalUniCoxPHRegr") || is(seq_des_inf, "DesignInferenceSurvivalUnivKKClaytonCopulaIVWC") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKClaytonCopulaIVWC") || is(seq_des_inf, "DesignInferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKLWACoxIVWC") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKStratCoxIVWC") || is(seq_des_inf, "DesignInferenceCountPoissonMultiKKCPoissonIVWC") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKWeibullFrailtyIVWC") || is(seq_des_inf, "DesignInferenceAllKKWilcoxRegrMultiIVWC") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKRankRegrIVWC") || is(seq_des_inf, "DesignInferenceCountPoissonMultiKKCPoissonCombinedLikelihood") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKLWACoxCombinedLikelihood") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKStratCoxCombinedLikelihood") || is(seq_des_inf, "DesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "DesignInferenceContinMultLin") || is(seq_des_inf, "DesignInferenceContinUnivQuantileRegr") || is(seq_des_inf, "DesignInferenceContinMultiQuantileRegr") || is(seq_des_inf, "DesignInferenceContinMultiKKLinIVWC") || is(seq_des_inf, "DesignInferenceContinMultiKKLinCombinedLikelihood")
	skip_bootstrap = is(seq_des_inf, "DesignInferenceAbstractKKGEE") || is(seq_des_inf, "DesignInferenceAbstractKKGLMM") || is(seq_des_inf, "DesignInferenceContinMultGLS") || is(seq_des_inf, "DesignInferenceAbstractKKClaytonCopulaIVWC") || is(seq_des_inf, "DesignInferenceAbstractKKClaytonCopulaCombinedLikelihood") || is(seq_des_inf, "DesignInferenceAbstractKKWeibullFrailtyIVWC") || is(seq_des_inf, "DesignInferenceAbstractKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "DesignInferenceAllKKWilcoxIVWC") || is(seq_des_inf, "DesignInferenceAbstractKKWilcoxRegrIVWC") || is(seq_des_inf, "DesignInferenceSurvivalUnivKKRankRegrIVWC") || is(seq_des_inf, "DesignInferenceIncidExactZhangAbstract") || is(seq_des_inf, "DesignInferenceAllSimpleWilcox") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKGEE") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKGLMM") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKGLMM") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKGLMMProbit") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKGLMMProbit") || is(seq_des_inf, "DesignInferenceOrdinalPairedSignTest") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKCondContRatioRegr") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKCondAdjCatLogitRegr") || is(seq_des_inf, "DesignInferenceOrdinalUniGCompMeanDiff") || is(seq_des_inf, "DesignInferenceOrdinalMultiGCompMeanDiff") || is(seq_des_inf, "DesignInferenceOrdinalMultiCLLRegr") || is(seq_des_inf, "DesignInferenceOrdinalUniOrderedProbitRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiOrderedProbitRegr") || is(seq_des_inf, "DesignInferenceOrdinalUniCauchitRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiCauchitRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiContRatioRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKGEE") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKCondContRatioRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKCondAdjCatLogitRegr")
	skip_rand      = is(seq_des_inf, "DesignInferenceAbstractKKGEE") || is(seq_des_inf, "DesignInferenceAbstractKKGLMM") || is(seq_des_inf, "DesignInferenceIncidExactZhangAbstract") || is(seq_des_inf, "DesignInferencePropUniGCompMeanDiff") || is(seq_des_inf, "DesignInferencePropMultiGCompMeanDiff") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKGEE") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKGLMM") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKGLMM") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKGLMMProbit") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKGLMMProbit") || is(seq_des_inf, "DesignInferenceOrdinalPairedSignTest") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKCondContRatioRegr") || is(seq_des_inf, "DesignInferenceOrdinalUnivKKCondAdjCatLogitRegr") || is(seq_des_inf, "DesignInferenceOrdinalUniGCompMeanDiff") || is(seq_des_inf, "DesignInferenceOrdinalMultiGCompMeanDiff") || is(seq_des_inf, "DesignInferenceOrdinalMultiCLLRegr") || is(seq_des_inf, "DesignInferenceOrdinalUniOrderedProbitRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiOrderedProbitRegr") || is(seq_des_inf, "DesignInferenceOrdinalUniCauchitRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiCauchitRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiContRatioRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKGEE") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKCondContRatioRegr") || is(seq_des_inf, "DesignInferenceOrdinalMultiKKCondAdjCatLogitRegr")
	skip_mle_pval  = is(seq_des_inf, "DesignInferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood")
	skip_rand_pval = is(seq_des_inf, "DesignInferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood") || is(seq_des_inf, "DesignInferenceContinMultGLS") || is(seq_des_inf, "DesignInferencePropUniGCompMeanDiff") || is(seq_des_inf, "DesignInferencePropMultiGCompMeanDiff")
	skip_ci_rand   = is(seq_des_inf, "DesignInferenceContinMultKKQuantileRegrIVWC") || is(seq_des_inf, "DesignInferencePropMultiKKQuantileRegrIVWC") || is(seq_des_inf, "DesignInferenceContinMultKKQuantileRegrCombinedLikelihood") || is(seq_des_inf, "DesignInferencePropMultiKKQuantileRegrCombinedLikelihood") || is(seq_des_inf, "DesignInferencePropUniGCompMeanDiff") || is(seq_des_inf, "DesignInferencePropMultiGCompMeanDiff")
	skip_ci_rand_custom = is(seq_des_inf, "DesignInferenceContinMultiKKGLMM") || is(seq_des_inf, "DesignInferenceContinUnivKKGLMM")
	
	skip_ci = beta_T == 1 && (
		is(seq_des_inf, "DesignInferenceIncidMultiLogRegr") ||
		is(seq_des_inf, "DesignInferencePropUniBetaRegr") ||
		is(seq_des_inf, "DesignInferencePropMultiBetaRegr") ||
		is(seq_des_inf, "DesignInferenceSurvivalUniCoxPHRegr") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiCoxPHRegr") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKLWACoxIVWC") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKStratCoxIVWC") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKClaytonCopulaIVWC") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKLWACoxCombinedLikelihood") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKStratCoxCombinedLikelihood") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKWeibullFrailtyIVWC") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood") ||
		is(seq_des_inf, "DesignInferenceAllKKWilcoxRegrMultiIVWC") ||
		is(seq_des_inf, "DesignInferenceSurvivalMultiKKRankRegrIVWC")
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
					 grepl("must implement", msg, fixed = TRUE) ||
					 grepl("Zhang incidence inference is only supported", msg, fixed = TRUE) ||
					 grepl("This type of inference is only available for incidence", msg, fixed = TRUE) ||
						 grepl("not enough discordant pairs", msg, ignore.case = TRUE) ||
						 grepl("Degenerate confidence interval", msg, fixed = TRUE) ||
						 grepl("inconsistent estimator units", msg, ignore.case = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned NA bounds", msg, fixed = TRUE) ||
						 					 grepl("Bootstrap confidence interval returned non-finite bounds", msg, fixed = TRUE) ||
						 					 						 grepl("Weibull regression failed to converge", msg, fixed = TRUE) ||
						 					 						 grepl("Invalid output detected", msg, fixed = TRUE) ||
						 					 						 grepl("missing value where TRUE/FALSE needed", msg, fixed = TRUE) ||
						 					 						 ((grepl("NA/NaN/Inf", msg, fixed = TRUE) || grepl("non-finite standard error", msg, fixed = TRUE) || grepl("could not compute a finite standard error", msg, fixed = TRUE)) &&
						 					 
						 													 (is(seq_des_inf, "DesignInferenceAbstractKKClogitIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKPoissonCPoissonIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKStratCoxIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKWeibullFrailtyIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKWilcoxRegrIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKSurvivalRankRegrIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKGEE") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKGLMM") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKRobustRegrIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKRobustRegrCombinedLikelihood") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKQuantileRegrIVWC") ||
						 													  is(seq_des_inf, "DesignInferenceAbstractKKQuantileRegrCombinedLikelihood") ||
						 													  													  is(seq_des_inf, "DesignInferencePropUniFractionalLogit") ||
						 													  													  is(seq_des_inf, "DesignInferenceIncidUnivRiskDiff") ||
						 													  													  is(seq_des_inf, "DesignInferenceIncidMultiRiskDiff") ||
						 													  													  is(seq_des_inf, "DesignInferenceIncidUnivGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "DesignInferenceIncidMultiGCompRiskDiff") ||
						 													  													  is(seq_des_inf, "DesignInferenceIncidUnivGCompRiskRatio") ||
						 													  													  is(seq_des_inf, "DesignInferenceIncidMultiGCompRiskRatio") ||
						 													  is(seq_des_inf, "DesignInferenceIncidUnivModifiedPoisson") ||
						 													  is(seq_des_inf, "DesignInferenceIncidMultiModifiedPoisson") ||
						 													  is(seq_des_inf, "DesignInferencePropUniGCompMeanDiff") ||
						 													  is(seq_des_inf, "DesignInferencePropMultiGCompMeanDiff") ||
						 													  is(seq_des_inf, "DesignInferencePropUniZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "DesignInferencePropMultiZeroOneInflatedBetaRegr") ||
						 													  is(seq_des_inf, "DesignInferenceContinMultGLS")))
						 													  										if (isTRUE(is_non_fatal)){
				message("Skipping ", label, " (non-fatal): ", e$message)
				duration_time_sec = unname(proc.time()[["elapsed"]]) - start_elapsed
				record_result(dataset_name, dataset_n_rows, dataset_n_cols, response_type, design_type, class(seq_des_inf)[1], label, NA_character_, status = "error", duration_time_sec = duration_time_sec, error_message = e$message)
			} else {
				stop(e$message)
			}
		})
	}

	if (is(seq_des_inf, "DesignInferenceOrdinalJonckheereTerpstraTest")){
		message("    Calling compute_exact_two_sided_pval_for_treatment_effect()")
		safe_call("compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		message("    Calling compute_treatment_estimate()")
		safe_call("compute_treatment_estimate", seq_des_inf$compute_treatment_estimate())
		return(invisible(NULL))
	}

	if (response_type == "incidence"){
		message("    Calling compute_exact_two_sided_pval_for_treatment_effect()")
		safe_call("compute_exact_two_sided_pval_for_treatment_effect", seq_des_inf$compute_exact_two_sided_pval_for_treatment_effect())
		message("    Calling compute_exact_confidence_interval()")
		safe_call("compute_exact_confidence_interval", seq_des_inf$compute_exact_confidence_interval(args_for_type = list(Zhang = list(combination_method = "Fisher", pval_epsilon = pval_epsilon))))
	}

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
		safe_call("compute_bootstrap_confidence_interval", seq_des_inf$compute_bootstrap_confidence_interval(B = r, na.rm = TRUE))
	}
	if (!skip_slow && !skip_bootstrap){
		message("    Calling compute_bootstrap_two_sided_pval()")
		safe_call("compute_bootstrap_two_sided_pval", seq_des_inf$compute_bootstrap_two_sided_pval(B = r, na.rm = TRUE))
	}
	if (!skip_slow && !skip_rand && !skip_rand_pval && response_type %in% c("continuous", "survival", "proportion")){
		message("    Calling compute_two_sided_pval_for_treatment_effect_rand()")
		safe_call("compute_two_sided_pval_for_treatment_effect_rand", seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(r = r, show_progress = FALSE))
		message("    Calling compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)")
		use_log   = response_type == "survival"
		use_logit = response_type == "proportion"
		transform_for_rand = if (use_log) "log" else if (use_logit) "logit" else "none"
		delta_for_rand = 0.5
		safe_call("compute_two_sided_pval_for_treatment_effect_rand(delta=0.5)",
				seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(r = r, delta = delta_for_rand, transform_responses = transform_for_rand, show_progress = FALSE))
	}

	if (!skip_slow && !skip_rand && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous", "proportion", "count")){
		message("    Calling compute_confidence_interval_rand()")
		safe_call("compute_confidence_interval_rand", seq_des_inf$compute_confidence_interval_rand(r = r, pval_epsilon = pval_epsilon, show_progress = FALSE))
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
			safe_call("compute_two_sided_pval_for_treatment_effect_rand(custom)", seq_des_inf$compute_two_sided_pval_for_treatment_effect_rand(r = r, show_progress = FALSE))
		}
		if (!skip_slow && !skip_ci && !skip_ci_rand && test_compute_confidence_interval_rand && response_type %in% c("continuous")){
			message("    Calling compute_confidence_interval_rand(custom)")
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

	seq_des_obj = tryCatch(switch(design_type,
		KK21 =         SeqDesignKK21$new(        response_type = response_type, n = n),
		KK21stepwise = SeqDesignKK21stepwise$new(response_type = response_type, n = n),
		KK14 =         SeqDesignKK14$new(        response_type = response_type, n = n),
		Bernoulli =    SeqDesignBernoulli$new(   response_type = response_type, n = n),
		Efron =        SeqDesignEfron$new(       response_type = response_type, n = n),
		Atkinson =     SeqDesignAtkinson$new(    response_type = response_type, n = n),
		iBCRD =        SeqDesigniBCRD$new(       response_type = response_type, n = n),
		SPBR =         SeqDesignSPBR$new(        strata_cols = names(X_design)[1:min(2, ncol(X_design))], block_size = 4, response_type = response_type, n = n),
		FixedBernoulli = FixedDesignBernoulli$new( response_type = response_type, n = n),
		FixediBCRD =     FixedDesigniBCRD$new(     response_type = response_type, n = n),
		FixedBlocking =  FixedDesignBlocking$new(  strata_cols = names(X_design)[1:min(2, ncol(X_design))], response_type = response_type, n = n),
		FixedCluster =   FixedDesignCluster$new(   cluster_col = names(X_design)[1], response_type = response_type, n = n),
		FixedBlockedCluster = FixedDesignBlockedCluster$new( strata_cols = names(X_design)[2:min(2, ncol(X_design))], cluster_col = names(X_design)[1], response_type = response_type, n = n),
		FixedBinaryMatch = FixedDesignBinaryMatch$new( response_type = response_type, n = n),
		FixedGreedy =    FixedDesignGreedy$new(    response_type = response_type, n = n),
		FixedRerandomization = FixedDesignRerandomization$new( response_type = response_type, n = n),
		FixedMatchingGreedy = FixedDesignMatchingGreedyPairSwitching$new( response_type = response_type, n = n),
		stop("Unsupported design_type: ", design_type)
	), error = function(e){ message("    Skipping design (creation error): ", e$message); NULL })
	if (is.null(seq_des_obj)) return(invisible(NULL))

	if (inherits(seq_des_obj, "SeqDesign")){
		for (t in 1 : n){
			w_t = seq_des_obj$add_subject_to_experiment_and_assign(X_design[t, , drop = FALSE])
			y_t = apply_treatment_effect_and_noise(y[t], w_t, response_type)
			seq_des_obj$add_subject_response(t, y_t, dead[t])
		}
	} else {
		# It is a FixedDesign but not a SeqDesign
		for (t in 1 : n){
			seq_des_obj$add_subject(X_design[t, , drop = FALSE])
		}
		randomize_ok = tryCatch({ seq_des_obj$randomize(); TRUE }, error = function(e){ message("    Skipping design: ", e$message); FALSE })
		if (!randomize_ok) return(invisible(NULL))
		w = seq_des_obj$get_w()
		for (t in 1 : n){
			y_t = apply_treatment_effect_and_noise(y[t], w[t], response_type)
			seq_des_obj$add_subject_response(t, y_t, dead[t])
		}
	}

	is_kk_design = design_type %in% c("KK21", "KK21stepwise", "KK14")
	if (response_type == "continuous"){
		inference_banner("DesignInferenceAllSimpleMeanDiff")
		run_inference_checks(DesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		inference_banner("DesignInferenceAllSimpleWilcox")
		run_inference_checks(DesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		if (is_kk_design){
			if (design_type == "KK14"){
				inference_banner("DesignInferenceBaiAdjustedTKK14")
				run_inference_checks(DesignInferenceBaiAdjustedTKK14$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			}
			if (design_type %in% c("KK21", "KK21stepwise")){
				inference_banner("DesignInferenceBaiAdjustedTKK21")
				run_inference_checks(DesignInferenceBaiAdjustedTKK21$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			}
			inference_banner("DesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(DesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultOLSKKCombinedLikelihood")
			run_inference_checks(DesignInferenceContinMultOLSKKCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultOLSKKIVWC")
			run_inference_checks(DesignInferenceContinMultOLSKKIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultiKKLinIVWC")
			run_inference_checks(DesignInferenceContinMultiKKLinIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultiKKLinCombinedLikelihood")
			run_inference_checks(DesignInferenceContinMultiKKLinCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultGLS")
			run_inference_checks(DesignInferenceContinMultGLS$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinUnivKKGLMM")
			run_inference_checks(DesignInferenceContinUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultiKKGLMM")
			run_inference_checks(DesignInferenceContinMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinUnivKKRobustRegrIVWC")
			run_inference_checks(DesignInferenceContinUnivKKRobustRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultiKKRobustRegrIVWC")
			run_inference_checks(DesignInferenceContinMultiKKRobustRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinUnivKKRobustRegrCombinedLikelihood")
			run_inference_checks(DesignInferenceContinUnivKKRobustRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultiKKRobustRegrCombinedLikelihood")
			run_inference_checks(DesignInferenceContinMultiKKRobustRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultKKQuantileRegrIVWC")
			run_inference_checks(DesignInferenceContinMultKKQuantileRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultKKQuantileRegrCombinedLikelihood")
			run_inference_checks(DesignInferenceContinMultKKQuantileRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		} else {
			inference_banner("DesignInferenceContinUnivRobustRegr")
			run_inference_checks(DesignInferenceContinUnivRobustRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultiRobustRegr")
			run_inference_checks(DesignInferenceContinMultiRobustRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinUnivQuantileRegr")
			run_inference_checks(DesignInferenceContinUnivQuantileRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultiQuantileRegr")
			run_inference_checks(DesignInferenceContinMultiQuantileRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultLin")
			run_inference_checks(DesignInferenceContinMultLin$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
			inference_banner("DesignInferenceContinMultOLS")
			run_inference_checks(DesignInferenceContinMultOLS$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, dataset_n_rows, dataset_n_cols)
		}
	}

	if (response_type == "incidence"){
		inference_banner("DesignInferenceAllSimpleMeanDiff")
		run_inference_checks(DesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("DesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(DesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKClogitCombinedLikelihood")
			run_inference_checks(DesignInferenceIncidUnivKKClogitCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiKKClogitCombinedLikelihood")
			run_inference_checks(DesignInferenceIncidMultiKKClogitCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKClogitIVWC")
			run_inference_checks(DesignInferenceIncidUnivKKClogitIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiKKClogitIVWC")
			run_inference_checks(DesignInferenceIncidMultiKKClogitIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKGEE")
			run_inference_checks(DesignInferenceIncidUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiKKGEE")
			run_inference_checks(DesignInferenceIncidMultiKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKNewcombeRiskDiff")
			run_inference_checks(DesignInferenceIncidUnivKKNewcombeRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKGCompRiskDiff")
			run_inference_checks(DesignInferenceIncidUnivKKGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiKKGCompRiskDiff")
			run_inference_checks(DesignInferenceIncidMultiKKGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKGCompRiskRatio")
			run_inference_checks(DesignInferenceIncidUnivKKGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiKKGCompRiskRatio")
			run_inference_checks(DesignInferenceIncidMultiKKGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKModifiedPoisson")
			run_inference_checks(DesignInferenceIncidUnivKKModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiKKModifiedPoisson")
			run_inference_checks(DesignInferenceIncidMultiKKModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivKKGLMM")
			run_inference_checks(DesignInferenceIncidUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiKKGLMM")
			run_inference_checks(DesignInferenceIncidMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("DesignInferenceIncidUnivLogRegr")
		run_inference_checks(DesignInferenceIncidUnivLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceIncidMultiLogRegr")
		run_inference_checks(DesignInferenceIncidMultiLogRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (!is_kk_design){
			inference_banner("DesignInferenceIncidUnivMiettinenNurminenRiskDiff")
			run_inference_checks(DesignInferenceIncidUnivMiettinenNurminenRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivNewcombeRiskDiff")
			run_inference_checks(DesignInferenceIncidUnivNewcombeRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivRiskDiff")
			run_inference_checks(DesignInferenceIncidUnivRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiRiskDiff")
			run_inference_checks(DesignInferenceIncidMultiRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivGCompRiskDiff")
			run_inference_checks(DesignInferenceIncidUnivGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiGCompRiskDiff")
			run_inference_checks(DesignInferenceIncidMultiGCompRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivGCompRiskRatio")
			run_inference_checks(DesignInferenceIncidUnivGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiGCompRiskRatio")
			run_inference_checks(DesignInferenceIncidMultiGCompRiskRatio$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivModifiedPoisson")
			run_inference_checks(DesignInferenceIncidUnivModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiModifiedPoisson")
			run_inference_checks(DesignInferenceIncidMultiModifiedPoisson$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivLogBinomial")
			run_inference_checks(DesignInferenceIncidUnivLogBinomial$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiLogBinomial")
			run_inference_checks(DesignInferenceIncidMultiLogBinomial$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidUnivBinomialIdentityRiskDiff")
			run_inference_checks(DesignInferenceIncidUnivBinomialIdentityRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceIncidMultiBinomialIdentityRiskDiff")
			run_inference_checks(DesignInferenceIncidMultiBinomialIdentityRiskDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
	}

	if (response_type == "proportion"){
		inference_banner("DesignInferenceAllSimpleMeanDiff")
		run_inference_checks(DesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceAllSimpleWilcox")
		run_inference_checks(DesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("DesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(DesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropUnivKKGEE")
			run_inference_checks(DesignInferencePropUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropMultiKKGEE")
			run_inference_checks(DesignInferencePropMultiKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropUnivKKGLMM")
			run_inference_checks(DesignInferencePropUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropMultiKKGLMM")
			run_inference_checks(DesignInferencePropMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropMultiKKQuantileRegrIVWC")
			run_inference_checks(DesignInferencePropMultiKKQuantileRegrIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropMultiKKQuantileRegrCombinedLikelihood")
			run_inference_checks(DesignInferencePropMultiKKQuantileRegrCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("DesignInferencePropUniBetaRegr")
		run_inference_checks(DesignInferencePropUniBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferencePropMultiBetaRegr")
		run_inference_checks(DesignInferencePropMultiBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (!is_kk_design){
			inference_banner("DesignInferencePropUniZeroOneInflatedBetaRegr")
			run_inference_checks(DesignInferencePropUniZeroOneInflatedBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropMultiZeroOneInflatedBetaRegr")
			run_inference_checks(DesignInferencePropMultiZeroOneInflatedBetaRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropUniGCompMeanDiff")
			run_inference_checks(DesignInferencePropUniGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropMultiGCompMeanDiff")
			run_inference_checks(DesignInferencePropMultiGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropUniFractionalLogit")
			run_inference_checks(DesignInferencePropUniFractionalLogit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferencePropMultiFractionalLogit")
			run_inference_checks(DesignInferencePropMultiFractionalLogit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
	}

	if (response_type == "count"){
		inference_banner("DesignInferenceAllSimpleMeanDiff")
		run_inference_checks(DesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceAllSimpleWilcox")
		run_inference_checks(DesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("DesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(DesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonUnivKKGEE")
			run_inference_checks(DesignInferenceCountPoissonUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonMultiKKGEE")
			run_inference_checks(DesignInferenceCountPoissonMultiKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonUnivKKCPoissonCombinedLikelihood")
			run_inference_checks(DesignInferenceCountPoissonUnivKKCPoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonMultiKKCPoissonCombinedLikelihood")
			run_inference_checks(DesignInferenceCountPoissonMultiKKCPoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivKKHurdlePoissonCombinedLikelihood")
			run_inference_checks(DesignInferenceCountUnivKKHurdlePoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiKKHurdlePoissonCombinedLikelihood")
			run_inference_checks(DesignInferenceCountMultiKKHurdlePoissonCombinedLikelihood$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivKKHurdlePoissonIVWC")
			run_inference_checks(DesignInferenceCountUnivKKHurdlePoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiKKHurdlePoissonIVWC")
			run_inference_checks(DesignInferenceCountMultiKKHurdlePoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonUnivKKCPoissonIVWC")
			run_inference_checks(DesignInferenceCountPoissonUnivKKCPoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonMultiKKCPoissonIVWC")
			run_inference_checks(DesignInferenceCountPoissonMultiKKCPoissonIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonUnivKKGLMM")
			run_inference_checks(DesignInferenceCountPoissonUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountPoissonMultiKKGLMM")
			run_inference_checks(DesignInferenceCountPoissonMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		if (!is_kk_design){
			inference_banner("DesignInferenceCountUnivPoissonRegr")
			run_inference_checks(DesignInferenceCountUnivPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiPoissonRegr")
			run_inference_checks(DesignInferenceCountMultiPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivRobustPoissonRegr")
			run_inference_checks(DesignInferenceCountUnivRobustPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiRobustPoissonRegr")
			run_inference_checks(DesignInferenceCountMultiRobustPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivQuasiPoissonRegr")
			run_inference_checks(DesignInferenceCountUnivQuasiPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiQuasiPoissonRegr")
			run_inference_checks(DesignInferenceCountMultiQuasiPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivZeroInflatedPoissonRegr")
			run_inference_checks(DesignInferenceCountUnivZeroInflatedPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiZeroInflatedPoissonRegr")
			run_inference_checks(DesignInferenceCountMultiZeroInflatedPoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivZeroInflatedNegBinRegr")
			run_inference_checks(DesignInferenceCountUnivZeroInflatedNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiZeroInflatedNegBinRegr")
			run_inference_checks(DesignInferenceCountMultiZeroInflatedNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivHurdlePoissonRegr")
			run_inference_checks(DesignInferenceCountUnivHurdlePoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiHurdlePoissonRegr")
			run_inference_checks(DesignInferenceCountMultiHurdlePoissonRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountUnivHurdleNegBinRegr")
			run_inference_checks(DesignInferenceCountUnivHurdleNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceCountMultiHurdleNegBinRegr")
			run_inference_checks(DesignInferenceCountMultiHurdleNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("DesignInferenceCountUnivNegBinRegr")
		run_inference_checks(DesignInferenceCountUnivNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceCountMultiNegBinRegr")
		run_inference_checks(DesignInferenceCountMultiNegBinRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}

	if (response_type == "survival"){
		if (is_kk_design){
			is_censoring_skip_error = function(msg){
				grepl("only available for uncensored", msg, fixed = TRUE) ||
				grepl("does not currently support censored", msg, fixed = TRUE) ||
				grepl("does not support censored", msg, fixed = TRUE)
			}
			for (kk_surv_class in list(DesignInferenceAllKKCompoundMeanDiff, DesignInferenceAllKKWilcoxIVWC, DesignInferenceAllKKWilcoxRegrUnivIVWC, DesignInferenceAllKKWilcoxRegrMultiIVWC)){
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
				DesignInferenceSurvivalUnivKKGEE,
				DesignInferenceSurvivalMultiKKGammaGEE,
				DesignInferenceSurvivalUnivKKGLMM,
				DesignInferenceSurvivalMultiKKGammaGLMM,
				DesignInferenceSurvivalUnivKKClaytonCopulaIVWC,
				DesignInferenceSurvivalMultiKKClaytonCopulaIVWC,
				DesignInferenceSurvivalUnivKKLWACoxIVWC,
				DesignInferenceSurvivalMultiKKLWACoxIVWC,
				DesignInferenceSurvivalUnivKKStratCoxIVWC,
				DesignInferenceSurvivalMultiKKStratCoxIVWC,
				DesignInferenceSurvivalUnivKKRankRegrIVWC,
				DesignInferenceSurvivalMultiKKRankRegrIVWC,
				DesignInferenceSurvivalUnivKKWeibullFrailtyIVWC,
				DesignInferenceSurvivalMultiKKWeibullFrailtyIVWC,
				DesignInferenceSurvivalUnivKKClaytonCopulaCombinedLikelihood,
				DesignInferenceSurvivalMultiKKClaytonCopulaCombinedLikelihood,
				DesignInferenceSurvivalUnivKKLWACoxCombinedLikelihood,
				DesignInferenceSurvivalMultiKKLWACoxCombinedLikelihood,
				DesignInferenceSurvivalUnivKKStratCoxCombinedLikelihood,
				DesignInferenceSurvivalMultiKKStratCoxCombinedLikelihood,
				DesignInferenceSurvivalUnivKKWeibullFrailtyCombinedLikelihood,
				DesignInferenceSurvivalMultiKKWeibullFrailtyCombinedLikelihood
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
		inference_banner("DesignInferenceAllSimpleWilcox")
		err_msg_sw = tryCatch({
			run_inference_checks(DesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			NULL
		}, error = function(e) if (length(e$message) == 0L) "" else e$message)
		if (!is.null(err_msg_sw)){
			if (grepl("does not support censored", err_msg_sw, fixed = TRUE)) message("  Skipping DesignInferenceAllSimpleWilcox (censored data): ", err_msg_sw)
			else stop(err_msg_sw)
		}
		inference_banner("DesignInferenceSurvivalGehanWilcox")
		run_inference_checks(DesignInferenceSurvivalGehanWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalLogRank")
		run_inference_checks(DesignInferenceSurvivalLogRank$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalRestrictedMeanDiff")
		run_inference_checks(DesignInferenceSurvivalRestrictedMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalKMDiff")
		run_inference_checks(DesignInferenceSurvivalKMDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalUniWeibullRegr")
		run_inference_checks(DesignInferenceSurvivalUniWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalMultiWeibullRegr")
		run_inference_checks(DesignInferenceSurvivalMultiWeibullRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalUniDepCensTransformRegr")
		run_inference_checks(DesignInferenceSurvivalUniDepCensTransformRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalMultiDepCensTransformRegr")
		run_inference_checks(DesignInferenceSurvivalMultiDepCensTransformRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalUniCoxPHRegr")
		run_inference_checks(DesignInferenceSurvivalUniCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalMultiCoxPHRegr")
		run_inference_checks(DesignInferenceSurvivalMultiCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalUniStratCoxPHRegr")
		run_inference_checks(DesignInferenceSurvivalUniStratCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceSurvivalMultiStratCoxPHRegr")
		run_inference_checks(DesignInferenceSurvivalMultiStratCoxPHRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
	}

	if (response_type == "ordinal"){
		inference_banner("DesignInferenceAllSimpleMeanDiff")
		run_inference_checks(DesignInferenceAllSimpleMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceAllSimpleWilcox")
		run_inference_checks(DesignInferenceAllSimpleWilcox$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		if (is_kk_design){
			inference_banner("DesignInferenceAllKKCompoundMeanDiff")
			run_inference_checks(DesignInferenceAllKKCompoundMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxRegrUnivIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrUnivIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceAllKKWilcoxRegrMultiIVWC")
			run_inference_checks(DesignInferenceAllKKWilcoxRegrMultiIVWC$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalUnivKKGEE")
			run_inference_checks(DesignInferenceOrdinalUnivKKGEE$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalUnivKKGLMM")
			run_inference_checks(DesignInferenceOrdinalUnivKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalMultiKKGLMM")
			run_inference_checks(DesignInferenceOrdinalMultiKKGLMM$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalUnivKKGLMMProbit")
			run_inference_checks(DesignInferenceOrdinalUnivKKGLMMProbit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalMultiKKGLMMProbit")
			run_inference_checks(DesignInferenceOrdinalMultiKKGLMMProbit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalUnivKKCondPropOddsRegr")
			run_inference_checks(DesignInferenceOrdinalUnivKKCondPropOddsRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr")
			run_inference_checks(DesignInferenceOrdinalUnivKKCondPropOddsCombinedRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalUnivKKCondContRatioRegr")
			run_inference_checks(DesignInferenceOrdinalUnivKKCondContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalMultiKKCondContRatioRegr")
			run_inference_checks(DesignInferenceOrdinalMultiKKCondContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalUnivKKCondAdjCatLogitRegr")
			run_inference_checks(DesignInferenceOrdinalUnivKKCondAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalMultiKKCondAdjCatLogitRegr")
			run_inference_checks(DesignInferenceOrdinalMultiKKCondAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
			inference_banner("DesignInferenceOrdinalPairedSignTest")
			run_inference_checks(DesignInferenceOrdinalPairedSignTest$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		}
		inference_banner("DesignInferenceOrdinalUniAdjCatLogitRegr")
		run_inference_checks(DesignInferenceOrdinalUniAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiAdjCatLogitRegr")
		run_inference_checks(DesignInferenceOrdinalMultiAdjCatLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniCumulProbitRegr")
		run_inference_checks(DesignInferenceOrdinalUniCumulProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiCumulProbitRegr")
		run_inference_checks(DesignInferenceOrdinalMultiCumulProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniOrderedProbitRegr")
		run_inference_checks(DesignInferenceOrdinalUniOrderedProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiOrderedProbitRegr")
		run_inference_checks(DesignInferenceOrdinalMultiOrderedProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniStereotypeLogitRegr")
		run_inference_checks(DesignInferenceOrdinalUniStereotypeLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiStereotypeLogitRegr")
		run_inference_checks(DesignInferenceOrdinalMultiStereotypeLogitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniStereotypeProbitRegr")
		run_inference_checks(DesignInferenceOrdinalUniStereotypeProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiStereotypeProbitRegr")
		run_inference_checks(DesignInferenceOrdinalMultiStereotypeProbitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniPropOddsRegr")
		run_inference_checks(DesignInferenceOrdinalUniPropOddsRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniPartialProportionalOddsRegr")
		run_inference_checks(DesignInferenceOrdinalUniPartialProportionalOddsRegr$new(seq_des_obj, verbose = FALSE, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiPartialProportionalOddsRegr")
		run_inference_checks(DesignInferenceOrdinalMultiPartialProportionalOddsRegr$new(seq_des_obj, verbose = FALSE, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalPartialProportionalOdds")
		run_inference_checks(DesignInferenceOrdinalPartialProportionalOdds$new(seq_des_obj, verbose = FALSE, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniCLLRegr")
		run_inference_checks(DesignInferenceOrdinalUniCLLRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiCLLRegr")
		run_inference_checks(DesignInferenceOrdinalMultiCLLRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniGCompMeanDiff")
		run_inference_checks(DesignInferenceOrdinalUniGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiGCompMeanDiff")
		run_inference_checks(DesignInferenceOrdinalMultiGCompMeanDiff$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalJonckheereTerpstraTest")
		run_inference_checks(DesignInferenceOrdinalJonckheereTerpstraTest$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalUniCauchitRegr")
		run_inference_checks(DesignInferenceOrdinalUniCauchitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiCauchitRegr")
		run_inference_checks(DesignInferenceOrdinalMultiCauchitRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalContRatioRegr")
		run_inference_checks(DesignInferenceOrdinalContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalMultiContRatioRegr")
		run_inference_checks(DesignInferenceOrdinalMultiContRatioRegr$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
		inference_banner("DesignInferenceOrdinalRidit")
		run_inference_checks(DesignInferenceOrdinalRidit$new(seq_des_obj, num_cores = NUM_CORES), response_type, design_type, dataset_name, nrow(D$X), ncol(D$X))
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
				for (design_type in c("Bernoulli", "iBCRD", "Efron", "KK14", "KK21", "KK21stepwise", "SPBR", "FixedBernoulli", "FixediBCRD", "FixedBlocking", "FixedCluster", "FixedBlockedCluster", "FixedBinaryMatch", "FixedGreedy", "FixedRerandomization", "FixedMatchingGreedy")) {
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
