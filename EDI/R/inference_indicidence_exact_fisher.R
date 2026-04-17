#' Exact Fisher Incidence Inference
#'
#' @export
InferenceIncidExactFisher = R6::R6Class("InferenceIncidExactFisher",
	lock_objects = FALSE,
	inherit = InferenceExact,
	public = list(
		#' @description
		#' Initialize exact Fisher inference for incidence outcomes.
		#' @param des_obj A completed design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new \code{InferenceIncidExactFisher} object.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Compute the Fisher exact treatment estimate on the log-odds scale.
		#' @param estimate_only Ignored for this estimator.
		#' @return The treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$get_exact_fisher_log_or_estimate()
		}
	),

	private = list(
		default_exact_type = "Fisher",

		resolve_exact_type = function(type){
			if (is.null(type)) type = private$default_exact_type
			if (should_run_asserts()) {
				assertChoice(type, c("Fisher"))
			}
			type
		},

		normalize_exact_inference_args = function(type, args_for_type = NULL){
			if (should_run_asserts()) {
				assertChoice(type, c("Fisher"))
				assertList(args_for_type, null.ok = TRUE)
			}
			utils::modifyList(setNames(list(list()), type), if (is.null(args_for_type)) list() else args_for_type)
		},

		assert_exact_inference_params = function(type, args_for_type){
			if (should_run_asserts()) {
				assertChoice(type, c("Fisher"))
				assertList(args_for_type)
				if (!(type %in% names(args_for_type))) stop("args_for_type must contain a list for ", type)
			}
			args = args_for_type[[type]]
			if (should_run_asserts()) {
				assertList(args)
				assertResponseType(private$des_obj$get_response_type(), "incidence")
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts()) {
				if (!private$design_supports_exact_fisher()) {
					stop("Fisher exact inference requires iBCRD, blocking, or matching designs.")
				}
			}
			invisible(args)
		},

		compute_exact_confidence_interval_by_type = function(type, alpha, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Fisher = private$ci_exact_fisher(alpha)
			)
		},

		compute_exact_two_sided_pval_for_treatment_effect_by_type = function(type, delta, args_for_type){
			if (should_run_asserts()) {
				assertNumeric(delta, len = 1)
				private$assert_exact_inference_params(type, args_for_type)
			}
			switch(type,
				Fisher = private$pval_exact_fisher(delta)
			)
		},

		design_supports_exact_fisher = function(){
				is(private$des_obj, "DesignSeqOneByOneiBCRD") ||
					is(private$des_obj, "FixedDesigniBCRD") ||
					private$is_supported_blocking_design(private$des_obj) ||
					private$has_match_structure
		},

		is_supported_blocking_design = function(des_obj){
			is(des_obj, "FixedDesignBlocking") ||
				is(des_obj, "DesignSeqOneByOneSPBR") ||
				is(des_obj, "DesignSeqOneByOneRandomBlockSize")
		},

		pval_exact_fisher = function(delta_0){
			as.numeric(private$get_exact_fisher_htest(delta_0 = delta_0)$p.value)
		},

		ci_exact_fisher = function(alpha){
			test = private$get_exact_fisher_htest(alpha = alpha, delta_0 = 0)
			ci = log(as.numeric(test$conf.int))
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		get_exact_fisher_log_or_estimate = function(){
			log(as.numeric(private$get_exact_fisher_htest()$estimate[[1]]))
		},

		get_exact_fisher_htest = function(alpha = 0.05, delta_0 = 0){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
				assertNumeric(delta_0, len = 1)
			}
			fisher_tables = private$get_exact_fisher_tables()
			conf_level = 1 - alpha
			if (fisher_tables$n_strata == 1L) {
				return(stats::fisher.test(
					fisher_tables$table,
					alternative = "two.sided",
					or = exp(delta_0),
					conf.level = conf_level
				))
			}
			if (should_run_asserts()) {
				if (abs(delta_0) > sqrt(.Machine$double.eps)) {
					stop("Stratified Fisher exact inference only supports delta = 0.")
				}
			}
			stats::mantelhaen.test(
				fisher_tables$array,
				alternative = "two.sided",
				exact = TRUE,
				conf.level = conf_level
			)
		},

		get_exact_fisher_tables = function(){
			if (!is.null(private$cached_values$incid_exact_fisher_tables)) {
				return(private$cached_values$incid_exact_fisher_tables)
			}
			fisher_tables =
					if (private$has_match_structure) {
					private$build_exact_fisher_tables_kk()
				} else if (private$is_supported_blocking_design(private$des_obj)) {
					private$build_exact_fisher_tables_blocking()
				} else {
					private$format_exact_fisher_tables(list(private$build_exact_fisher_2x2_table(seq_len(private$n))))
				}
			private$cached_values$incid_exact_fisher_tables = fisher_tables
			fisher_tables
		},

		build_exact_fisher_tables_blocking = function(){
			strata_cols = private$des_obj_priv_int$strata_cols
			if (is.null(strata_cols) || length(strata_cols) == 0L) {
				return(private$format_exact_fisher_tables(list(private$build_exact_fisher_2x2_table(seq_len(private$n)))))
			}
			Xraw = private$des_obj_priv_int$Xraw
			strata_keys = vapply(seq_len(private$n), function(i){
				private$get_exact_fisher_strata_key(Xraw[i, ], strata_cols)
			}, character(1))
			strata_indices = split(seq_len(private$n), strata_keys)
			private$format_exact_fisher_tables(lapply(strata_indices, private$build_exact_fisher_2x2_table))
		},

		get_exact_fisher_strata_key = function(x_row, strata_cols){
			paste(vapply(strata_cols, function(col){
				val = x_row[[col]]
				if (is.na(val)) "NA" else as.character(val)
			}, character(1)), collapse = "|")
		},

		build_exact_fisher_tables_kk = function(){
			m_vec = private$des_obj_priv_int$m
			if (should_run_asserts()) {
				if (is.null(m_vec)) {
					stop("Matching structure is unavailable for Fisher exact inference.")
				}
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			table_list = list()
			matched_ids = sort(unique(m_vec[m_vec > 0L]))
			for (match_id in matched_ids) {
				table_list[[length(table_list) + 1L]] = private$build_exact_fisher_2x2_table(which(m_vec == match_id))
			}
			reservoir_indices = which(m_vec == 0L)
			if (length(reservoir_indices) > 0L) {
				table_list[[length(table_list) + 1L]] = private$build_exact_fisher_2x2_table(reservoir_indices)
			}
			private$format_exact_fisher_tables(table_list)
		},

		build_exact_fisher_2x2_table = function(indices){
			i_t = private$w[indices] == 1L
			i_c = private$w[indices] == 0L
			y_idx = private$y[indices]
			matrix(
				c(
					sum(y_idx[i_t] == 1L, na.rm = TRUE),
					sum(y_idx[i_t] == 0L, na.rm = TRUE),
					sum(y_idx[i_c] == 1L, na.rm = TRUE),
					sum(y_idx[i_c] == 0L, na.rm = TRUE)
				),
				nrow = 2,
				byrow = TRUE,
				dimnames = list(c("treated", "control"), c("case", "noncase"))
			)
		},

		format_exact_fisher_tables = function(table_list){
			table_list = Filter(function(tab) sum(tab[1, ]) > 0L && sum(tab[2, ]) > 0L, table_list)
			if (should_run_asserts()) {
				if (length(table_list) == 0L) {
					stop("Cannot compute Fisher exact inference: no informative strata are available.")
				}
			}
			if (length(table_list) == 1L) {
				return(list(n_strata = 1L, table = table_list[[1]]))
			}
			table_array = array(
				0,
				dim = c(2L, 2L, length(table_list)),
				dimnames = c(dimnames(table_list[[1]]), list(NULL))
			)
			for (k in seq_along(table_list)) {
				table_array[, , k] = table_list[[k]]
			}
			list(n_strata = length(table_list), array = table_array)
		}
	)
)
