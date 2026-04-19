#' An Abstract Experimental Design
#'
#' @name Design
#' @description Internal method.
#' An abstract R6 Class encapsulating the data and functionality for an experimental design.
#' This class takes care of data storage and response handling.
#'
#' @keywords internal
#' @export
Design = R6::R6Class("Design",
	public = list(
		#' @description
		#' Initialize an experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size (if fixed).
		#' @param verbose	Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates when building the
		#'   model matrix for inference. One of:
		#'   \describe{
		#'     \item{\code{"impute"} (default)}{Missing values are filled in using random-forest
		#'       imputation (\code{missRanger}, falling back to \code{missForest} on failure).
		#'       The response vector is included as an auxiliary predictor when available.
		#'       This preserves all covariates and all subjects but introduces imputed values
		#'       that influence inference.}
		#'     \item{\code{"drop_column"}}{Any covariate column that contains at least one
		#'       missing value is dropped entirely from the model matrix before inference.
		#'       No values are invented; the remaining complete columns are used as-is.
		#'       This is conservative but transparent.}
		#'     \item{\code{"error"}}{An error is thrown as soon as any missing value is
		#'       detected in the covariate matrix. Use this when you want to guarantee that
		#'       inference runs on exactly the data you supplied, with no silent modification.}
		#'   }
		#'
		#' @return 			A new `Design` object
		initialize = function(
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = FALSE,
				n = NULL,
				verbose = FALSE,
				missingness_method = "impute"
			) {
			if (should_run_asserts()) {
				assertChoice(response_type, c("continuous", "incidence", "proportion", "count", "survival", "ordinal"))
				assertNumeric(prob_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
				assertFlag(include_is_missing_as_a_new_feature)
				assertFlag(verbose)
				assertCount(n, null.ok = TRUE)
				assertChoice(missingness_method, c("impute", "drop_column", "error"))
			}

			if (is.null(n)){
				private$fixed_sample = FALSE
			} else {
				n = as.integer(n)
				private$n = n
				private$fixed_sample = TRUE
			}

			private$prob_T = prob_T
			private$response_type = response_type
			private$include_is_missing_as_a_new_feature = include_is_missing_as_a_new_feature
			private$missingness_method = missingness_method

			# Ensure budget is respected among openmp and other packages

			private$verbose = verbose
			if (private$fixed_sample){
				private$y = 	rep(NA_real_, n)
				private$w = 	rep(NA_real_, n)
				private$dead =  rep(NA_real_, n)
			}

			if (private$verbose){
				cat(paste0("Initialized a ",
				class(self)[1],
				" experiment with response type ",
				response_type,
				" and ",
				ifelse(private$fixed_sample, "fixed sample", "not fixed sample"),
				 ".\n"))
			}
		},

		#' @description
		#' For CARA designs, add a single subject response.
		#'
		#' @param t The subject index.
		#' @param y The response value.
		#' @param dead If the response is censored (0 for survival).
		add_one_subject_response = function(t, y, dead = 1) {
			if (should_run_asserts()) {
				assertNumeric(t, len = 1)
				assertNumeric(y, len = 1)
				assertNumeric(dead, len = 1)
				assertChoice(dead, c(0, 1))
				assertCount(t, positive = TRUE)
				if (t > private$t){
					stop(paste("You cannot add response for subject", t, "when the most recent subjects' record added is", private$t))
				}
			}
			if (length(private$y) >= t & !is.na(private$y[t])){
				warning(paste("Overwriting previous response for t =", t, "y[t] =", private$y[t]))
			}
			if (private$response_type == "continuous"){
				if (should_run_asserts()) {
					assertNumeric(y, any.missing = FALSE)
				}
			} else if (private$response_type == "incidence"){
				if (should_run_asserts()) {
					assertChoice(y, c(0, 1))
				}
			} else if (private$response_type == "proportion"){
				if (should_run_asserts()) {
					assertNumeric(y, any.missing = FALSE, lower = 0, upper = 1)
				}
			} else if (private$response_type == "count"){
				if (should_run_asserts()) {
					assertCount(y, na.ok = FALSE)
				}
			} else if (private$response_type == "survival"){
				if (should_run_asserts()) {
					assertNumeric(y, any.missing = FALSE, lower = 0)
				}
				if (y == 0){
					warning("0 survival responses not allowed --- recording .Machine$double.eps as its value instead")
					y = .Machine$double.eps
				}
			} else if (private$response_type == "ordinal"){
				if (is.factor(y)){
					if (should_run_asserts()) {
						assertFactor(y, ordered = TRUE, any.missing = FALSE)
					}
					y = as.integer(y)
				} else {
					if (should_run_asserts()) {
						assertCount(y, positive = TRUE, na.ok = FALSE)
					}
				}
			}
			if (should_run_asserts()) {
				if (dead == 0 & private$response_type != "survival"){
					stop("censored observations are only available for survival response types")
				}
			}
			if (private$fixed_sample | t <= length(private$y)){
				private$y[t] = y
				private$dead[t] = dead
			} else if (t == length(private$y) + 1){
				private$y = c(private$y, y)
				private$dead = c(private$dead, dead)
			} else {
				if (should_run_asserts()) {
					stop("You cannot add a response for a subject that has not yet arrived when the sample size is not fixed in advance.")
				}
			}
			private$y_i_t_i[[t]] = private$t
		},

		#' @description
		#' For non-CARA designs, add all subject responses.
		#'
		#' @param ys The responses as a numeric vector.
		#' @param deads The binary vector indicating if dead/censored.
		add_all_subject_responses = function(ys, deads = NULL) {
			if (is.null(deads)){
				deads = rep(1, private$t)
			}
			if (should_run_asserts()) {
				if (private$response_type == "ordinal" && is.factor(ys)){
					assertFactor(ys, len = private$t, ordered = TRUE, any.missing = FALSE)
				} else {
					assertNumeric(ys, len = private$t)
				}
				assertNumeric(deads, len = private$t)
				
				if (private$response_type != "survival" && any(deads == 0)){
					stop("censored observations are only available for survival response types")
				}
			}
			
			if (private$response_type == "ordinal" && is.factor(ys)){
				ys = as.integer(ys)
			}

			private$y = as.numeric(ys)
			private$dead = as.numeric(deads)
			private$y_i_t_i = as.list(seq_len(private$t))
		},

		#' @description
		#' For analysis on already-completed experimental data
		#'
		#' @param w The binary responses.
		overwrite_all_subject_assignments = function(w) {
			if (should_run_asserts()) {
				assertIntegerish(w, lower = 0, upper = 1, any.missing = FALSE, len = private$t)
			}
			private$w = w
		},

		#' @description
		#' Check if this design was initialized with a fixed sample size n
		#'
		#' @return TRUE if fixed.
		is_fixed_sample_size = function(){
			private$fixed_sample
		},

		#' @description
		#' Asserts if all subjects arrived.
		assert_all_subjects_arrived = function(){
			if (should_run_asserts()) {
				if (private$fixed_sample & private$t < private$n){
					stop("This experiment is incomplete as all n subjects haven't arrived yet.")
				}
			}
		},

		#' @description
		#' Asserts if all responses are recorded.
		assert_all_responses_recorded = function(){
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
				if (sum(!is.na(private$y)) != length(private$w)){
					stop("This experiment is incomplete as all responses aren't recorded yet.")
				}
			}
		},

		#' @description
		#' Checks if the experiment is completed.
		#'
		#' @return	\code{TRUE} if experiment is complete, \code{FALSE} otherwise.
		check_experiment_completed = function(){
			if (private$fixed_sample & private$t < private$n){
				FALSE
			} else if (sum(!is.na(private$y)) != length(private$w)){
				FALSE
			} else {
				TRUE
			}
		},

		#' @description
		#' Checks if the experiment has a 50-50 allocation.
		assert_even_allocation = function(){
			if (should_run_asserts()) {
				if (private$prob_T != 0.5){
					stop("This type of design currently only works with even treatment allocation, i.e. you must set prob_T = 0.5 upon initialization")
				}
			}
		},

		#' @description
		#' Checks whether this design is a blocking design.
		assert_blocking_design = function(){
			if (should_run_asserts()) {
				if (!private$design_is_supported_blocking(self)){
					stop("This design requires a blocking design.")
				}
			}
		},

		#' @description
		#' Checks whether all blocks have the same number of subjects.
		#'
		#' @return `TRUE` invisibly if the block sizes are equal.
		assert_equal_block_sizes = function(){
			self$assert_blocking_design()
			block_ids = self$get_block_ids()
			block_sizes = as.integer(table(block_ids))
			if (should_run_asserts()) {
				if (length(block_sizes) > 1L && any(block_sizes != block_sizes[1L])) {
					stop("All blocks must have the same number of subjects.")
				}
			}
			invisible(TRUE)
		},

		#' @description
		#' Checks if the experiment has a fixed sample size.
		assert_fixed_sample = function(){
			if (should_run_asserts()) {
				if (!private$fixed_sample){
					stop("This type of design currently only works with fixed sample, i.e., you must specify n upon initialization")
				}
			}
		},

		#' @description
		#' Checks if the experiment has any censored responses
		#'
		#' @return	\code{TRUE} if any censored.
		any_censoring = function(){
			sum(private$dead) < length(private$dead)
		},

		#' @description Get t
		#'
		#' @return 			The current number of subjects.
		get_t = function(){
			private$t
		},

		#' @description Get raw X information
		#'
		#' @return 			A data frame of subject data.
		get_X_raw = function(){
			private$Xraw
		},

		#' @description Get imputed X information
		#'
		#' @return 		Same as \code{Xraw} except with imputations.
		get_X_imp = function(){
			private$Ximp
		},

		#' @description Get X matrix
		#'
		#' @return 			A numeric matrix of subject data.
		get_X = function(){
			private$X
		},

		#' @description Get y
		#'
		#' @return 			A numeric vector of subject responses.
		get_y = function(){
			private$y
		},

		#' @description Get w
		#'
		#' @return 			A binary vector of subject assignments.
		get_w = function(){
			private$w
		},

		#' @description If the design is a block design, get block identifiers (otherwise halts)
		#'
		#' @return 			An integer vector of block identifiers.
		get_block_ids = function(){
			if (!is.null(private$m) && length(private$m) == length(private$y)) {
				return(private$m)
			}
			block_ids = private$m
			strata_cols = private$strata_cols
			Xraw = private$Xraw
			if (is.null(block_ids) && is(self, "FixedDesignOptimalBlocks")) {
				block_ids = private$get_or_compute_block_ids()
			}
			if (is.null(block_ids) && !is.null(strata_cols) && length(strata_cols) > 0L &&
					nrow(Xraw) == length(private$y)) {
				strata_keys = private$get_strata_keys()
				if (length(strata_keys) == length(private$y)) {
					block_ids = match(strata_keys, unique(strata_keys))
				}
			}
			if (should_run_asserts()) {
				if (is.null(block_ids)) {
					stop("Block identifiers are undefined for this design.")
				}
			}
			block_ids = as.integer(block_ids)
			if (should_run_asserts()) {
				if (length(block_ids) != length(private$y)) {
					stop("Block identifiers are improperly sized for this design.")
				}
			}
			# Cache it back to private$m for next time if it wasn't already there
			private$m = block_ids
			block_ids
		},

		#' @description
		#' Summarize the covariates within blocks.
		#'
		#' @param block_ids A vector of block identifiers to summarize. If NULL, defaults to all blocks.
		#'
		#' @return A list of data.tables, one for each block.
		summarize_blocks = function(block_ids = NULL) {
			if (should_run_asserts()) {
				self$assert_blocking_design()
			}
			all_block_ids = self$get_block_ids()
			if (is.null(block_ids)) {
				block_ids = sort(unique(all_block_ids))
			}
			
			Xraw = self$get_X_raw()
			res = list()
			
			for (bid in block_ids) {
				idx = which(all_block_ids == bid)
				if (length(idx) == 0) {
					warning(paste("Block ID", bid, "not found in the design."))
					next
				}
				X_block = Xraw[idx, ]
				
				# Summary for each column
				col_summaries = list()
				for (col_name in names(X_block)) {
					col_data = X_block[[col_name]]
					if (is.numeric(col_data)) {
						col_summaries[[col_name]] = data.table::data.table(
							variable = col_name,
							level = NA_character_,
							mean_or_pct = mean(col_data, na.rm = TRUE),
							sd = stats::sd(col_data, na.rm = TRUE)
						)
					} else {
						# Categorical (factor or character)
						counts = table(col_data, useNA = "no")
						pcts = prop.table(counts) * 100
						if (length(pcts) == 0) {
							col_summaries[[col_name]] = data.table::data.table(
								variable = col_name,
								level = "N/A",
								mean_or_pct = NA_real_,
								sd = NA_real_
							)
						} else {
							col_summaries[[col_name]] = data.table::data.table(
								variable = col_name,
								level = names(pcts),
								mean_or_pct = as.numeric(pcts),
								sd = NA_real_
							)
						}
					}
				}
				res[[as.character(bid)]] = data.table::rbindlist(col_summaries)
			}
			res
		},

		#' @description Get n, the sample size
		#'
		#' @return 			The number of subjects.
		get_n = function(){
			ifelse(private$fixed_sample, private$n, private$t)
		},

		#' @description Get dead
		#'
		#' @return 			A binary vector of whether the subject is dead.
		get_dead = function(){
			private$dead
		},

		#' @description Get probability of treatment
		#'
		#' @return 			The specified probability.
		get_prob_T = function(){
			private$prob_T
		},

		#' @description
		#' Checks if the design contains any covariates.
		#'
		#' @return \code{TRUE} if there are covariates, \code{FALSE} otherwise.
		has_covariates = function(){
			!is.null(private$Xraw) && ncol(private$Xraw) > 0L
		},

		#' @description Get response type
		#'
		#' @return 			The specified response type.
		get_response_type = function(){
			private$response_type
		},

		#' @description Get the missingness method
		#'
		#' @return 			The missingness handling method: \code{"impute"}, \code{"drop_column"},
		#'   or \code{"error"}.
		get_missingness_method = function(){
			private$missingness_method
		},

		#' @description
		#' Duplicate this design object
		#'
		#' @param verbose 	A flag for verbosity.
		#' @return 			A new `Design` object with the same data
		duplicate = function(verbose = FALSE){
			if (should_run_asserts()) {
				self$assert_all_responses_recorded() #can't duplicate without the experiment being done
			}
			# Use the built-in R6 clone method (shallow by default) to bypass $new() logic.
			d = self$clone()
			d$.__enclos_env__$private$verbose = verbose
			d
		}

	),

	active = list(
		#' @field num_cores Current number of cores in the global budget.
		num_cores = function() get_num_cores()
	),

	private = list(
		all_subject_data_cache = list(),
		t = 0L,
		n = NULL,
		Xraw = data.table(),
		p_raw_t = NULL,
		Ximp = data.table(),
		X = NULL,
		w = numeric(),
		y = numeric(),
		dead = numeric(),
		m = NULL,
		kk_xm_structural     = NULL,
		kk_xm_m_vec          = NULL,
		kk_lin_xm_structural = NULL,
		kk_lin_xm_m_vec      = NULL,
		kk_cluster_id        = NULL,
		kk_cluster_id_m_vec  = NULL,
		permutations_cache   = list(),
		lin_centered_covariates = NULL,
		kk_boot_pair_rows    = NULL,
		kk_boot_i_reservoir  = NULL,
		kk_boot_n_reservoir  = NULL,
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			list(i_b = sample.int(private$n, private$n, replace = TRUE), m_vec_b = NULL)
		},
		strata_cols = NULL,
		prob_T = NULL,
		response_type = NULL,
		fixed_sample = NULL,
		include_is_missing_as_a_new_feature = NULL,
		missingness_method = "impute",
		verbose = NULL,
		design_is_supported_blocking = function(des_obj){
			is(des_obj, "DesignSeqOneByOneKK14") ||
				is(des_obj, "FixedDesigniBCRD") ||
				is(des_obj, "FixedDesignBinaryMatch") ||
				is(des_obj, "FixedDesignBlocking") ||
				is(des_obj, "FixedDesignOptimalBlocks") ||
				is(des_obj, "DesignSeqOneByOneSPBR") ||
				is(des_obj, "DesignSeqOneByOneRandomBlockSize")
		},
		y_i_t_i = list(),	 #at what point during the experiment are the subjects recorded?
		uses_covariates = FALSE, #does this design use the covariates to make assignments? The default is FALSE
		resample_assignment = function(){
			n = private$n
			i_b = sample(n, n, replace = TRUE)
			private$w    = private$w[i_b]
			private$y    = private$y[i_b]
			private$dead = private$dead[i_b]
			if (!is.null(private$m)){
				private$m = private$m[i_b]
			}
			invisible(self)
		},

		covariate_impute_if_necessary_and_then_create_model_matrix = function(){
			#make a copy... sometimes the raw will be the same as the imputed if there are no imputations
			private$Ximp = copy(private$Xraw)

			column_has_missingness = columns_have_missingness_cpp(private$Xraw)
			if (any(column_has_missingness)){
				if (private$missingness_method == "error"){
					if (should_run_asserts()) {
						missing_names = names(private$Xraw)[column_has_missingness]
						stop("Missing values detected in covariate(s): ",
							paste(missing_names, collapse = ", "),
							". Set missingness_method = \"impute\" or \"drop_column\" to handle missing data automatically.")
					}
				} else if (private$missingness_method == "drop_column"){
					cols_to_keep = which(!column_has_missingness)
					private$Ximp = private$Ximp[, ..cols_to_keep]
				} else {
					# "impute": random-forest imputation (missRanger with missForest fallback)

					#deal with include_is_missing_as_a_new_feature here
					if (private$include_is_missing_as_a_new_feature){
						missing_cols_idx = which(column_has_missingness)
						if (length(missing_cols_idx) > 0){
							# Use C++ function to create missingness indicators efficiently
							missingness_indicators = create_missingness_indicators_cpp(private$Ximp, missing_cols_idx)

							# Add the new columns to Ximp
							for (col_name in names(missingness_indicators)) {
								private$Ximp[[col_name]] = missingness_indicators[[col_name]]
							}
						}
					}

					#we need to convert characters into factor for the imputation to work
					col_types = get_column_types_cpp(private$Ximp)
					idx_cols_to_convert_to_factor = which(col_types == "character")
					private$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]

					#now do the imputation here by using missRanger (fast but fragile) and if that fails, use missForest (slow but more robust)
					private$Ximp = tryCatch({
											if (any(!is.na(private$y))){
												suppressWarnings(missRanger(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), verbose = FALSE, num.threads = self$num_cores)[, 1 : ncol(private$Ximp)])
											} else {
												suppressWarnings(missRanger(private$Ximp, verbose = FALSE, num.threads = self$num_cores))
											}
										}, error = function(e){
											if (any(!is.na(private$y))){
												suppressWarnings(missForest(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), num.threads = self$num_cores)$ximp[, 1 : ncol(private$Ximp)])
											} else {
												suppressWarnings(missForest(private$Ximp, num.threads = self$num_cores)$ximp)
											}
										}
									)
				}
			}

			analysis_col_names = names(private$Ximp)[!startsWith(names(private$Ximp), ".assignment_only_")]
			Ximp_for_model = if (length(analysis_col_names)) {
				private$Ximp[, ..analysis_col_names]
			} else {
				private$Ximp[, .SD, .SDcols = integer(0)]
			}

			#now let's drop any columns that don't have any variation
			num_unique_values_per_column = count_unique_values_cpp(Ximp_for_model)
			Ximp_for_model = Ximp_for_model[, .SD, .SDcols = which(num_unique_values_per_column > 1)]

			#for nonblank data frames...
			if (ncol(Ximp_for_model) > 0){
				#now we need to update the numeric model matrix which may have expanded due to new factors, new missingness cols, etc
				private$X = model.matrix(~ ., data = Ximp_for_model)[, -1, drop = FALSE]
				# Ensure it is a numeric matrix (not character)
				if (should_run_asserts()) {
					if (is.character(private$X)){
						stop("model.matrix returned a character matrix - this should not happen.")
					}
				}
				private$X = drop_linearly_dependent_cols(private$X)$M

				if (should_run_asserts()) {
					if (nrow(private$X) != nrow(private$Xraw) | nrow(private$X) != nrow(private$Ximp) | nrow(private$Ximp) != nrow(private$Xraw)){
						stop("improper sizing for the internal X representation")
					}
				}
			} else {
				#blank covariate matrix
				private$X = matrix(NA, nrow = nrow(private$Xraw), ncol = 0)
			}

		},

		compute_all_subject_data = function(){
			i_present_y = which(!is.na(private$y))
			i_all = 1 : private$t
			i_all_y_present = intersect(i_all, i_present_y)
			
			# Cache lookup
			# Since covariates are fixed and NA positions in y are fixed during randomization,
			# the set of subjects with responses up to t is constant for a given t.
			cache_key = as.character(private$t)
			if (!is.null(private$all_subject_data_cache[[cache_key]])) {
				cpp_result = private$all_subject_data_cache[[cache_key]]
			} else {
				# Call consolidated C++ function for all matrix computations
				cpp_result = compute_all_subject_data_cpp(
					as.matrix(private$X[1:private$t, , drop = FALSE]),
					private$t,
					as.integer(i_all_y_present)
				)
				# Restore column names
				X_names = colnames(private$X)

				if (length(cpp_result$cols_prev) > 0) {
					nms = X_names[cpp_result$cols_prev]
					colnames(cpp_result$X_prev) = nms
					names(cpp_result$xt_prev) = nms
				}

				if (length(cpp_result$cols_all) > 0) {
					colnames(cpp_result$X_all) = X_names[cpp_result$cols_all]
				}

				if (length(cpp_result$cols_all_scaled) > 0) {
					nms = X_names[cpp_result$cols_all_scaled]
					colnames(cpp_result$X_all_scaled) = nms
					names(cpp_result$xt_all_scaled) = nms
				}

				if (length(cpp_result$cols_all_with_y_scaled) > 0) {
					colnames(cpp_result$X_all_with_y_scaled) = X_names[cpp_result$cols_all_with_y_scaled]
				}
				
				if (is.null(private$all_subject_data_cache)) private$all_subject_data_cache = list()
				private$all_subject_data_cache[[cache_key]] = cpp_result
			}

			# Add the simple array slices that don't need C++ optimization
			# These MUST NOT be cached because w and y change during randomization!
			cpp_result$w_all_with_y_scaled = private$w[i_all_y_present]
			cpp_result$y_all = private$y[i_all_y_present]
			cpp_result$dead_all = private$dead[i_all_y_present]

			cpp_result
		},

		assign_wt_Bernoulli = function(){
			rbinom(1, 1, private$prob_T)
		},

		has_private_method = function(method_name) {
			exists(method_name, envir = self$.__enclos_env__$private, inherits = FALSE)
		}
	)
)
