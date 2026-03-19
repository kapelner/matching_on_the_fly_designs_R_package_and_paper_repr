#' A stratified permuted block Sequential Design (SPBR)
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a stratified permuted block sequential experimental design.
#' This design ensures balance within specified strata using permuted blocks.
#'
#' @export
SeqDesignSPBR = R6::R6Class("SeqDesignSPBR",
	inherit = SeqDesign,
	public = list(
		#'
		#' @description
		#' Initialize a stratified permuted block sequential experimental design
		#'
		#' @param strata_cols A character vector of column names to use for stratification.
		#' @param block_size The size of the permuted blocks. Must be a multiple of the inverse of \code{prob_T} to ensure integer treatment/control counts.
		#' @param	response_type 	The data type of response values which must be one of the following:
		#' 								"continuous" (the default),
		#' 								"incidence",
		#' 								"proportion",
		#' 								"count",
		#' 								"survival",
		#' 								"ordinal".
		#' @param	prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		#' @param include_is_missing_as_a_new_feature     If missing data is present in a variable,
		#'   should we include another dummy variable for its missingness? Default is \code{TRUE}.
		#' @param	n			The sample size (if fixed). Default is \code{NULL} for not fixed.
		#' @param num_cores The number of CPU cores to use to parallelize the sampling.
		#' @param verbose A flag indicating whether messages should be displayed. Default is \code{FALSE}.
		#' @return	A new `SeqDesignSPBR` object
		#'
		initialize = function(
						strata_cols,
						block_size = 4,
						response_type = "continuous",
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_cores = 1,
						verbose = FALSE
					) {
			assertCharacter(strata_cols, min.len = 1)
			assertCount(block_size, positive = TRUE)
			
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			
			private$strata_cols = strata_cols
			private$block_size = as.integer(block_size)
			private$uses_covariates = TRUE
			private$strata_states = new.env(parent = emptyenv())
			
			# Validation for block size and prob_T
			if (abs(block_size * prob_T - round(block_size * prob_T)) > 1e-10) {
				stop("block_size * prob_T must be an integer.")
			}
		}
	),
	private = list(
		strata_cols = NULL,
		block_size = NULL,
		strata_states = NULL, # hash map of stratum -> vector of remaining assignments

		assign_wt = function(){
			cat("      assign_wt start for t =", private$t, "\n")
			x_new = private$Xraw[private$t, ]
			key = private$get_strata_key(x_row = x_new)
			cat("      strata key =", key, "\n")
			
			if (is.null(private$strata_states[[key]]) || length(private$strata_states[[key]]) == 0) {
				cat("      refilling block for key =", key, "\n")
				# Refill block
				n_T = round(private$block_size * private$prob_T)
				n_C = private$block_size - n_T
				new_block = sample(c(rep(1, n_T), rep(0, n_C)))
				private$strata_states[[key]] = new_block
			}
			
			# Pop one
			block = private$strata_states[[key]]
			w_t = block[1]
			private$strata_states[[key]] = block[-1]
			cat("      assigned w_t =", w_t, "\n")
			w_t
		},

		get_strata_key = function(x_row) {
			# Concatenate strata column values into a key string
			vals = vapply(private$strata_cols, function(col) {
				val = x_row[[col]]
				if (is.na(val)) "NA" else as.character(val)
			}, character(1))
			paste(vals, collapse = "|")
		},

		redraw_w_according_to_design = function(){
			# Re-simulating from scratch using C++ speedup
			# First, we need the keys for all subjects added so far
			strata_keys = vapply(1:private$t, function(i) {
				private$get_strata_key(private$Xraw[i, ])
			}, character(1))
			
			private$w[1:private$t] = spbr_redraw_w_cpp(as.character(unname(strata_keys)), as.integer(private$block_size), as.numeric(private$prob_T))
		},

		# Option B: Fast stratified bootstrap for SPBR
		get_bootstrap_indices = function() {
			# We need the strata keys for all subjects added so far (private$t)
			strata_keys = vapply(1:private$t, function(i) {
				private$get_strata_key(private$Xraw[i, ])
			}, character(1))
			
			stratified_bootstrap_indices_cpp(as.character(unname(strata_keys)))
		}
	)
)
