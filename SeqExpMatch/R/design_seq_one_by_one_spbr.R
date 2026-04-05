#' A stratified permuted block Sequential Design (SPBR)
#'
#' An R6 Class encapsulating the data and functionality for a stratified permuted
#' block sequential experimental design.
#' This design ensures balance within specified strata using blocks of a fixed size.
#'
#' @export
DesignSeqOneByOneSPBR = R6::R6Class("DesignSeqOneByOneSPBR",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description
		#' Initialize a stratified permuted block sequential experimental design
		#'
		#' @param strata_cols A character vector of column names to use for stratification.
		#' @param block_size The size of the permuted blocks.
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `DesignSeqOneByOneSPBR` object
		initialize = function(
						strata_cols,
						block_size = 4,
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						
						verbose = FALSE
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$strata_cols = strata_cols
			private$block_size = as.integer(block_size)
			private$uses_covariates = TRUE
			private$strata_states = new.env(parent = emptyenv())
			
			if (abs(block_size * prob_T - round(block_size * prob_T)) > 1e-10) {
				stop("block_size must result in an integer number of treatment assignments (block_size * prob_T).")
			}
		},

		#' @description
		#' Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			x_new = private$Xraw[private$t, ]
			key = private$get_strata_key(x_row = x_new)
			
			if (is.null(private$strata_states[[key]]) || length(private$strata_states[[key]]) == 0) {
				n_T = round(private$block_size * private$prob_T)
				n_C = private$block_size - n_T
				new_block = sample(c(rep(1, n_T), rep(0, n_C)))
				private$strata_states[[key]] = new_block
			}
			
			block = private$strata_states[[key]]
			w_t = block[1]
			private$strata_states[[key]] = block[-1]
			w_t
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to SPBR design.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			strata_keys = vapply(1:private$t, function(i) {
				private$get_strata_key(private$Xraw[i, ])
			}, character(1))
			
			generate_permutations_spbr_cpp(
				as.character(unname(strata_keys)),
				as.integer(private$block_size),
				as.numeric(private$prob_T),
				as.integer(r)
			)$w_mat
		}
	),
	private = list(
		strata_cols = NULL,
		block_size = NULL,
		strata_states = NULL,

		get_strata_key = function(x_row) {
			vals = vapply(private$strata_cols, function(col) {
				val = x_row[[col]]
				if (is.na(val)) "NA" else as.character(val)
			}, character(1))
			paste(vals, collapse = "|")
		}

	)
)
