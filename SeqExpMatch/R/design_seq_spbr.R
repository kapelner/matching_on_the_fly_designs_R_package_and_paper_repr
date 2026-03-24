#' A stratified permuted block Sequential Design (SPBR)
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a stratified permuted block sequential experimental design.
#'
#' @export
SeqDesignSPBR = R6::R6Class("SeqDesignSPBR",
	inherit = SeqDesign,
	public = list(
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
		},

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
		},

		redraw_w_according_to_design = function(){
			strata_keys = vapply(1:private$t, function(i) {
				private$get_strata_key(private$Xraw[i, ])
			}, character(1))
			private$w[1:private$t] = spbr_redraw_w_cpp(as.character(unname(strata_keys)), as.integer(private$block_size), as.numeric(private$prob_T))
		},

		get_bootstrap_indices = function() {
			strata_keys = vapply(1:private$t, function(i) {
				private$get_strata_key(private$Xraw[i, ])
			}, character(1))
			stratified_bootstrap_indices_cpp(as.character(unname(strata_keys)))
		}
	)
)
