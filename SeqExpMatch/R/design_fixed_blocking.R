#' A stratified blocking Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed stratified blocking experimental design.
#' This design ensures balance within specified strata using the \pkg{randomizr} package.
#'
#' @export
FixedDesignBlocking = R6::R6Class("FixedDesignBlocking",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a fixed stratified blocking experimental design
		#'
		#' @param strata_cols A character vector of column names to use for stratification.
		#' @param	response_type 	The data type of response values.
		#' @param	prob_T	The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores.
		#' @param verbose A flag for verbosity.
		#' @return	A new `FixedDesignBlocking` object
		#'
		initialize = function(
						strata_cols,
						response_type = "continuous",
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_cores = 1,
						verbose = FALSE
					) {
			assertCharacter(strata_cols, min.len = 1)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$strata_cols = strata_cols
			private$uses_covariates = TRUE
		},

		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			
			strata_keys = vapply(1:private$t, function(i) {
				vals = vapply(private$strata_cols, function(col) {
					val = private$Xraw[i, ][[col]]
					if (is.na(val)) "NA" else as.character(val)
				}, character(1))
				paste(vals, collapse = "|")
			}, character(1))
			
			# Use randomizr::block_ra for canonical stratified blocking
			# We want numeric 0/1. block_ra returns a factor/numeric vector.
			w_mat = replicate(r, as.numeric(as.character(randomizr::block_ra(blocks = strata_keys, prob = private$prob_T))))
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),
	private = list(
		strata_cols = NULL
	)
)
