#' A stratified blocking Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed stratified blocking experimental design.
#'
#' @export
FixedDesignBlocking = R6::R6Class("FixedDesignBlocking",
	inherit = FixedDesign,
	public = list(
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
			
			unique_keys = unique(strata_keys)
			strata_indices = lapply(unique_keys, function(key) which(strata_keys == key))
			
			generate_permutations_blocking_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T),
				strata_indices
			)$w_mat
		}
	),
	private = list(
		strata_cols = NULL
	)
)
