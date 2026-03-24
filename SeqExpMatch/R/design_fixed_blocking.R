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
		},

		redraw_w_according_to_design = function(){
			strata_keys = vapply(1:private$t, function(i) {
				vals = vapply(private$strata_cols, function(col) {
					val = private$Xraw[i, ][[col]]
					if (is.na(val)) "NA" else as.character(val)
				}, character(1))
				paste(vals, collapse = "|")
			}, character(1))
			
			unique_keys = unique(strata_keys)
			new_w = rep(NA_real_, private$t)
			
			for (key in unique_keys) {
				idxs = which(strata_keys == key)
				m = length(idxs)
				n_T = round(m * private$prob_T)
				new_w[idxs] = sample(c(rep(1, n_T), rep(0, m - n_T)))
			}
			private$w[1:private$t] = new_w
		}
	),
	private = list(
		strata_cols = NULL
	)
)
