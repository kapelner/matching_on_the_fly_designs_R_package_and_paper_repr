#' A stratified blocking Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed stratified blocking experimental design.
#'
#' @export
FixedDesignBlocking = R6::R6Class("FixedDesignBlocking",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a fixed stratified blocking experimental design
		#'
		#' @param	response_type 	"continuous", "incidence", "proportion", "count", "survival", or "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `FixedDesignBlocking` object
		initialize = function(
						response_type = "continuous",
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_cores = 1,
						verbose = FALSE
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to stratified blocking.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r){
			# Logic for stratified blocking permutations
			# Get strata from Ximp
			strata = apply(private$Ximp, 1, paste, collapse = "|")
			strata_levels = unique(strata)
			strata_indices = lapply(strata_levels, function(lev) which(strata == lev))
			
			generate_permutations_blocking_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T),
				strata_indices
			)$w_mat
		}
	),
	private = list(
		redraw_w_according_to_design = function(){
			# Manual redraw logic
			strata = apply(private$Ximp, 1, paste, collapse = "|")
			strata_levels = unique(strata)
			new_w = rep(NA_real_, private$t)
			
			for (lev in strata_levels) {
				idxs = which(strata == lev)
				m = length(idxs)
				n_T = round(m * private$prob_T)
				# Balanced shuffle within stratum
				new_w[idxs] = sample(c(rep(1, n_T), rep(0, m - n_T)))
			}
			
			private$w[1:private$t] = new_w
		}
	)
)
