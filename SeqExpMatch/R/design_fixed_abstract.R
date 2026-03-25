#' A Fixed Design
#'
#' @description
#' An abstract R6 Class encapsulating the data and functionality for a fixed experimental design.
#' This class takes care of whole-experiment randomization.
#'
#' @keywords internal
#' @export
FixedDesign = R6::R6Class("FixedDesign",
	inherit = Design,
	public = list(
		#' @description
		#' Initialize a fixed experimental design
		#'
		#' @param	response_type 	"continuous", "incidence", "proportion", "count", "survival", or "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `FixedDesign` object
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
		#' Randomize treatment assignments for the entire experiment.
		randomize = function(){
			private$draw_one_w()
		},

		#' @description
		#' Check if the design supports resampling.
		#'
		#' @return 	TRUE if supported.
		supports_resampling = function(){
			class(self)[1] != "FixedDesign"
		},

		#' @description
		#' Redraw treatment assignments according to the design.
		redraw_w_according_to_design = function(){
			stop("Must be implemented by subclass.")
		},

		#' @description
		#' Draw multiple treatment assignment vectors.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			stop("Must be implemented by subclass.")
		}
	),

	private = list(
		draw_one_w = function() {
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		}
	)
)
