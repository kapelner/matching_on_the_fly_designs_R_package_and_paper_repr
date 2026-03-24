# A Fixed Design
#
# @description
# An abstract R6 Class encapsulating the data and functionality for a fixed experimental design.
# This class takes care of whole-experiment randomization.
#
# @keywords internal
#' @export
FixedDesign = R6::R6Class("FixedDesign",
	inherit = Design,
	public = list(
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

		randomize = function(){
			private$draw_one_w()
		},

		supports_resampling = function(){
			class(self)[1] != "FixedDesign"
		},

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
