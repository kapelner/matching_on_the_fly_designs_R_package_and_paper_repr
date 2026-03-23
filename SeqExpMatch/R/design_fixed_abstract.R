# A Fixed Design
#
# @description
# An abstract R6 Class encapsulating the data and functionality for a fixed experimental design.
# This class takes care of data storage, response handling, and whole-experiment randomization.
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
			self$redraw_w_according_to_design()
		},

		# @description
		# For those who wish to use this package for analysis on already-completed experimental data
		#
		# @param w 		The binary responses as a numeric vector of length equal to the number of subjects in the study
		#
		add_all_subject_assignments = function(w) {
			assertIntegerish(w, lower = 0, upper = 1, any.missing = FALSE, len = private$t)
			private$w = w
		},

		redraw_w_according_to_design = function(){
			stop("Must be implemented by subclass.")
		},

		draw_ws_according_to_design = function(r = 100){
			# Abstract for multiple draws
			# Default implementation loops redraw_w_according_to_design
			# Most designs override this with a fast C++ path
			w_mat = matrix(NA_real_, nrow = self$get_n(), ncol = r)
			for (j in 1 : r){
				self$redraw_w_according_to_design()
				w_mat[, j] = self$get_w()
			}
			w_mat
		}
	)
)
