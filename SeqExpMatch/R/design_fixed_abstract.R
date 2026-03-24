# A Fixed Design
#
# @description
# A concrete R6 Class encapsulating the data and functionality for a completed fixed experimental design.
# This class can store completed-study covariates, assignments, and responses for downstream analysis.
# Redraw-based randomization is only available in concrete subclasses that implement an assignment mechanism.
#
# @keywords internal
#' @export
FixedDesign = R6::R6Class("FixedDesign",
	inherit = Design,
	public = list(
		#' @description
		#' Initialize a completed fixed experimental design container.
		#'
		#' Unlike \code{Design}, this class requires a fixed sample size \code{n}.
		#' The base \code{FixedDesign} class is intended for representing already-completed
		#' experiments for analysis. It does not itself define how assignments are redrawn.
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n,
				num_cores = 1,
				verbose = FALSE
			) {
			assertCount(n, positive = TRUE)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
		},

		# @description Whether this design supports redraw-based resampling
		#
		# @return 			\code{TRUE} for concrete subclasses with an assignment mechanism and \code{FALSE} for plain \code{FixedDesign}.
		supports_resampling = function(){
			class(self)[1] != "FixedDesign"
		},

		randomize = function(){
			if (!self$supports_resampling()){
				stop("Plain FixedDesign objects do not support randomization because the base class does not define an assignment redraw mechanism.")
			}
			self$redraw_w_according_to_design()
		},

		# @description
		# For those who wish to use this package for analysis on already-completed experimental data
		#
		# @param m 		The integer vector of matched pair IDs of length equal to the number of subjects in the study
		#
		add_all_subject_matched_pair_ids = function(m) {
			assertIntegerish(m, any.missing = FALSE, len = private$t)
			private$m = m
		},

		redraw_w_according_to_design = function(){
			stop("Plain FixedDesign objects cannot redraw treatment assignments. Use a concrete design subclass for randomization-based procedures.")
		},

		draw_ws_according_to_design = function(r = 100){
			if (!self$supports_resampling()){
				stop("Plain FixedDesign objects cannot draw randomization replicates because the base class does not define an assignment mechanism.")
			}
			# Default implementation loops redraw_w_according_to_design.
			# Most concrete designs override this with a faster path.
			w_mat = matrix(NA_real_, nrow = private$n, ncol = r)
			for (j in 1 : r){
				self$redraw_w_according_to_design()
				w_mat[, j] = self$get_w()
			}
			w_mat
		}
	)
)
