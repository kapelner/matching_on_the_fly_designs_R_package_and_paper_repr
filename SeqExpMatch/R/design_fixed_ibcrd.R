#' A balanced completely randomized Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed balanced completely randomized experimental design.
#'
#' @export
FixedDesigniBCRD = R6::R6Class("FixedDesigniBCRD",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a fixed balanced completely randomized experimental design
		#'
		#' @param	response_type 	"continuous", "incidence", "proportion", "count", "survival", or "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `FixedDesigniBCRD` object
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
		#' Draw multiple treatment assignment vectors according to balanced randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r){
			generate_permutations_ibcrd_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		},

		#' @description
		#' Redraw treatment assignments according to the balanced completely randomized design.
		redraw_w_according_to_design = function(){
			n_T_total = round(private$t * private$prob_T)
			if (abs(n_T_total - private$t * private$prob_T) > 1e-10) {
				warning("prob_T * t is not an integer; rounding to the nearest integer for balanced allocation.")
			}
			private$w[1:private$t] = shuffle_cpp(c(rep(1, n_T_total), rep(0, private$t - n_T_total)))
		}
	),
	private = list(
	)
)
