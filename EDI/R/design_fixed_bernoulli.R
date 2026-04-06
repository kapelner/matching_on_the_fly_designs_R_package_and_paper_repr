#' A completely randomized / Bernoulli Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed Bernoulli experimental design.
#'
#' @export
FixedDesignBernoulli = R6::R6Class("FixedDesignBernoulli",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a fixed Bernoulli experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `FixedDesignBernoulli` object
		initialize = function(
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						
						verbose = FALSE
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to Bernoulli randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r){
			generate_permutations_bernoulli_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		}
	),
	private = list(
	)
)
