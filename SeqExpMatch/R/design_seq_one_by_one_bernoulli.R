#' A completely randomized / Bernoulli Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#'
#' @export
DesignSeqOneByOneBernoulli = R6::R6Class("DesignSeqOneByOneBernoulli",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description
		#' Initialize a Bernoulli sequential experimental design
		#'
		#' @param	response_type 	The data type of response values which must be one of the following:
		#' 								"continuous" (the default),
		#' 								"incidence",
		#' 								"proportion",
		#' 								"count",
		#' 								"survival",
		#' 								"ordinal".
		#' @param	prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		#' @param include_is_missing_as_a_new_feature     If missing data is present in a variable,
		#'   should we include another dummy variable for its missingness? The default is \code{TRUE}.
		#' @param	n			The sample size (if fixed). Default is \code{NULL}.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @return	A new `DesignSeqOneByOneBernoulli` object
		#'
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
		#' Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			rbinom(1, 1, private$prob_T)
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to Bernoulli randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			generate_permutations_bernoulli_cpp(as.integer(private$t), as.integer(r), as.numeric(private$prob_T))$w_mat
		}
	),
	private = list(
	)
)
