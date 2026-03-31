#' Efron's (1971) Biased Coin Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for Efron's biased coin
#' sequential experimental design.
#'
#' @export
DesignSeqOneByOneEfron = R6::R6Class("DesignSeqOneByOneEfron",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description
		#' Initialize an Efron biased coin sequential experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param verbose A flag for verbosity.
		#' @param weighted_coin_prob The probability of assigning to the under-represented group.
		#'
		#' @return	A new `DesignSeqOneByOneEfron` object
		initialize = function(
						response_type = "continuous",
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						
						verbose = FALSE,
						weighted_coin_prob = 2/3
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$weighted_coin_prob = weighted_coin_prob
		},

		#' @description
		#' Assign the next subject to a treatment group
		#'
		#' @return 	The treatment assignment (0 or 1)
		assign_wt = function(){
			#if it's the first subject or if balance is equal, then Bernoulli
			nT = sum(private$w == 1, na.rm = TRUE)
			nC = sum(private$w == 0, na.rm = TRUE)
			
			if (nT == nC){
				private$assign_wt_Bernoulli()
			} else {
				#assign to the group with fewer subjects with probability weighted_coin_prob
				if (nT > nC){
					rbinom(1, 1, 1 - private$weighted_coin_prob)
				} else {
					rbinom(1, 1, private$weighted_coin_prob)
				}
			}
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to Efron's design.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			generate_permutations_efron_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T),
				as.numeric(private$weighted_coin_prob)
			)$w_mat
		}
	),
	private = list(
		weighted_coin_prob = NULL

	)
)
