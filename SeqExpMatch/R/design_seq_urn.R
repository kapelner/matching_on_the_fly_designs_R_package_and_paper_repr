#' Wei's Urn Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for Wei's Urn sequential experimental design.
#' This design uses an adaptive urn model where the probability of assignment to a group
#' decreases as the number of subjects in that group increases.
#'
#' @export
SeqDesignUrn = R6::R6Class("SeqDesignUrn",
	inherit = SeqDesign,
	public = list(
		#'
		#' @description
		#' Initialize an Urn sequential experimental design
		#'
		#' @param alpha The initial number of balls of each type (Treatment/Control) in the urn.
		#' @param beta The number of balls of the opposite type to add to the urn after an assignment.
		#' @param	response_type 	The data type of response values.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag for verbosity.
		#' @return	A new `SeqDesignUrn` object
		#'
		initialize = function(
						alpha = 1,
						beta = 1,
						response_type = "continuous",
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_cores = 1,
						verbose = FALSE
					) {
			assertNumber(alpha, lower = 0)
			assertNumber(beta, lower = 0)
			
			super$initialize(response_type, 0.5, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			
			private$alpha = alpha
			private$beta = beta
		},

		assign_wt = function(){
			# Probability of Treatment based on current counts
			# In Wei's Urn, P(T) = (alpha + beta * nC) / (2 * alpha + beta * (nT + nC))
			# where nT and nC are current counts of assigned subjects.
			nT = sum(private$w == 1, na.rm = TRUE)
			nC = sum(private$w == 0, na.rm = TRUE)
			
			prob_T = (private$alpha + private$beta * nC) / (2 * private$alpha + private$beta * (nT + nC))
			
			rbinom(1, 1, prob_T)
		}
	),
	private = list(
		alpha = NULL,
		beta = NULL,

		redraw_w_according_to_design = function(){
			# Sequential simulation of the urn process
			new_w = rep(NA_real_, private$t)
			nT = 0
			nC = 0
			for (i in 1:private$t) {
				prob_T = (private$alpha + private$beta * nC) / (2 * private$alpha + private$beta * (nT + nC))
				w_i = rbinom(1, 1, prob_T)
				new_w[i] = w_i
				if (w_i == 1) nT = nT + 1 else nC = nC + 1
			}
			private$w[1:private$t] = new_w
		}
	)
)
