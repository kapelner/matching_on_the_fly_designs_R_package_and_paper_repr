#' Atkinson's (1982) Covariate-Adjusted Biased Coin Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments.
#'
#' @export
DesignSeqOneByOneAtkinson = R6::R6Class("DesignSeqOneByOneAtkinson",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description
		#' Initialize an Atkinson sequential experimental design
		#'
		#' @param	response_type 	The data type of response values.
		#' @param	prob_T	The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag for verbosity.
		#'
		#' @return 			A new `DesignSeqOneByOneAtkinson` object
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
			#if it's too early in the trial or if all the assignments are the same, then randomize
			if (private$t <= ncol(private$Xraw) + 2 + 1){
				rbinom(1, 1, private$prob_T)
			} else {
				all_subject_data = private$compute_all_subject_data()
				tryCatch({
					atkinson_assign_weight_cpp(
						private$w[1 : (private$t - 1)],
						all_subject_data$X_prev,
						all_subject_data$xt_prev,
						all_subject_data$rank_prev,
						private$t
					)
				}, error = function(e){
					rbinom(1, 1, private$prob_T)
				})
			}
		},

		#' @description
		#' Draw multiple treatment assignment vectors.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			generate_permutations_atkinson_cpp(
				as.matrix(private$X[1:private$t, , drop = FALSE]),
				as.integer(private$t),
				as.integer(ncol(private$Xraw)),
				as.numeric(private$prob_T),
				as.integer(r)
			)$w_mat
		}

	),
	private = list(
	)
)
