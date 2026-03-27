#' An incomplete / balanced completely randomized Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments.
#'
#' @export
DesignSeqOneByOneiBCRD = R6::R6Class("DesignSeqOneByOneiBCRD",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description
		#' Initialize a balanced sequential experimental design
		#'
		#' @param	response_type 	"continuous", "incidence", "proportion", "count", "survival", or "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `DesignSeqOneByOneiBCRD` object
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
			nT = sum(private$w == 1, na.rm = TRUE)
			nC = sum(private$w == 0, na.rm = TRUE)
			
			if (is.null(private$n)){
				#if n is not fixed, we cannot really ensure balance at the end, 
				#so we just use Bernoulli
				private$assign_wt_Bernoulli()
			} else {
				#if n is fixed, we use the remaining slots
				nT_rem = round(private$n * private$prob_T) - nT
				nC_rem = (private$n - round(private$n * private$prob_T)) - nC
				
				if (nT_rem <= 0) return(0)
				if (nC_rem <= 0) return(1)
				
				rbinom(1, 1, nT_rem / (nT_rem + nC_rem))
			}
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to balanced randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			generate_permutations_ibcrd_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		}
	),
	private = list(
		redraw_w_according_to_design = function(){
			n_T_total = round(private$t * private$prob_T) 
			private$w[1:private$t] = shuffle_cpp(c(rep(1, n_T_total), rep(0, private$t - n_T_total)))
		}
	)
)
