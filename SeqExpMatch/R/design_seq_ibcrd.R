#' An incomplete / balanced completely randomized Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#'
#' @export
SeqDesigniBCRD = R6::R6Class("SeqDesigniBCRD",
	inherit = SeqDesign,
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
			self$assert_fixed_sample()
		},

		assign_wt = function(){
			n_T_total = round(private$n * private$prob_T) 
			nT = sum(private$w == 1, na.rm = TRUE)
			nC = sum(private$w == 0, na.rm = TRUE)
			sample(c(rep(1, n_T_total - nT), rep(0, private$n - n_T_total - nC)), 1)
		}
	),
	private = list(
		redraw_w_according_to_design = function(){
			n_T_total = round(private$t * private$prob_T) 
			private$w[1:private$t] = shuffle_cpp(c(rep(1, n_T_total), rep(0, private$t - n_T_total)))
		}
	)
)
