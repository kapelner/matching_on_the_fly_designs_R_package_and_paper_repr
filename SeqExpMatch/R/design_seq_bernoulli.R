#' A completely randomized / Bernoulli Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#'
#' @export
SeqDesignBernoulli = R6::R6Class("SeqDesignBernoulli",
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
		},

		assign_wt = function(){
			rbinom(1, 1, private$prob_T)
		},

		draw_ws_according_to_design = function(r = 100){
			generate_permutations_bernoulli_cpp(as.integer(private$t), as.integer(r), as.numeric(private$prob_T))$w_mat
		}
	),
	private = list(
		redraw_w_according_to_design = function(){
			private$w[1:private$t] = rbinom(private$t, 1, private$prob_T)
		}
	)
)
