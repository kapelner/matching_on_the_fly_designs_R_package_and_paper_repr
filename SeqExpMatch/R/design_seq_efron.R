#' Efron's (1971) Biased Coin Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments.
#'
#' @export
SeqDesignEfron = R6::R6Class("SeqDesignEfron",
	inherit = SeqDesign,
	public = list(
		initialize = function(
						response_type = "continuous",
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_cores = 1,
						verbose = FALSE,
						weighted_coin_prob = 2/3
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$weighted_coin_prob = weighted_coin_prob
		},

		assign_wt = function(){
			n_T = sum(private$w, na.rm = TRUE)
			n_C = private$t - n_T
			if (n_T * private$prob_T > n_C * (1 - private$prob_T)){
				rbinom(1, 1, 1 - private$weighted_coin_prob)
			} else if (n_T * private$prob_T < n_C * (1 - private$prob_T)){
				rbinom(1, 1, private$weighted_coin_prob)
			} else {
				rbinom(1, 1, private$prob_T)
			}
		},

		draw_ws_according_to_design = function(r = 100){
			generate_permutations_efron_cpp(as.integer(private$t), as.integer(r), as.numeric(private$prob_T), as.numeric(private$weighted_coin_prob))$w_mat
		}
	),
	private = list(
		weighted_coin_prob = NULL,

		redraw_w_according_to_design = function(){
			private$w[1:private$t] = efron_redraw_cpp(private$t, private$prob_T, private$weighted_coin_prob)
		}
	)
)
