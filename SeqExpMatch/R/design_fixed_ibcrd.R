#' A balanced completely randomized Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed balanced completely randomized experimental design.
#'
#' @export
FixedDesigniBCRD = R6::R6Class("FixedDesigniBCRD",
	inherit = FixedDesign,
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

		draw_ws_according_to_design = function(r = 100){
			generate_permutations_ibcrd_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		}
	)
)
