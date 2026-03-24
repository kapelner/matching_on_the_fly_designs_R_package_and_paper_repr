#' A balanced completely randomized Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed balanced completely randomized experimental design.
#' This design ensures balanced allocation using the \pkg{randomizr} package.
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
			self$assert_all_subjects_arrived()
			# Use randomizr::complete_ra for canonical complete randomization
			w_mat = replicate(r, as.numeric(as.character(randomizr::complete_ra(N = self$get_n(), prob = private$prob_T))))
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	)
)
