#' A Greedy Search Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed greedy experimental design.
#' Uses the \pkg{GreedyExperimentalDesign} package to search for balanced allocations.
#'
#' @export
FixedDesignGreedy = R6::R6Class("FixedDesignGreedy",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a greedy search fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment. Must be 0.5.
		#' @param objective 	The objective function to use. Default is "mahal_dist".
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignGreedy` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				objective = "mahal_dist",
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				
				verbose = FALSE
			) {
			if (prob_T != 0.5){
				stop("Greedy designs currently only support even treatment allocation (prob_T = 0.5)")
			}
			assert_greedy_experimental_design_installed("FixedDesignGreedy")
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$objective = objective
			private$uses_covariates = TRUE
		},


		#' @description
		#' Draw multiple treatment assignment vectors.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			assertCount(r, positive = TRUE)
			assert_greedy_experimental_design_installed("FixedDesignGreedy")
			self$assert_all_subjects_arrived()
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X = private$X[1:n, , drop = FALSE]

			search_obj = GreedyExperimentalDesign::initGreedyExperimentalDesignObject(
				X          = X,
				max_designs = r,
				objective  = private$objective,
				wait       = TRUE,
				start      = TRUE,
				num_cores  = self$num_cores,
				verbose    = private$verbose
			)
			w_mat = GreedyExperimentalDesign::resultsGreedySearch(search_obj, max_vectors = r, form = "one_zero")$ending_indicTs
			storage.mode(w_mat) = "numeric"
			w_mat[, seq_len(min(r, ncol(w_mat))), drop = FALSE]
		}
	),

	private = list(
		objective = NULL
	)
)
