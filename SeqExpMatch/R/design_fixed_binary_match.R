#' A Binary Match Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed binary match experimental design.
#' This class is a thin wrapper around the \pkg{GreedyExperimentalDesign}
#' binary-match search API.
#'
#' @export
FixedDesignBinaryMatch = R6::R6Class("FixedDesignBinaryMatch",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a binary match fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment. Must be 0.5.
		#' @param mahal_match 	Match using Mahalanobis distance. Default is \code{FALSE} (Euclidean).
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignBinaryMatch` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				mahal_match = FALSE,
				include_is_missing_as_a_new_feature = TRUE,
				n,
				num_cores = 1,
				verbose = FALSE
			) {
			if (prob_T != 0.5){
				stop("Binary match designs only support even treatment allocation (prob_T = 0.5)")
			}
			assert_greedy_experimental_design_installed("FixedDesignBinaryMatch")
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$mahal_match = mahal_match
			private$uses_covariates = TRUE
		},

		redraw_w_according_to_design = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		draw_ws_according_to_design = function(r = 100){
			assertCount(r, positive = TRUE)
			assert_greedy_experimental_design_installed("FixedDesignBinaryMatch")
			self$assert_experiment_completed()

			n = self$get_n()
			if (n %% 2 != 0){
				stop("Binary match designs require an even number of subjects.")
			}

			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X = private$X[1:n, , drop = FALSE]
			binary_match_structure = GreedyExperimentalDesign::computeBinaryMatchStructure(
				X = X,
				mahal_match = private$mahal_match
			)
			search_obj = GreedyExperimentalDesign::initBinaryMatchExperimentalDesignSearchObject(
				binary_match_structure = binary_match_structure,
				max_designs = r,
				wait = TRUE,
				start = TRUE,
				num_cores = private$num_cores,
				verbose = private$verbose
			)
			w_mat = GreedyExperimentalDesign::resultsBinaryMatchSearch(
				search_obj,
				form = "one_zero"
			)
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < r) {
				stop("resultsBinaryMatchSearch returned an unexpected allocation matrix shape.")
			}
			storage.mode(w_mat) = "numeric"
			w_mat[, seq_len(r), drop = FALSE]
		}
	),

	private = list(
		mahal_match = NULL
	)
)
