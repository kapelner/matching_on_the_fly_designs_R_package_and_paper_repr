#' A Fixed Design Combining Binary Matching and Greedy Pair Switching
#'
#' An R6 class encapsulating a fixed experimental design that first computes
#' binary matches and then improves the matched design with greedy pair switching
#' using the \pkg{GreedyExperimentalDesign} package.
#'
#' @export
FixedDesignMatchingGreedyPairSwitching = R6::R6Class("FixedDesignMatchingGreedyPairSwitching",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a fixed design that performs binary matching followed by greedy pair switching.
		#'
		#' @param response_type The data type of response values.
		#' @param prob_T The probability of treatment assignment. Must be \code{0.5}.
		#' @param include_is_missing_as_a_new_feature Flag for missingness indicators.
		#' @param n The sample size.
		#' @param verbose A flag for verbosity.
		#' @param max_designs The number of searched designs to retain.
		#' @param objective The imbalance objective passed to \pkg{GreedyExperimentalDesign}.
		#' @param wait If \code{TRUE}, wait for the search to finish before returning.
		#' @param diff_method Passed to the upstream binary-match-then-greedy initializer.
		#'
		#' @return A new \code{FixedDesignMatchingGreedyPairSwitching} object.
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n,
				
				verbose = FALSE,
				max_designs = 100,
				objective = "mahal_dist",
				wait = TRUE,
				diff_method = FALSE
			) {
			if (prob_T != 0.5) {
				stop("FixedDesignMatchingGreedyPairSwitching only supports balanced designs (prob_T = 0.5).")
			}
			assert_greedy_experimental_design_installed("FixedDesignMatchingGreedyPairSwitching")

			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)

			private$max_designs = max_designs
			private$objective = objective
			private$wait = wait
			private$diff_method = diff_method
			private$uses_covariates = TRUE
		},


		#' @description
		#' Draw multiple treatment assignment vectors according to binary match followed by
		#' greedy switching.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			assertCount(r, positive = TRUE)
			assert_greedy_experimental_design_installed("FixedDesignMatchingGreedyPairSwitching")
			self$assert_all_subjects_arrived()

			n = self$get_n()
			if (n %% 2 != 0) {
				stop("Binary-match-then-greedy designs require an even number of subjects.")
			}

			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X = private$X[1:n, , drop = FALSE]
			search_budget = max(as.integer(r), as.integer(private$max_designs))

			search_obj = GreedyExperimentalDesign::initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject(
				X = X,
				diff_method = private$diff_method,
				max_designs = search_budget,
				objective = private$objective,
				wait = private$wait,
				start = TRUE,
				num_cores = self$num_cores,
				verbose = private$verbose
			)

			w_mat = GreedyExperimentalDesign::resultsBinaryMatchThenGreedySearch(
				search_obj,
				max_vectors = r,
				form = "one_zero"
			)

			w_mat = private$extract_allocation_matrix(w_mat, n = n, r = r)
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),

	private = list(
		max_designs = NULL,
		objective = NULL,
		wait = NULL,
		diff_method = NULL,
		
		draw_one_w = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		extract_allocation_matrix = function(res, n, r){
			if (is.matrix(res)) {
				w_mat = res
			} else if (is.list(res)) {
				if (!is.null(res$ending_indicTs)) {
					w_mat = res$ending_indicTs
				} else if (!is.null(res$designs)) {
					w_mat = res$designs
				} else if (!is.null(res$indicTs)) {
					w_mat = res$indicTs
				} else {
					stop("resultsBinaryMatchThenGreedySearch returned an unsupported result structure.")
				}
			} else if (is.vector(res)) {
				w_mat = matrix(res, nrow = n, ncol = 1)
			} else {
				stop("resultsBinaryMatchThenGreedySearch returned an unsupported result type.")
			}

			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < r) {
				stop("resultsBinaryMatchThenGreedySearch returned an unexpected allocation matrix shape.")
			}
			storage.mode(w_mat) = "numeric"
			w_mat[, seq_len(r), drop = FALSE]
		}
	)
)
