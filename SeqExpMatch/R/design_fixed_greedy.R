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
			if (n %% 2 != 0){
				stop("FixedDesignGreedy requires an even number of subjects.")
			}
			if (is.null(private$X) || ncol(private$X) == 0){
				w_mat = replicate(r, sample(c(rep(1, n / 2), rep(0, n / 2))))
				return(private$validate_allocation_matrix(w_mat, n = n, r = r))
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
			private$validate_allocation_matrix(w_mat, n = n, r = r)
		}
	),

	private = list(
		objective = NULL,

		validate_allocation_matrix = function(w_mat, n, r){
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < r) {
				stop("GreedyExperimentalDesign returned an unexpected allocation matrix shape.")
			}

			w_mat = w_mat[, seq_len(r), drop = FALSE]
			storage.mode(w_mat) = "numeric"

			if (any(!is.finite(w_mat)) || any(is.na(w_mat))) {
				stop("GreedyExperimentalDesign returned non-finite treatment assignments.")
			}
			if (any(!(w_mat %in% c(0, 1)))) {
				stop("GreedyExperimentalDesign returned an invalid treatment assignment matrix.")
			}

			treated_counts = colSums(w_mat)
			if (any(treated_counts != n / 2)) {
				stop("GreedyExperimentalDesign returned an unbalanced allocation.")
			}

			w_mat
		}
	)
)
