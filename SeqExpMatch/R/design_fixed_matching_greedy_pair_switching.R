#' A Fixed Design Combining Binary Matching and Greedy Pair Switching
#'
#' @description
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
		#' @param num_cores The number of CPU cores to use.
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
				num_cores = 1,
				verbose = FALSE,
				max_designs = 100,
				objective = "mahal_dist",
				wait = TRUE,
				diff_method = FALSE
			) {
			if (prob_T != 0.5) {
				stop("FixedDesignMatchingGreedyPairSwitching only supports balanced designs (prob_T = 0.5).")
			}
			if (!requireNamespace("GreedyExperimentalDesign", quietly = TRUE)) {
				stop("Package 'GreedyExperimentalDesign' is required for FixedDesignMatchingGreedyPairSwitching.")
			}

			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)

			private$max_designs = max_designs
			private$objective = objective
			private$wait = wait
			private$diff_method = diff_method
			private$uses_covariates = TRUE
		},

		redraw_w_according_to_design = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		draw_ws_according_to_design = function(r = 100){
			assertCount(r, positive = TRUE)
			self$assert_experiment_completed()

			n = self$get_n()
			if (n %% 2 != 0) {
				stop("Binary-match-then-greedy designs require an even number of subjects.")
			}

			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X = private$X[1:n, , drop = FALSE]

			if (is.null(X) || ncol(X) == 0) {
				w_mat = replicate(r, sample(c(rep(1, n / 2), rep(0, n / 2))))
				return(matrix(as.numeric(w_mat), nrow = n, ncol = r))
			}

			search_obj = GreedyExperimentalDesign::initBinaryMatchFollowedByGreedyExperimentalDesignSearchObject(
				X = X,
				diff_method = private$diff_method,
				max_designs = r,
				objective = private$objective,
				wait = private$wait,
				start = TRUE,
				num_cores = private$num_cores,
				verbose = private$verbose
			)

			w_mat = tryCatch(
				GreedyExperimentalDesign::resultsBinaryMatchThenGreedySearch(
					search_obj,
					max_vectors = r,
					form = "one_zero"
				),
				error = function(e) {
					GreedyExperimentalDesign::resultsBinaryMatchThenGreedySearch(search_obj, max_vectors = r)
				}
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
			if (is.matrix(w_mat) && nrow(w_mat) != n && ncol(w_mat) == n) {
				w_mat = t(w_mat)
			}
			if (!is.matrix(w_mat) || nrow(w_mat) != n) {
				stop("resultsBinaryMatchThenGreedySearch returned an unexpected allocation matrix shape.")
			}

			if (all(w_mat %in% c(-1, 1), na.rm = TRUE)) {
				w_mat = (w_mat + 1) / 2
			}

			if (ncol(w_mat) < r) {
				stop("resultsBinaryMatchThenGreedySearch returned fewer allocation vectors than requested.")
			}
			w_mat[, seq_len(r), drop = FALSE]
		}
	)
)
