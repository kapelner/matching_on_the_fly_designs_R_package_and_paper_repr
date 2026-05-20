#' A Fixed Design Combining Binary Matching and Greedy Pair Switching
#'
#' An R6 class encapsulating a fixed experimental design that first computes
#' binary matches and then improves the matched design with greedy pair switching
#' using a native C++ implementation.
#'
#' @examples
#' \dontrun{
#' des = DesignFixedMatchingGreedyPairSwitching$new(n = 10, response_type = 'continuous')
#' }
#' @export
DesignFixedMatchingGreedyPairSwitching = R6::R6Class("DesignFixedMatchingGreedyPairSwitching",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a fixed design that performs binary matching followed by greedy pair switching.
		#'
		#' @param response_type The data type of response values.
		#' @param prob_T The probability of treatment assignment. Must be \code{0.5}.
		#' @param include_is_missing_as_a_new_feature Flag for missingness indicators.
		#' @param n The sample size.
		#' @param verbose A flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param model_formula A formula object.
		#' @param objective The imbalance objective. Either \code{"mahal_dist"} (default) or \code{"abs_sum_diff"}.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return A new \code{DesignFixedMatchingGreedyPairSwitching} object.
		initialize = function(
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n,
				verbose = FALSE,
				objective = "mahal_dist",
				missingness_method = "impute",
				model_formula = ~ .,
				seed = NULL
			) {
			if (should_run_asserts()) {
				if (prob_T != 0.5) {
					stop("DesignFixedMatchingGreedyPairSwitching only supports balanced designs (prob_T = 0.5).")
				}
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, model_formula, seed = seed)
			private$objective = objective
			private$uses_covariates = TRUE
		},
		#' @description Draw multiple treatment assignment vectors according to binary match followed by
		#' greedy pair switching.
		#'
		#' @param r The number of designs to draw.
		#'
		#' @return A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			if (should_run_asserts()) {
				if (n %% 4 != 0) {
					stop("DesignFixedMatchingGreedyPairSwitching requires n divisible by 4.")
				}
			}
			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X = private$X[1:n, , drop = FALSE]
			if (is.null(private$bms)) {
				private$bms = compute_binary_match_structure(X, mahal_match = (private$objective == "mahal_dist"))
			}
			pairs_mat = private$bms$indicies_pairs
			storage.mode(pairs_mat) = "integer"
			n_iter = max(500L, 500L * as.integer(n))
			w_mat = greedy_design_search_cpp(
				X_raw          = X,
				r              = as.integer(r),
				objective      = private$objective,
				n_iter         = n_iter,
				indicies_pairs = pairs_mat
			)
			storage.mode(w_mat) = "numeric"
			private$validate_allocation_matrix(w_mat, n = n, r = r)
		}
	),
	private = list(
		objective = NULL,
		bms = NULL,
		validate_allocation_matrix = function(w_mat, n, r){
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (should_run_asserts()) {
				if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < r) {
					stop("DesignFixedMatchingGreedyPairSwitching returned an unexpected allocation matrix shape.")
				}
			}
			w_mat = w_mat[, seq_len(r), drop = FALSE]
			storage.mode(w_mat) = "numeric"
			if (should_run_asserts()) {
				if (any(!is.finite(w_mat)) || any(is.na(w_mat))) {
					stop("DesignFixedMatchingGreedyPairSwitching returned non-finite treatment assignments.")
				}
				if (any(!(w_mat %in% c(0, 1)))) {
					stop("DesignFixedMatchingGreedyPairSwitching returned an invalid treatment assignment matrix.")
				}
				treated_counts = colSums(w_mat)
				if (any(treated_counts != n / 2)) {
					stop("DesignFixedMatchingGreedyPairSwitching returned an unbalanced allocation.")
				}
			}
			w_mat
		}
	)
)
