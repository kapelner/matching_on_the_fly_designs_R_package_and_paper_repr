#' A Binary Match Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed binary match
#' experimental design.
#' This design pairs subjects based on covariate distances and randomizes within pairs.
#' Uses the \pkg{GreedyExperimentalDesign} package for distance computation and random allocation.
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
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignBinaryMatch` object
		#'
		initialize = function(
				response_type,
				prob_T = 0.5,
				mahal_match = FALSE,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,				
				verbose = FALSE
			) {
			if (should_run_asserts()) {
				if (prob_T != 0.5){
					stop("Binary match designs only support even treatment allocation (prob_T = 0.5)")
				}
				assert_greedy_experimental_design_installed("FixedDesignBinaryMatch")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$mahal_match = mahal_match
			private$uses_covariates = TRUE
		},


		#' @description
		#' Draw multiple treatment assignment vectors according to binary matching.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
				self$assert_all_subjects_arrived()
			}
			private$ensure_bms_computed()
			n = self$get_n()
			if (is.null(private$m)){
				# No covariates: fall back to BCRD
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			# Use the in-house Rcpp search instead of GreedyExperimentalDesign Java search
			w_mat = draw_binary_match_assignments_cpp(
				private$bms$indicies_pairs,
				as.integer(n),
				as.integer(r),
				as.integer(self$num_cores)
			)
			private$validate_allocation_matrix(w_mat, n = n, r = r)
		}
	),

	private = list(
		mahal_match = NULL,
		bms = NULL,

		draw_bootstrap_indices = function(bootstrap_type = NULL){
			private$ensure_bms_computed()
			if (is.null(private$m)){
				n = self$get_n()
				return(list(i_b = sample.int(n, n, replace = TRUE), m_vec_b = NULL))
			}
			.draw_kk_bootstrap_indices(private)
		},

		validate_allocation_matrix = function(w_mat, n, r){
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			}
			if (should_run_asserts()) {
				if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < 1L) {
					stop("resultsBinaryMatchSearch returned an unexpected allocation matrix shape.")
				}
			}

			storage.mode(w_mat) = "numeric"
			if (should_run_asserts()) {
				if (any(!is.finite(w_mat)) || any(is.na(w_mat))) {
					stop("resultsBinaryMatchSearch returned non-finite treatment assignments.")
				}
				if (any(!(w_mat %in% c(0, 1)))) {
					stop("resultsBinaryMatchSearch returned an invalid treatment assignment matrix.")
				}

				treated_counts = colSums(w_mat)
				if (any(treated_counts != n / 2)) {
					stop("resultsBinaryMatchSearch returned an unbalanced allocation.")
				}
			}

			if (ncol(w_mat) < r){
				w_mat = w_mat[, rep(seq_len(ncol(w_mat)), length.out = r), drop = FALSE]
			}
			w_mat[, seq_len(r), drop = FALSE]
		},

		ensure_bms_computed = function(){
			n = self$get_n()
			if (is.null(private$bms) && !is.null(private$X) && ncol(private$X) > 0){
				X = private$X[1:n, , drop = FALSE]
				private$bms = GreedyExperimentalDesign::computeBinaryMatchStructure(X, mahal_match = private$mahal_match)

				# Build pair-ID vector m where m[i] = pair index for subject i
				m_vec = integer(n)
				pairs = private$bms$indicies_pairs
				for (i in seq_len(nrow(pairs))){
					m_vec[pairs[i, 1]] = i
					m_vec[pairs[i, 2]] = i
				}
				private$m = m_vec
			}
		}
	)
)
