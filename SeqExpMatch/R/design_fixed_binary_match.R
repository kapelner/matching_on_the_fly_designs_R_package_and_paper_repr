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
				response_type = "continuous",
				prob_T = 0.5,
				mahal_match = FALSE,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				
				verbose = FALSE
			) {
			if (prob_T != 0.5){
				stop("Binary match designs only support even treatment allocation (prob_T = 0.5)")
			}
			assert_greedy_experimental_design_installed("FixedDesignBinaryMatch")
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
			assertCount(r, positive = TRUE)
			assert_greedy_experimental_design_installed("FixedDesignBinaryMatch")
			self$assert_all_subjects_arrived()
			private$ensure_bms_computed()
			n = self$get_n()
			if (is.null(private$m)){
				# No covariates: fall back to BCRD
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			# Cap at the total number of unique balanced allocations for small n
			max_designs_cap = if (n/2 < 30) floor(2^(n/2)) else r
			max_designs = min(r, max_designs_cap)

			search_obj = GreedyExperimentalDesign::initBinaryMatchExperimentalDesignSearchObject(
				private$bms,
				max_designs = max_designs,
				wait       = TRUE,
				start      = TRUE,
				num_cores  = self$num_cores,
				verbose    = private$verbose
			)
			w_mat = GreedyExperimentalDesign::resultsBinaryMatchSearch(search_obj, form = "one_zero")
			# resultsBinaryMatchSearch returns num_designs x n, we need n x num_designs
			w_mat = t(w_mat)
			storage.mode(w_mat) = "numeric"

			# If fewer unique designs exist than requested, recycle columns
			if (ncol(w_mat) < r){
				w_mat = w_mat[, rep(seq_len(ncol(w_mat)), length.out = r), drop = FALSE]
			}
			w_mat[, seq_len(r), drop = FALSE]
		}
	),

	private = list(
		mahal_match = NULL,
		bms = NULL,

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
