#' A D-optimal Search Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed D-optimal experimental design.
#' This design searches for an allocation that maximizes the determinant of the information matrix
#' (equivalent to minimizing the variance of the parameter estimates).
#'
#' @export
FixedDesignDOptimal = R6::R6Class("FixedDesignDOptimal",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a D-optimal search fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment. Must be 0.5 for exchange search.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignDOptimal` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				
				verbose = FALSE
			) {
			if (prob_T != 0.5){
				stop("D-optimal exchange search currently only supports even treatment allocation (prob_T = 0.5)")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$uses_covariates = TRUE
		},


		#' @description
		#' Draw treatment assignments according to the D-optimal search fixed design.
		#'
		#' @param r Number of designs to draw.
		#'
		#' @return A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			# Precompute projection matrix P = Z0 (Z0' Z0)^-1 Z0'
			if (is.null(private$P)){
				X = private$X[1:n, , drop = FALSE]
				Z0 = cbind(1, X)
				Z0_qr = qr(Z0)
				Q = qr.Q(Z0_qr)
				private$P = Q %*% t(Q)
			}

			# Use C++ speedup
			# Note: we don't parallelize here because C++ handles the loop over r
			# and it's already very fast. If r is huge, we could parallelize outside.
			# But generate_permutations... style usually does it in one call.
			res = d_optimal_search_cpp(private$P, as.integer(r), as.integer(round(n * private$prob_T)))
			w_mat = res
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),

	private = list(
		P = NULL # Projection matrix
	)
)
