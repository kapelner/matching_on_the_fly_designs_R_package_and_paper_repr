#' A Factorial Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed factorial experimental design.
#' This design handles multiple treatment factors and balances assignments across
#' all factor combinations.
#'
#' @export
FixedDesignFactorial = R6::R6Class("FixedDesignFactorial",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a factorial fixed experimental design
		#'
		#' @param factors         A list where names are factor names and values are number of
		#'   levels (e.g. list(A=2, B=2)).
		#' @param response_type 	The data type of response values.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignFactorial` object
		initialize = function(
				factors,
				response_type,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				
				verbose = FALSE
			) {
			assertList(factors, types = "numeric", min.len = 1)
			# We don't use prob_T in the standard way here, as we have multiple factors
			# But base Design needs it. We'll set it to 0.5.
			super$initialize(response_type, 0.5, include_is_missing_as_a_new_feature, n, verbose)
			private$factors = factors
			
			# Precompute all combinations
			private$combinations = expand.grid(lapply(factors, function(l) 1:l))
			private$num_combinations = nrow(private$combinations)
		},


		#' @description
		#' Draw multiple treatment assignment vectors according to balanced factorial randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			n = self$get_n()
			
			# Standard balanced factorial: each combination appears n / num_combinations times
			base_alloc = rep(1:private$num_combinations, length.out = n)
			
			w_mat = matrix(NA_integer_, nrow = n, ncol = r)
			for (j in 1:r){
				w_mat[, j] = sample(base_alloc)
			}
			w_mat
		},

		#' @description
		#' Get the data frame of factor assignments for each subject.
		#'
		#' @return A data frame with n rows and columns corresponding to factors.
		get_w_factorial = function(){
			w_idx = self$get_w()
			if (length(w_idx) == 0 || any(is.na(w_idx))) return(NULL)
			private$combinations[w_idx, , drop = FALSE]
		}
	),

	private = list(
		factors = NULL,
		combinations = NULL,
		num_combinations = NULL
	)
)
