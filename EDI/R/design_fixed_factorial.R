#' A Factorial Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed factorial experimental design.
#' This design handles one or more treatment factors and balances assignments across
#' all factor combinations.
#'
#' Currently restricted to exactly two total factor-level combinations (i.e. two arms),
#' e.g. a single two-level factor. \code{w} follows the same \{0,1\} internal /
#' \{-1,+1\} public convention as every other \code{Design} subclass, so this design
#' works unmodified with every \code{Inference} class. Support for more than two
#' combinations is tracked separately (see \code{package_metadata/new_feature_plans/multi_arm_designs.md}).
#'
#' @examples
#' des = DesignFixedFactorial$new(n = 12, response_type = 'continuous', factors = list(treatment = 2))
#' des$add_all_subjects_to_experiment(data.frame(x=1:12))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixedFactorial = R6::R6Class("DesignFixedFactorial",
	inherit = DesignFixed,
	public = list(
		#' @description Initialize a factorial fixed experimental design
		#'
		#' @param factors         A list where names are factor names and values are number of
		#'   levels (e.g. list(treatment = 2)). The product of levels across all factors must
		#'   currently equal exactly 2 (two-arm only).
		#' @param response_type 	The data type of response values.
		#' @param include_is_missing_as_a_new_feature  Flag for missingness indicators.
		#' @param n  		The sample size.
		#' @param verbose  Flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param design_formula A formula object.
		#' @param seed Integer seed for reproducibility.
		#'
		#' @return 			A new `DesignFixedFactorial` object
		initialize = function(
				factors,
				response_type,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,

				verbose = FALSE,
				missingness_method = "impute",
				design_formula = ~ .,
				seed = NULL
			) {
			if (should_run_asserts()) {
				assertList(factors, types = "numeric", min.len = 1)
			}
			# We don't use prob_T in the standard way here, as we have multiple factors
			# But base Design needs it. We'll set it to 0.5.
			super$initialize(response_type, 0.5, include_is_missing_as_a_new_feature, n, verbose, missingness_method, design_formula, seed = seed)
			private$factors = factors

			# Precompute all combinations
			private$combinations = expand.grid(lapply(factors, function(l) 1:l))
			private$num_combinations = nrow(private$combinations)
			if (should_run_asserts()) {
				if (private$num_combinations != 2L) {
					stop(
						"DesignFixedFactorial currently only supports exactly two total factor-level ",
						"combinations (two-arm only), e.g. a single two-level factor. The supplied ",
						"'factors' argument implies ", private$num_combinations, " combinations."
					)
				}
			}
		},
		#' @description Get the data frame of factor assignments for each subject.
		#'
		#' @return A data frame with n rows and columns corresponding to factors.
		get_w_factorial = function(){
			w_idx = private$w
			if (length(w_idx) == 0 || any(is.na(w_idx))) return(NULL)
			private$combinations[w_idx + 1L, , drop = FALSE]
		}
	),
	private = list(
		factors = NULL,
		combinations = NULL,
		num_combinations = NULL,
		draw_ws_raw = function(r = 100){
			private$maybe_set_seed()
			if (should_run_asserts()) {
				self$assert_all_subjects_arrived()
			}
			n = self$get_n()
			# Standard balanced two-arm factorial: each of the two combinations appears n/2 times.
			# {0,1} storage matches the shared Design/DesignFixed convention (see get_w()/
			# assign_w_to_all_subjects() on the base classes, which this class now inherits
			# unmodified) so every Inference class works against this design unchanged.
			base_alloc = rep(0:(private$num_combinations - 1L), length.out = n)
			w_mat = matrix(NA_integer_, nrow = n, ncol = r)
			for (j in 1:r){
				w_mat[, j] = sample(base_alloc)
			}
			w_mat
		}
	)
)
