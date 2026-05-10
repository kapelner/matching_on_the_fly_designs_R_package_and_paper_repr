#' A balanced completely randomized Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed balanced
#' completely randomized experimental design.
#'
#' @examples
#' des = DesignFixediBCRD$new(n = 10, response_type = 'continuous')
#' des$add_all_subjects_to_experiment(data.frame(x1 = rnorm(10)))
#' des$assign_w_to_all_subjects()
#' @export
DesignFixediBCRD = R6::R6Class("DesignFixediBCRD",
	inherit = DesignFixed,
	public = list(
		#' @description
		#' Initialize a fixed balanced completely randomized experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param  prob_T  Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param  n  		The sample size.
		#' @param verbose A flag for verbosity.
		#' @param missingness_method How to handle missing values in covariates.
		#' @param model_formula A formula object.
		#'
		#' @return  A new `DesignFixediBCRD` object
		initialize = function(
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						
						verbose = FALSE,
				missingness_method = "impute",
				model_formula = ~ .
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose, missingness_method, model_formula)
			if (!is.null(n)) {
				private$m = rep(1L, as.integer(n))
			}
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to balanced randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r){
			generate_permutations_ibcrd_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T)
			)$w_mat
		}

	)
)
