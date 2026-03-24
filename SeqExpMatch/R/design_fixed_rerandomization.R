#' A Rerandomization Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed rerandomization
#' experimental design. This class is a thin wrapper around the
#' \pkg{GreedyExperimentalDesign} rerandomization search API.
#'
#' @export
FixedDesignRerandomization = R6::R6Class("FixedDesignRerandomization",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a rerandomization fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment.
		#' @param obj_val_cutoff 	The maximum allowable objective value.
		#' @param objective 	The objective function to use. Default is "mahal_dist".
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignRerandomization` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				obj_val_cutoff = NULL,
				objective = "mahal_dist",
				include_is_missing_as_a_new_feature = TRUE,
				n,
				num_cores = 1,
				verbose = FALSE
			) {
			assert_greedy_experimental_design_installed("FixedDesignRerandomization")
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$obj_val_cutoff = obj_val_cutoff
			private$objective = objective
			private$uses_covariates = TRUE
		},

		redraw_w_according_to_design = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		draw_ws_according_to_design = function(r = 100){
			assertCount(r, positive = TRUE)
			assert_greedy_experimental_design_installed("FixedDesignRerandomization")
			self$assert_experiment_completed()

			n = self$get_n()
			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			X = private$X[1:n, , drop = FALSE]
			search_obj = GreedyExperimentalDesign::initRerandomizationExperimentalDesignObject(
				X = X,
				obj_val_cutoff_to_include = private$obj_val_cutoff,
				max_designs = r,
				objective = private$objective,
				wait = TRUE,
				start = TRUE,
				num_cores = private$num_cores,
				verbose = private$verbose
			)
			res = GreedyExperimentalDesign::resultsRerandomizationSearch(search_obj)
			if (!is.list(res) || is.null(res$designs)) {
				stop("resultsRerandomizationSearch returned an unsupported result structure.")
			}
			w_mat = res$designs
			if (is.vector(w_mat)) {
				w_mat = matrix(w_mat, nrow = n, ncol = 1)
			} else if (nrow(w_mat) != n && ncol(w_mat) == n) {
				w_mat = t(w_mat)
			}
			if (all(w_mat %in% c(-1, 1), na.rm = TRUE)) {
				w_mat = (w_mat + 1) / 2
			}
			if (!is.matrix(w_mat) || nrow(w_mat) != n || ncol(w_mat) < r) {
				stop("resultsRerandomizationSearch returned an unexpected allocation matrix shape.")
			}
			storage.mode(w_mat) = "numeric"
			w_mat[, seq_len(r), drop = FALSE]
		}
	),

	private = list(
		obj_val_cutoff = NULL,
		objective = NULL
	)
)
