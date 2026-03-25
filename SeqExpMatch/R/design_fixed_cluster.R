#' A Cluster Randomized Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed cluster randomized experimental design.
#' This design randomizes entire groups (clusters) of subjects together using the \pkg{randomizr} package.
#'
#' @export
FixedDesignCluster = R6::R6Class("FixedDesignCluster",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a cluster randomized fixed experimental design
		#'
		#' @param cluster_col 	The column name in the data that identifies the cluster for each subject.
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment for each cluster.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignCluster` object
		#'
		initialize = function(
				cluster_col,
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			assertCharacter(cluster_col, len = 1)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$cluster_col = cluster_col
			private$uses_covariates = TRUE
		},

		#' @description
		#' Redraw treatment assignments according to the cluster randomized design.
		redraw_w_according_to_design = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to cluster randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			cluster_ids = as.character(private$Xraw[[private$cluster_col]])
			if (any(is.na(cluster_ids))){
				stop("Cluster IDs cannot be missing.")
			}
			
			# Use randomizr::cluster_ra for canonical cluster randomization
			w_mat = replicate(r, as.numeric(as.character(randomizr::cluster_ra(clusters = cluster_ids, prob = private$prob_T))))
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),

	private = list(
		cluster_col = NULL
	)
)
