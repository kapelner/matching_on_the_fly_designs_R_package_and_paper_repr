#' A Blocked and Cluster Randomized Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed blocked and cluster randomized experimental design.
#' This design randomizes clusters within specified blocks using the \pkg{randomizr} package.
#'
#' @export
FixedDesignBlockedCluster = R6::R6Class("FixedDesignBlockedCluster",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a blocked and cluster randomized fixed experimental design
		#'
		#' @param strata_cols 	A character vector of column names to use for stratification (blocks).
		#' @param cluster_col 	The column name in the data that identifies the cluster for each subject.
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment for each cluster.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignBlockedCluster` object
		#'
		initialize = function(
				strata_cols,
				cluster_col,
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			assertCharacter(strata_cols, min.len = 1)
			assertCharacter(cluster_col, len = 1)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$strata_cols = strata_cols
			private$cluster_col = cluster_col
			private$uses_covariates = TRUE
		},

		#' @description
		#' Redraw treatment assignments according to the blocked cluster randomized design.
		redraw_w_according_to_design = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to blocked cluster randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			
			strata_keys = vapply(1:private$t, function(i) {
				vals = vapply(private$strata_cols, function(col) {
					val = private$Xraw[i, ][[col]]
					if (is.na(val)) "NA" else as.character(val)
				}, character(1))
				paste(vals, collapse = "|")
			}, character(1))
			
			cluster_ids = as.character(private$Xraw[[private$cluster_col]])
			
			# Use randomizr::block_and_cluster_ra for canonical blocked and clustered randomization
			w_mat = replicate(r, as.numeric(as.character(randomizr::block_and_cluster_ra(
				blocks = strata_keys, 
				clusters = cluster_ids, 
				prob = private$prob_T
			))))
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),

	private = list(
		strata_cols = NULL,
		cluster_col = NULL
	)
)
