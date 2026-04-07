#' A Blocked and Cluster Randomized Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed blocked and
#' cluster randomized experimental design.
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
		#' @param num_bins_for_continuous_covariate The number of quantile bins to use for continuous strata. Default is 2.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignBlockedCluster` object
		#'
		initialize = function(
				strata_cols,
				cluster_col,
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_bins_for_continuous_covariate = 2,
				verbose = FALSE
			) {
			assertCharacter(strata_cols, min.len = 1)
			assertCharacter(cluster_col, len = 1)
			assertCount(num_bins_for_continuous_covariate, positive = TRUE)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$strata_cols = strata_cols
			private$cluster_col = cluster_col
			private$num_bins_for_continuous_covariate = num_bins_for_continuous_covariate
			private$uses_covariates = TRUE
		},


		#' @description
		#' Draw multiple treatment assignment vectors according to blocked cluster randomization.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			
			strata_keys = private$get_strata_keys()
			
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
		cluster_col = NULL,

		draw_bootstrap_indices = function(bootstrap_type = NULL){
			n = private$t
			strata_keys = private$get_strata_keys()
			cluster_ids = as.character(private$Xraw[1:n, ][[private$cluster_col]])

			if (is.null(bootstrap_type) || bootstrap_type == "within_blocks") {
				# Resample clusters within each stratum
				unique_strata = unique(strata_keys)
				i_b = unlist(lapply(unique_strata, function(stratum) {
					stratum_idx = which(strata_keys == stratum)
					stratum_clusters = unique(cluster_ids[stratum_idx])
					sampled_clusters = sample(stratum_clusters, length(stratum_clusters), replace = TRUE)
					unlist(lapply(sampled_clusters, function(cl) which(cluster_ids == cl)), use.names = FALSE)
				}), use.names = FALSE)
			} else {
				# Resample blocks (strata) themselves
				unique_strata = unique(strata_keys)
				sampled_strata = sample(unique_strata, length(unique_strata), replace = TRUE)
				i_b = unlist(lapply(sampled_strata, function(stratum) which(strata_keys == stratum)), use.names = FALSE)
			}
			list(i_b = as.integer(i_b), m_vec_b = NULL)
		}
	)
)
