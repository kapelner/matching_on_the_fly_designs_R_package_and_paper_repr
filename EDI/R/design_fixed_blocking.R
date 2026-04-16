#' A stratified blocking Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed stratified
#' blocking experimental design.
#'
#' @export
FixedDesignBlocking = R6::R6Class("FixedDesignBlocking",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a fixed stratified blocking experimental design
		#'
		#' @param strata_cols A character vector of column names to use for stratification.
		#'   If `NULL` (the default), all available covariate columns are used.
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_bins_for_continuous_covariate The number of quantile bins to use for continuous strata. Default is 2.
		#' @param B_preferred The desired number of blocks. Columns from `strata_cols`
		#'   are added greedily in order, each column being included only if it does not push
		#'   the total number of unique blocks beyond this target. For categorical covariates
		#'   their natural levels are used; for continuous covariates
		#'   `num_bins_for_continuous_covariate` quantile bins are used. Earlier columns
		#'   are always preferred over later ones. The default is `floor(sqrt(n))` when `n`
		#'   is known at construction time, or is resolved to `floor(sqrt(n))` when subjects
		#'   are added. Set to `NULL` to disable the target and use all `strata_cols` columns
		#'   unconditionally (the original behaviour).
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `FixedDesignBlocking` object
		initialize = function(
						strata_cols = NULL,
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_bins_for_continuous_covariate = 2,
						B_preferred = if (!is.null(n)) max(1L, floor(sqrt(n))) else NA_integer_,
						verbose = FALSE
					) {
			if (!is.null(strata_cols)) assertCharacter(strata_cols, min.len = 1)
			assertCount(num_bins_for_continuous_covariate, positive = TRUE)
			if (!is.null(B_preferred) && !is.na(B_preferred)) assertCount(B_preferred, positive = TRUE)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$strata_cols = strata_cols
			private$num_bins_for_continuous_covariate = num_bins_for_continuous_covariate
			private$B_preferred = B_preferred
			private$uses_covariates = TRUE
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to stratified blocking.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			
			strata_keys = private$get_strata_keys()
			
			# Use randomizr::block_ra for canonical stratified blocking if available,
			# or fallback to our C++ implementation.
			if (check_package_installed("randomizr")) {
				w_mat = replicate(r, as.numeric(as.character(randomizr::block_ra(blocks = strata_keys, prob = private$prob_T))))
				storage.mode(w_mat) = "numeric"
				return(w_mat)
			}
			
			unique_keys = unique(strata_keys)
			strata_indices = lapply(unique_keys, function(key) which(strata_keys == key))
			
			res = generate_permutations_blocking_cpp(
				as.integer(self$get_n()),
				as.integer(r),
				as.numeric(private$prob_T),
				strata_indices
			)
			w_mat = res$w_mat
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),
	private = list(
		draw_bootstrap_indices = function(bootstrap_type = NULL){
			strata_keys = private$get_strata_keys()
			if (is.null(bootstrap_type) || bootstrap_type == "within_blocks") {
				list(i_b = stratified_bootstrap_indices_cpp(as.character(strata_keys)), m_vec_b = NULL)
			} else {
				unique_keys = unique(strata_keys)
				sampled_keys = sample(unique_keys, length(unique_keys), replace = TRUE)
				i_b = unlist(lapply(sampled_keys, function(key) which(strata_keys == key)), use.names = FALSE)
				list(i_b = as.integer(i_b), m_vec_b = NULL)
			}
		}
	)
)
