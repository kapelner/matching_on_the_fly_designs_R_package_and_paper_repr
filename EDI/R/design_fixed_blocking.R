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
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `FixedDesignBlocking` object
		initialize = function(
						strata_cols,
						response_type,
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						verbose = FALSE
					) {
			assertCharacter(strata_cols, min.len = 1)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			private$strata_cols = strata_cols
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
			
			strata_keys = vapply(1:private$t, function(i) {
				vals = vapply(private$strata_cols, function(col) {
					val = private$Xraw[i, ][[col]]
					if (is.na(val)) "NA" else as.character(val)
				}, character(1))
				paste(vals, collapse = "|")
			}, character(1))
			
			# Use randomizr::block_ra for canonical stratified blocking if available,
			# or fallback to our C++ implementation.
			if (requireNamespace("randomizr", quietly = TRUE)) {
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
		strata_cols = NULL,

		get_strata_keys = function(){
			n = private$t
			vapply(1:n, function(i) {
				vals = vapply(private$strata_cols, function(col) {
					val = private$Xraw[i, ][[col]]
					if (is.na(val)) "NA" else as.character(val)
				}, character(1))
				paste(vals, collapse = "|")
			}, character(1))
		},

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
