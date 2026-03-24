#' Pocock & Simon's Minimization Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a Pocock & Simon sequential experimental design.
#' This design minimizes the imbalance across treatments for multiple covariates.
#'
#' @export
SeqDesignPocockSimon = R6::R6Class("SeqDesignPocockSimon",
	inherit = SeqDesign,
	public = list(
		#' @description
		#' Initialize a Pocock & Simon sequential experimental design
		#'
		#' @param strata_cols 	The names of the covariates to be used for minimization. These must be factor or categorical variables.
		#' @param weights 		A numeric vector of weights for each covariate. Defaults to 1 for all.
		#' @param p_best 		The probability of assigning the treatment that minimizes the imbalance. Defaults to 0.8.
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `SeqDesignPocockSimon` object
		#'
		initialize = function(
				strata_cols,
				weights = NULL,
				p_best = 0.8,
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			assertCharacter(strata_cols, min.len = 1)
			assertNumeric(p_best, lower = 0.5, upper = 1)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			
			private$strata_cols = strata_cols
			private$p_best = p_best
			private$uses_covariates = TRUE
			
			if (is.null(weights)){
				private$weights = rep(1, length(strata_cols))
			} else {
				assertNumeric(weights, len = length(strata_cols), lower = 0)
				private$weights = weights
			}
		},

		assign_wt = function(){
			private$ensure_factor_metadata()
			subject_levels_idx = private$get_subject_levels_idx(private$Xraw[private$t, ])
			
			if (is.null(private$counts)){
				private$counts = matrix(0, nrow = private$num_levels_total, ncol = 2)
			}
			
			# Call Rcpp function that assigns and updates counts in-place
			pocock_simon_assign_and_update_cpp(
				private$counts, 
				as.integer(subject_levels_idx), 
				private$weights, 
				private$p_best, 
				private$prob_T
			)
		}
	),

	private = list(
		strata_cols = NULL,
		weights = NULL,
		p_best = NULL,
		factor_levels = list(), 
		num_levels_total = NULL,
		level_offsets = NULL,   
		counts = NULL, # State: num_levels_total x 2

		ensure_factor_metadata = function(){
			if (length(private$factor_levels) == 0){
				offset = 0
				private$level_offsets = numeric(length(private$strata_cols))
				
				for (i in 1 : length(private$strata_cols)){
					col = private$strata_cols[i]
					vals = private$Xraw[[col]]
					lvls = levels(as.factor(vals))
					if (length(lvls) == 0) lvls = "NA"
					private$factor_levels[[col]] = lvls
					private$level_offsets[i] = offset
					offset = offset + length(lvls)
				}
				private$num_levels_total = offset
			}
		},

		get_subject_levels_idx = function(x_row){
			indices = numeric(length(private$strata_cols))
			for (i in 1 : length(private$strata_cols)){
				col = private$strata_cols[i]
				val = as.character(x_row[[col]])
				if (is.na(val)) val = "NA"
				
				# Use match() for faster lookup
				lvl_idx = match(val, private$factor_levels[[col]])
				if (is.na(lvl_idx)) lvl_idx = 1
				
				indices[i] = private$level_offsets[i] + lvl_idx
			}
			indices
		},

		redraw_w_according_to_design = function(){
			private$ensure_factor_metadata()
			# Reset incremental counts
			private$counts = matrix(0, nrow = private$num_levels_total, ncol = 2)
			
			n = private$t
			x_levels_matrix = matrix(NA_integer_, nrow = n, ncol = length(private$strata_cols))
			for (i in 1 : n){
				x_levels_matrix[i, ] = private$get_subject_levels_idx(private$Xraw[i, ])
			}
			
			# Redraw entire vector using optimized C++ loop
			private$w[1:n] = pocock_simon_redraw_w_cpp(
				x_levels_matrix,
				as.integer(private$num_levels_total),
				private$weights,
				private$p_best,
				private$prob_T
			)
			
			# Re-calculate final incremental counts state
			for (i in 1 : n){
				w_i = private$w[i]
				idx = x_levels_matrix[i, ]
				for (j in idx) private$counts[j, w_i + 1] = private$counts[j, w_i + 1] + 1
			}
		},

		get_bootstrap_indices = function() {
			sample_int_replace_cpp(private$t, private$t)
		}
	)
)
