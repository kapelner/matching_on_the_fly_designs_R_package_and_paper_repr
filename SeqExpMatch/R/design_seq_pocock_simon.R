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
		}
	),

	private = list(
		strata_cols = NULL,
		weights = NULL,
		p_best = NULL,
		factor_levels = list(), # list of levels for each strata_col
		num_levels_total = NULL,
		level_offsets = NULL,   # offset in the counts matrix for each covariate

		assign_wt = function(){
			# Ensure we have factor levels precomputed
			private$ensure_factor_metadata()
			
			# Current subject data
			x_new = private$Xraw[private$t, ]
			subject_levels_idx = private$get_subject_levels_idx(x_new)
			
			# Current counts
			counts = private$compute_current_counts()
			
			pocock_simon_assign_cpp(
				counts, 
				as.integer(subject_levels_idx), 
				private$weights, 
				private$p_best, 
				private$prob_T
			)
		},

		ensure_factor_metadata = function(){
			if (length(private$factor_levels) == 0){
				offset = 0
				private$level_offsets = numeric(length(private$strata_cols))
				
				for (i in 1 : length(private$strata_cols)){
					col = private$strata_cols[i]
					vals = private$Xraw[[col]]
					lvls = levels(as.factor(vals))
					if (length(lvls) == 0) lvls = "NA" # Fallback
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
				
				lvl_idx = which(private$factor_levels[[col]] == val)
				if (length(lvl_idx) == 0) lvl_idx = 1 # Should not happen if ensure_factor_metadata works
				
				indices[i] = private$level_offsets[i] + lvl_idx
			}
			indices
		},

		compute_current_counts = function(){
			# Matrix of num_levels_total x 2
			counts = matrix(0, nrow = private$num_levels_total, ncol = 2)
			if (private$t <= 1) return(counts)
			
			# We only count subjects 1 to t-1
			for (i in 1 : (private$t - 1)){
				w_i = private$w[i]
				if (is.na(w_i)) next
				
				idx = private$get_subject_levels_idx(private$Xraw[i, ])
				for (j in idx){
					counts[j, w_i + 1] = counts[j, w_i + 1] + 1
				}
			}
			counts
		},

		redraw_w_according_to_design = function(){
			private$ensure_factor_metadata()
			n = private$t
			x_levels_matrix = matrix(NA_integer_, nrow = n, ncol = length(private$strata_cols))
			for (i in 1 : n){
				x_levels_matrix[i, ] = private$get_subject_levels_idx(private$Xraw[i, ])
			}
			
			private$w[1:n] = pocock_simon_redraw_w_cpp(
				x_levels_matrix,
				as.integer(private$num_levels_total),
				private$weights,
				private$p_best,
				private$prob_T
			)
		},

		get_bootstrap_indices = function() {
			sample_int_replace_cpp(private$t, private$t)
		}
	)
)
