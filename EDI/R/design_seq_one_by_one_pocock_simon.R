#' Pocock & Simon's Minimization Sequential Design
#'
#' An R6 Class encapsulating the data and functionality for a Pocock & Simon
#' sequential experimental design.
#' This design minimizes the imbalance across treatments for multiple covariates.
#'
#' @export
DesignSeqOneByOnePocockSimon = R6::R6Class("DesignSeqOneByOnePocockSimon",
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description
		#' Initialize a Pocock & Simon sequential experimental design
		#'
		#' @param strata_cols     The names of the covariates to be used for minimization. These
		#'   must be factor or categorical variables.
		#' @param weights 		A numeric vector of weights for each covariate. Defaults to 1 for all.
		#' @param p_best          The probability of assigning the treatment that minimizes the
		#'   imbalance. Defaults to 0.8.
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `DesignSeqOneByOnePocockSimon` object
		#'
		initialize = function(
				strata_cols,
				weights = NULL,
				p_best = 0.8,
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				
				verbose = FALSE
			) {
			assertCharacter(strata_cols, min.len = 1)
			assertNumeric(p_best, lower = 0.5, upper = 1)
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
			
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

		#' @description
		#' Assign the next subject to a treatment group using minimization.
		#'
		#' @return 	The treatment assignment (0 or 1)
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
		},

		#' @description
		#' Draw multiple treatment assignment vectors according to Pocock & Simon design.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			private$ensure_factor_metadata()
			n = self$get_n()
			x_levels_matrix = matrix(NA_integer_, nrow = n, ncol = length(private$strata_cols))
			for (i in 1 : n){
				x_levels_matrix[i, ] = private$get_subject_levels_idx(private$Xraw[i, ])
			}

			generate_permutations_pocock_simon_cpp(
				x_levels_matrix,
				as.integer(private$num_levels_total),
				private$weights,
				private$p_best,
				private$prob_T,
				as.integer(r)
			)$w_mat
		}

	),

	private = list(
		p_best = NULL,
		weights = NULL,
		counts = NULL,
		draw_bootstrap_indices = function() {
			list(i_b = sample_int_replace_cpp(private$t, private$t), m_vec_b = NULL)
		}
	)
)
