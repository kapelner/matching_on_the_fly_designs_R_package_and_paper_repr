#' Atkinson's (1982) Covariate-Adjusted Biased Coin Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments.
#'
#' @export
SeqDesignAtkinson = R6::R6Class("SeqDesignAtkinson",
	inherit = SeqDesign,
	public = list(
		initialize = function(
						response_type = "continuous",
						prob_T = 0.5,
						include_is_missing_as_a_new_feature = TRUE,
						n = NULL,
						num_cores = 1,
						verbose = FALSE
					) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
		},

		assign_wt = function(){
			#if it's too early in the trial or if all the assignments are the same, then randomize
			if (private$t <= ncol(private$Xraw) + 2 + 1 | length(unique(private$t)) == 1){
				rbinom(1, 1, private$prob_T)
			} else {
				all_subject_data = private$compute_all_subject_data()
				tryCatch({
					atkinson_assign_weight_cpp(
						private$w[1 : (private$t - 1)],
						all_subject_data$X_prev,
						all_subject_data$xt_prev,
						all_subject_data$rank_prev,
						private$t
					)
				}, error = function(e){
					rbinom(1, 1, private$prob_T)
				})
			}
		}
	),
	private = list(
		redraw_w_according_to_design = function(){
			private$w[1:private$t] = atkinson_redraw_batch_cpp(
				as.matrix(private$X[1:private$t, , drop = FALSE]),
				private$t,
				ncol(private$Xraw),
				private$prob_T
			)
		}
	)
)
