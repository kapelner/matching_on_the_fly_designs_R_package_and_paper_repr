#' Sequential One-by-One Experimental Design
#'
#' Abstract R6 class encapsulating data and functionality for a sequential one-by-
#' one experimental design.
#'
#' @keywords internal
#' @export
DesignSeqOneByOne = R6::R6Class("DesignSeqOneByOne",
	inherit = Design,
	public = list(
		#' @description
		#' Initialize a sequential one-by-one design.
		#' @param response_type The data type of response values.
		#' @param prob_T The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature If missing data is present, include a dummy
		#'   variable for it.
		#' @param n The sample size (if fixed).
		#' @param num_cores Number of CPU cores for parallel processing.
		#' @param verbose Whether to print progress messages.
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

		#' @description
		#' Checks if the design supports resampling.
		#' @return Always TRUE for sequential designs.
		supports_resampling = function(){
			TRUE
		},

		#' @description
		#' Adds a subject and assigns treatment.
		#' @param x_new A data frame with one row representing the new subject's covariates.
		#' @return The treatment assignment (0 or 1).
		add_subject_to_experiment_and_assign = function(x_new){
			self$add_subject(x_new)
			w_t = self$assign_wt()
			if (private$fixed_sample){
				private$w[private$t] = w_t
			} else {
				if (length(private$w) < private$t) {
					private$w = c(private$w, w_t)
				} else {
					private$w[private$t] = w_t
				}
			}
			private$w[private$t]
		},

		#' @description
		#' Assigns treatment to the current subject.
		#' @return The treatment assignment (0 or 1).
		assign_wt = function(){
			stop("Must be implemented by subclass.")
		},

		#' @description
		#' Draws treatment vectors according to the design.
		#' @param r Number of vectors to draw.
		#' @return A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			w_mat = matrix(NA_real_, nrow = private$t, ncol = r)
			for (j in 1 : r){
				if (private$fixed_sample) {
					private$w = rep(NA_real_, private$n)
				} else {
					private$w = rep(NA_real_, private$t)
				}
				for (t in 1 : private$t){
					private$w[t] = self$assign_wt()
				}
				w_mat[, j] = private$w[1:private$t]
			}
			w_mat
		},

		#' @description
		#' Prints the current subject's assignment.
		print_current_subject_assignment = function(){
			cat("Subject number", private$t, "is assigned to", ifelse(private$w[private$t] == 1, "TREATMENT", "CONTROL"), "via design", class(self)[1], "\n")
		}
	),

	private = list()
)
