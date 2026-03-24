# A Sequential Design
#
# @description
# An abstract R6 Class encapsulating the data and functionality for a sequential experimental design.
# This class takes care of sequential assignments.
#
# @keywords internal
SeqDesign = R6::R6Class("SeqDesign",
	inherit = Design,
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

		supports_resampling = function(){
			TRUE
		},

		add_subject_to_experiment_and_assign = function(x_new){
			self$add_subject(x_new)
			w_t = self$assign_wt()
			if (private$fixed_sample){
				private$w[private$t] = w_t
			} else {
				private$w = c(private$w, w_t)
			}
			private$w[private$t]
		},

		assign_wt = function(){
			stop("Must be implemented by subclass.")
		},

		print_current_subject_assignment = function(){
			cat("Subject number", private$t, "is assigned to", ifelse(private$w[private$t] == 1, "TREATMENT", "CONTROL"), "via design", class(self)[1], "\n")
		}
	),

	private = list(
		redraw_w_according_to_design = function(){
			if (private$fixed_sample) {
				private$w = rep(NA_real_, private$n)
			} else {
				private$w = rep(NA_real_, private$t)
			}
			for (t in 1 : private$t){
				# This is a bit of a hack since private$t is at the end, 
				# but most subclasses override redraw_w_according_to_design anyway.
				private$w[t] = private$assign_wt() 
			}
		}
	)
)
