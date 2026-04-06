#' A Fixed Design
#'
#' An abstract R6 Class encapsulating the data and functionality for a fixed experimental design.
#' This class takes care of whole-experiment randomization.
#'
#' @keywords internal
#' @export
FixedDesign = R6::R6Class("FixedDesign",
	inherit = Design,
	public = list(
		#' @description
		#' Initialize a fixed experimental design
		#'
		#' @param response_type   "continuous", "incidence", "proportion", "count", "survival", or
		#'   "ordinal".
		#' @param	prob_T	Probability of treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param verbose A flag for verbosity.
		#'
		#' @return	A new `FixedDesign` object
		initialize = function(
				response_type,
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				verbose = FALSE
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
		},

		#' @description
		#' Assign treatment to all subjects in the fixed experiment.
		assign_w_to_all_subjects = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		#' @description
		#' Add all subjects' covariates to a fixed design at once.
		#'
		#' @param X_all A data frame containing the full covariate matrix.
		#' @return Invisibly returns the design object.
		add_all_subjects_to_experiment = function(X_all){
			assertClass(X_all, "data.frame")
			self$assert_fixed_sample()
			n_all = nrow(X_all)
			n_expected = self$get_n()
			if (n_all != n_expected){
				stop("X_all must have exactly ", n_expected, " rows for this fixed design.")
			}
			if (private$t > 0L || nrow(private$Xraw) > 0L){
				stop("Subjects have already been added to this design.")
			}
			private$Xraw = copy(as.data.table(X_all))
			private$p_raw_t = ncol(private$Xraw)
			private$t = n_all
			private$covariate_impute_if_necessary_and_then_create_model_matrix()
			invisible(self)
		},

		#' @description
		#' Add all subject responses for a fixed design.
		#'
		#' @param ys The responses as a numeric vector.
		#' @param deads The binary vector indicating if dead/censored.
		add_all_subject_responses = function(ys, deads = NULL){
			self$assert_fixed_sample()
			if (is.null(deads)){
				deads = rep(1, private$t)
			}
			if (private$response_type == "ordinal" && is.factor(ys)){
				assertFactor(ys, len = private$t, ordered = TRUE, any.missing = FALSE)
			} else {
				assertNumeric(ys, len = private$t)
			}
			assertNumeric(deads, len = private$t)

			for (t in 1 : private$t){
				self$add_one_subject_response(t, ys[t], deads[t])
			}
		},

		#' @description
		#' Overwrite all subject assignments for a fixed design.
		#'
		#' @param w The binary responses.
		overwrite_all_subject_assignments = function(w){
			assertIntegerish(w, lower = 0, upper = 1, any.missing = FALSE, len = private$t)
			private$w = w
		},

		#' @description
		#' Check if the design supports resampling.
		#'
		#' @return 	TRUE if supported.
		supports_resampling = function(){
			class(self)[1] != "FixedDesign"
		},


		#' @description
		#' Draw multiple treatment assignment vectors.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			stop("Must be implemented by subclass.")
		}
	)
)
