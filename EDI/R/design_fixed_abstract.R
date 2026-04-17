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
			n_all = nrow(X_all)
			n_expected = self$get_n()
			if (should_run_asserts()) {
				assertClass(X_all, "data.frame")
				self$assert_fixed_sample()
				if (n_all != n_expected){
					stop("X_all must have exactly ", n_expected, " rows for this fixed design.")
				}
				if (private$t > 0L || nrow(private$Xraw) > 0L){
					stop("Subjects have already been added to this design.")
				}
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
			if (is.null(deads)){
				deads = rep(1, private$t)
			}
			if (should_run_asserts()) {
				self$assert_fixed_sample()
				if (private$response_type == "ordinal" && is.factor(ys)){
					assertFactor(ys, len = private$t, ordered = TRUE, any.missing = FALSE)
				} else {
					assertNumeric(ys, len = private$t)
				}
				assertNumeric(deads, len = private$t)
				
				if (private$response_type != "survival" && any(deads == 0)){
					stop("censored observations are only available for survival response types")
				}
			}
			
			if (private$response_type == "ordinal" && is.factor(ys)){
				ys = as.integer(ys)
			}

			private$y = as.numeric(ys)
			private$dead = as.numeric(deads)
			private$y_i_t_i = as.list(seq_len(private$t))
		},

		#' @description
		#' Overwrite all subject assignments for a fixed design.
		#'
		#' @param w The binary responses.
		overwrite_all_subject_assignments = function(w){
			if (should_run_asserts()) {
				assertIntegerish(w, lower = 0, upper = 1, any.missing = FALSE, len = private$t)
			}
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
			if (should_run_asserts()) {
				stop("Must be implemented by subclass.")
			}
		}
	),
	private = list(
		strata_cols = NULL,
		num_bins_for_continuous_covariate = NULL,
		B_preferred = NULL,

		get_strata_keys = function(){
			n = private$t
			if (n == 0) return(character(0))
			strata_cols = if (is.null(private$strata_cols)) names(private$Xraw) else private$strata_cols

			col_to_str = function(col) {
				vec = private$Xraw[[col]]
				if (is.numeric(vec)) {
					probs = seq(0, 1, length.out = private$num_bins_for_continuous_covariate + 1)
					breaks = unique(stats::quantile(vec, probs = probs, na.rm = TRUE))
					s = if (length(breaks) > 1) as.character(cut(vec, breaks = breaks, include.lowest = TRUE)) else as.character(vec)
				} else {
					s = as.character(vec)
				}
				s[is.na(s)] = "NA"
				s
			}

			append_key = function(keys, col_str) {
				if (all(nchar(keys) == 0L)) col_str else paste(keys, col_str, sep = "|")
			}

			keys = rep("", n)

			if (!is.null(private$B_preferred)) {
				target = if (is.na(private$B_preferred)) floor(sqrt(n)) else as.integer(private$B_preferred)
				for (col in strata_cols) {
					new_keys = append_key(keys, col_to_str(col))
					if (length(unique(new_keys)) <= target) {
						keys = new_keys
					}
				}
			} else {
				for (col in strata_cols) {
					keys = append_key(keys, col_to_str(col))
				}
			}

			num_blocks = length(unique(keys))
			if (should_run_asserts()) {
				if (num_blocks > n) {
					stop("Number of blocks (", num_blocks, ") exceeds sample size (", n, "). Reduce the number of strata columns or use fewer bins.")
				}
			}
			keys
		}
	)
)
