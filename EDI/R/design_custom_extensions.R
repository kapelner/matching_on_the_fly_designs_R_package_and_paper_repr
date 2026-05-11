#' Internal base for user-defined fixed-design extensions
#'
#' \code{DesignFixedCustom} is intentionally not exported. Subclasses implement
#' \code{draw_assignments(r = 1)} and return an \code{n x r} 0/1 assignment
#' matrix. EDI handles subject storage, responses, and validation through
#' \code{DesignFixed}.
#'
#' @keywords internal
DesignFixedCustom = R6::R6Class("DesignFixedCustom",
	lock_objects = FALSE,
	inherit = DesignFixed,
	public = list(
		#' @description Draw assignments from the custom design.
		#' @param r Number of assignment vectors to draw.
		#' @return An \code{n x r} matrix of 0/1 assignments.
		draw_assignments = function(r = 1){
			stop("Custom fixed-design subclasses must implement public$draw_assignments(r).")
		},
		#' @description Internal redraw utility.
		#' @param r Number of assignment vectors.
		#' @return An \code{n x r} matrix.
		draw_ws_according_to_design = function(r = 100){
			if (should_run_asserts()) {
				assertCount(r, positive = TRUE)
			}
			w_mat = self$draw_assignments(r = as.integer(r))
			w_mat = as.matrix(w_mat)
			if (should_run_asserts()) {
				if (nrow(w_mat) != self$get_n() || ncol(w_mat) != as.integer(r)) {
					stop("draw_assignments(r) must return an n x r assignment matrix.", call. = FALSE)
				}
				if (any(!w_mat %in% c(0, 1))) {
					stop("draw_assignments(r) must return only 0/1 assignments.", call. = FALSE)
				}
			}
			w_mat
		}
	)
)
#' Internal base for user-defined sequential-design extensions
#'
#' \code{DesignCustomSequential} is intentionally not exported. Subclasses
#' implement \code{assignment_rule()} and return a scalar 0/1 assignment for the
#' current subject. EDI handles subject storage, responses, and redraws through
#' \code{DesignSeqOneByOne}.
#'
#' @keywords internal
DesignCustomSequential = R6::R6Class("DesignCustomSequential",
	lock_objects = FALSE,
	inherit = DesignSeqOneByOne,
	public = list(
		#' @description User-defined assignment rule.
		#' @return A binary treatment assignment.
		assignment_rule = function(){
			stop("Custom sequential-design subclasses must implement public$assignment_rule().")
		},
		#' @description Standard internal assignment entry point.
		#' @return A binary treatment assignment.
		assign_wt = function(){
			w_t = self$assignment_rule()
			if (should_run_asserts()) {
				assertChoice(w_t, c(0, 1))
			}
			as.numeric(w_t)
		}
	)
)
