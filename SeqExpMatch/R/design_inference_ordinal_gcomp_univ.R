#' Univariate G-Computation for Ordinal Responses
#'
#' @description
#' Fits a proportional odds model with only the treatment indicator, then
#' standardizes to estimate the marginal mean difference.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- DesignInferenceOrdinalUniGCompMeanDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer$
#'   compute_treatment_estimate()
DesignInferenceOrdinalUniGCompMeanDiff = R6::R6Class("DesignInferenceOrdinalUniGCompMeanDiff",
	inherit = DesignInferenceOrdinalGCompAbstract,
	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		}
	)
)
