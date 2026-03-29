#' Survival Transformation Regression with Dependent Censoring
#'
#' Fits a lognormal transformation model for survival responses that jointly models
#' the event and censoring times. Dependence between the two transformed times is
#' represented with a Gaussian correlation parameter, allowing the censoring process
#' to be informative rather than independent. The treatment effect is reported on the
#' transformed log-time scale.
#'
#' @details
#' This implementation uses a parametric joint likelihood:
#' the transformed event time and transformed censoring time are each modeled as
#' Gaussian linear regressions with separate scale parameters, and their dependence is
#' captured by a correlation coefficient. For the univariate class, both regressions
#' include only the treatment indicator. The resulting treatment coefficient is an
#' adjusted log-time-ratio-style effect for the event process.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- InferenceSurvivalUniDepCensTransformRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalUniDepCensTransformRegr = R6::R6Class("InferenceSurvivalUniDepCensTransformRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMSummaryTable,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			X = matrix(private$w, ncol = 1)
			colnames(X) = "treatment"
			X
		},

		generate_mod = function(estimate_only = FALSE){
			mod = .fit_dep_cens_transform_model(
				y = private$y,
				dead = private$dead,
				Xmm = private$build_design_matrix()
			)
			if (!is.null(mod)) return(mod)
			stop("Dependent-censoring transformation model failed to converge.")
		}
	)
)
