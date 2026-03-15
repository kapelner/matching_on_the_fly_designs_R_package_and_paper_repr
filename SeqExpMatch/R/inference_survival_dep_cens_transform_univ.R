#' Survival Transformation Regression with Dependent Censoring
#'
#' @description
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
SeqDesignInferenceSurvivalUniDepCensTransformRegr = R6::R6Class("SeqDesignInferenceSurvivalUniDepCensTransformRegr",
	inherit = SeqDesignInferenceMLEorKMSummaryTable,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize the sampling during randomization-based inference and bootstrap resampling.
		#' @param verbose A flag indicating whether messages should be displayed to the user. Default is \code{FALSE}.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			assertResponseType(seq_des_obj$get_response_type(), "survival")
		}
	),

	private = list(
		build_design_matrix = function(){
			X = matrix(private$w, ncol = 1)
			colnames(X) = "treatment"
			X
		},

		generate_mod = function(){
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
