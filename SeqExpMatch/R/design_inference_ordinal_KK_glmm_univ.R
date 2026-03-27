#' Univariate GLMM Inference for KK Designs with Ordinal Response
#'
#' @description
#' Fits a Generalized Linear Mixed Model (GLMM) (via \code{ordinal::clmm})
#' for ordinal responses under a KK matching-on-the-fly design using only
#' the treatment indicator as a fixed effect predictor (intercept + treatment)
#' and a random intercept for each matched pair (cluster). Reservoir subjects
#' each form their own singleton cluster. Inference is based on the Wald test.
#'
#' @details
#' This class requires the \pkg{ordinal} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{ordinal} before using this class.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$
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
#' infer <- DesignInferenceOrdinalUnivKKGLMM$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceOrdinalUnivKKGLMM = R6::R6Class("DesignInferenceOrdinalUnivKKGLMM",
	inherit = DesignInferenceAbstractKKOrdinalCLMM,
	public = list(

		#' @description
		#' Initialize a univariate GLMM inference object for a completed KK design
		#' with an ordinal response.
		#' @param	des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose ordinal response y is recorded.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param	verbose			Whether to print progress messages. Default is \code{FALSE}.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		}
	),
	private = list(
		clmm_link = function() "logit"
	)
)
