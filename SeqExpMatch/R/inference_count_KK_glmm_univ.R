#' Univariate GLMM Inference for KK Designs with Count Response
#'
#' @description
#' Fits a Generalized Linear Mixed Model (GLMM) via \code{glmmTMB::glmmTMB} for count
#' responses under a KK matching-on-the-fly design using only the treatment indicator
#' as a fixed-effect predictor (intercept + treatment). A Poisson log-link working
#' model is used. The matched-pair strata enter the model as a subject-level random
#' intercept \code{(1 | group_id)}, which accounts for within-pair correlation.
#' Reservoir subjects each receive a unique singleton group. The treatment estimate
#' is the conditional log incidence-rate ratio (log-IRR); inference uses the Wald
#' Z-statistic from the GLMM fixed-effects table.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed under \code{Suggests}
#' and is not installed automatically with \pkg{SeqExpMatch}. Install it manually with
#' \code{install.packages("glmmTMB")} before using this class.
#'
#' @export
SeqDesignInferenceCountUnivKKGLMM = R6::R6Class("SeqDesignInferenceCountUnivKKGLMM",
	inherit = SeqDesignInferenceAbstractKKGLMM,
	public = list(

		#' @description
		#' Initialize a univariate GLMM inference object for a completed KK design
		#' with a count response.
		#' @param seq_des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose count response y is recorded.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param verbose			Whether to print progress messages. Default is \code{FALSE}.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignKK14$new(n = 20, response_type = "count")
		#' for (i in 1 : 20){
		#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		#' }
		#' seq_des$add_all_subject_responses(rpois(20, 3))
		#'
		#' seq_des_inf = SeqDesignInferenceCountUnivKKGLMM$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	private = list(
		glmm_response_type  = function() "count",
		glmm_family         = function() poisson(link = "log"),
		glmm_predictors_df  = function() data.frame(w = private$w)
	)
)
