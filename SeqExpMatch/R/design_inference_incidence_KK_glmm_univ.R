#' Univariate GLMM Inference for KK Designs with Binary Response
#'
#' @description
#' Fits a Generalized Linear Mixed Model (GLMM) using the \pkg{glmmTMB} fitter for binary
#' (incidence) responses under a KK matching-on-the-fly design using only the
#' treatment indicator as a fixed-effect predictor (intercept + treatment). The
#' matched-pair strata enter the model as a subject-level random intercept
#' \code{(1 | group_id)}, which accounts for within-pair correlation. Reservoir
#' subjects each receive a unique singleton group. The treatment estimate is the
#' conditional (subject-specific) log-odds ratio; inference uses the Wald Z-statistic
#' from the GLMM fixed-effects table.
#'
#' @details
#' This class requires the \pkg{glmmTMB} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{glmmTMB} before using this class.
#'
#' @export
DesignInferenceIncidUnivKKGLMM = R6::R6Class("DesignInferenceIncidUnivKKGLMM",
	inherit = DesignInferenceAbstractKKGLMM,
	public = list(

		#' @description
		#' Initialize a univariate GLMM inference object for a completed KK design
		#' with a binary (incidence) response.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose binary response y is recorded.
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
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignKK14$new(n = 20, response_type = "incidence")
		#' for (i in 1 : 20){
		#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		#' }
		#' seq_des$add_all_subject_responses(rbinom(20, 1, 0.5))
		#'
		#' seq_des_inf = DesignInferenceIncidUnivKKGLMM$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			super$compute_asymp_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_asymp_two_sided_pval_for_treatment_effect(delta = delta)
		}
	),
	private = list(
		glmm_response_type  = function() "incidence",
		glmm_family         = function() binomial(link = "logit"),
		glmm_predictors_df  = function() data.frame(w = private$w)
	)
)
