#' Multivariate GEE Inference for KK Designs with Survival Response (no censoring)
#'
#' @description
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for uncensored survival (time-to-event) responses under a KK matching-on-the-fly
#' design using the treatment indicator and all recorded covariates as predictors.
#' A Gamma working model with a log link is used, which is appropriate for positive
#' continuous survival times. The sandwich-robust standard errors from GEE are valid
#' regardless of Gamma misspecification. Matched pairs are treated as clusters (with
#' exchangeable correlation); reservoir subjects each form their own singleton cluster.
#' The treatment estimate is the log ratio of mean survival times (log-MTR); inference
#' uses Wald Z-statistics based on the robust SE.
#'
#' @details
#' Censored observations are not supported; an error is raised at initialization if
#' any censoring is detected. This class requires the \pkg{geepack} package, which is
#' listed under \code{Suggests} and is not installed automatically with
#' \pkg{SeqExpMatch}. Install \pkg{geepack} manually
#' before using this class.
#'
#' @export
SeqDesignInferenceSurvivalMultiKKGammaGEE = R6::R6Class("SeqDesignInferenceSurvivalMultiKKGammaGEE",
	inherit = SeqDesignInferenceAbstractKKGEE,
	public = list(

		#' @description
		#' Initialize a multivariate GEE inference object for a completed KK design
		#' with an uncensored survival response.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose uncensored survival response y is recorded.
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
		#' seq_des = SeqDesignKK14$new(n = 20, response_type = "survival")
		#' for (i in 1 : 20){
		#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		#' }
		#' seq_des$add_all_subject_responses(rexp(20, rate = 0.5), rep(1, 20))
		#'
		#' seq_des_inf = SeqDesignInferenceSurvivalMultiKKGammaGEE$new(seq_des)
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
		gee_response_type = function() "survival",
		gee_family        = function() Gamma(link = "log")
	)
)
