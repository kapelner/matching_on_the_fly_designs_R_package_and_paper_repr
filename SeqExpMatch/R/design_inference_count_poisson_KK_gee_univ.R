#' Univariate GEE Inference for KK Designs with Count Response
#'
#' @description
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{geepack})
#' for count responses under a KK matching-on-the-fly design using only the
#' treatment indicator as a predictor (intercept + treatment). A Poisson working
#' model with a log link is used; the sandwich-robust standard errors remain valid
#' even under overdispersion. Matched pairs are treated as clusters (with exchangeable
#' correlation); reservoir subjects each form their own singleton cluster. The
#' treatment estimate is the log incidence-rate ratio (log-IRR).
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{geepack} before using this class.
#'
#' @export
DesignInferenceCountPoissonUnivKKGEE = R6::R6Class("DesignInferenceCountPoissonUnivKKGEE",
	inherit = DesignInferenceAbstractKKGEE,
	public = list(

		#' @description
		#' Initialize a univariate GEE inference object for a completed KK design
		#' with a count response.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose count response y is recorded.
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
		#' seq_des = SeqDesignKK14$new(n = 20, response_type = "count")
		#' for (i in 1 : 20){
		#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x = rnorm(1)))
		#' }
		#' seq_des$add_all_subject_responses(rpois(20, 3))
		#'
		#' seq_des_inf = DesignInferenceCountPoissonUnivKKGEE$new(seq_des)
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
		gee_response_type  = function() "count",
		gee_family         = function() poisson(link = "log"),
		gee_predictors_df  = function() data.frame(w = private$w)
	)
)
