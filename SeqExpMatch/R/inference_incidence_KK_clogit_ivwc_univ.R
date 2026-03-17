#' Univariate Conditional Logistic Regression Inference for KK Designs with Binary Response
#'
#' @description
#' Fits a compound estimator for KK matching-on-the-fly designs with binary (incidence)
#' responses using only the treatment indicator (no additional covariates). For matched
#' pairs, a conditional logistic regression model is used (via the internal
#' \code{clogit_helper}, which processes discordant pairs and fits logistic regression
#' on within-pair signed differences). For reservoir subjects, a standard logistic
#' regression is used. The two estimates are combined via a variance-weighted linear
#' combination, analogous to \code{SeqDesignInferenceContinMultOLSKK}. If clogit fails
#' (e.g. no discordant pairs), the estimator falls back to logistic regression on the
#' reservoir only.
#'
#' @export
SeqDesignInferenceIncidUnivKKClogitIVWC = R6::R6Class("SeqDesignInferenceIncidUnivKKClogitIVWC",
	inherit = SeqDesignInferenceAbstractKKClogitIVWC,
	public = list(

		#' @description
		#' Initialize a univariate conditional logistic regression inference object for a
		#' completed KK design with a binary (incidence) response.
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
		#' seq_des_inf = SeqDesignInferenceIncidUnivKKClogitIVWC$new(seq_des)
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
		include_covariates = function() FALSE
	)
)
