#' Univariate GEE Inference for KK Designs with Binary Response
#'
#' @description
#' Fits a Generalized Estimating Equations (GEE) model (via \code{geepack::geeglm})
#' for binary (incidence) responses under a KK matching-on-the-fly design using only
#' the treatment indicator as a predictor (intercept + treatment). Matched pairs are
#' treated as clusters (with exchangeable correlation structure); reservoir subjects
#' each form their own singleton cluster. Unlike
#' \code{SeqDesignInferenceAbstractKKClogit}, all subjects (matched and reservoir) are
#' included. Inference is based on sandwich-robust standard errors, so the test
#' statistic is Z-distributed.
#'
#' @details
#' This class requires the \pkg{geepack} package, which is listed under \code{Suggests}
#' and is not installed automatically with \pkg{SeqExpMatch}. Install it manually with
#' \code{install.packages("geepack")} before using this class.
#'
#' @export
SeqDesignInferenceIncidUnivKKGEE = R6::R6Class("SeqDesignInferenceIncidUnivKKGEE",
	inherit = SeqDesignInferenceAbstractKKGEE,
	public = list(

		#' @description
		#' Initialize a univariate GEE inference object for a completed KK design
		#' with a binary (incidence) response.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose binary response y is recorded.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
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
		#' seq_des_inf = SeqDesignInferenceIncidUnivKKGEE$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	private = list(
		gee_response_type = function() "incidence",
		gee_family        = function() binomial(link = "logit"),
		gee_predictors_df = function() data.frame(w = private$w)
	)
)
