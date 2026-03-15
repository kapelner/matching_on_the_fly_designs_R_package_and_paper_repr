#' Multivariate Quasi-Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a quasi-Poisson log-link regression for count responses using the treatment
#' indicator and all recorded covariates. The treatment effect is reported on the
#' log-rate scale and inference uses the model-based quasi-Poisson variance with an
#' estimated dispersion parameter.
#'
#' @export
SeqDesignInferenceCountMultiQuasiPoissonRegr = R6::R6Class("SeqDesignInferenceCountMultiQuasiPoissonRegr",
	inherit = SeqDesignInferenceCountUnivQuasiPoissonRegr,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param	verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			Xmm = private$create_design_matrix()
			colnames(Xmm) = c("(Intercept)", "treatment", if (ncol(Xmm) > 2L) paste0("x", seq_len(ncol(Xmm) - 2L)) else NULL)
			Xmm
		}
	)
)
