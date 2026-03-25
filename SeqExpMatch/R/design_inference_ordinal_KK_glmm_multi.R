#' Multivariate GLMM Inference for KK Designs with Ordinal Response
#'
#' @description
#' Fits a cumulative-link mixed model via \code{ordinal::clmm} for ordinal
#' responses under a KK matching-on-the-fly design using the treatment indicator
#' and all recorded covariates as fixed effects. The matched-pair strata enter
#' the model as a random intercept \code{(1 | group_id)}. Reservoir subjects each
#' receive a unique singleton group. Inference is based on the Wald Z-statistic.
#'
#' @details
#' This class requires the \pkg{ordinal} package, which is listed under
#' \code{Suggests} and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{ordinal} manually before using this class.
#'
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
#' infer <- DesignInferenceOrdinalMultiKKGLMM$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceOrdinalMultiKKGLMM = R6::R6Class("DesignInferenceOrdinalMultiKKGLMM",
	inherit = DesignInferenceAbstractKKOrdinalCLMM,
	public = list(
		#' @description
		#' Initialize a multivariate ordinal CLMM inference object for a completed
		#' KK design.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design) whose entire
		#' 							n subjects are assigned and whose ordinal response y is
		#' 							recorded.
		#' @param	num_cores			The number of CPU cores to use to parallelize the
		#' 							sampling during randomization-based inference and bootstrap
		#' 							resampling.
		#' @param	verbose			Whether to print progress messages. Default is
		#' 							\code{FALSE}.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),
	private = list(
		clmm_link = function() "logit"
	)
)
