#' Univariate GLMM Inference for KK Designs with Ordinal Response
#'
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
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUnivKKGLMM$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalUnivKKGLMM = R6::R6Class("InferenceOrdinalUnivKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKOrdinalCLMM,
	public = list(

	),
	private = list(
		clmm_link = function() "logit"
	)
)
