#' Conditional Logistic Plus GLMM IVWC Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKClogitPlusGLMMOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceIncidKKClogitPlusGLMMIVWC = R6::R6Class("InferenceIncidKKClogitPlusGLMMIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitPlusGLMM,
	public = list(),
	private = list(
		combine_reservoir_into_glmm = function() FALSE
	)
)

#' Conditional Logistic Plus GLMM Combined-Likelihood Inference for KK Designs
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The treatment effect is estimated by the conditional-logistic
#' component; the GLMM component includes only the intercept and covariates.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKClogitPlusGLMMOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKClogitPlusGLMMOneLik = R6::R6Class("InferenceIncidKKClogitPlusGLMMOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitPlusGLMM,
	public = list(),
	private = list(
		combine_reservoir_into_glmm = function() TRUE
	)
)
