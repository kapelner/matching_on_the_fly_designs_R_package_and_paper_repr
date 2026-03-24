#' Univariate Stratified Cox PH Regression Inference for Survival Responses
#'
#' @description
#' Fits an all-subject stratified Cox proportional hazards model for survival
#' responses using only the treatment indicator in the linear predictor. The strata
#' are chosen automatically from low-cardinality observed covariates. If no suitable
#' stratification covariates are available, the method falls back to the standard
#' univariate Cox PH model.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- SeqDesignInferenceSurvivalUniStratCoxPHRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
SeqDesignInferenceSurvivalUniStratCoxPHRegr = R6::R6Class("SeqDesignInferenceSurvivalUniStratCoxPHRegr",
	inherit = SeqDesignInferenceSurvivalStratCoxPHAbstract,
	private = list(
		include_covariates = function() FALSE
	)
)
