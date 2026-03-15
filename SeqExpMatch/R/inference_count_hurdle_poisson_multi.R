#' Multivariate Hurdle Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a hurdle Poisson regression for count responses using the treatment
#' indicator and all recorded covariates in both the count and hurdle components.
#' The reported treatment effect is the treatment coefficient from the conditional
#' truncated-Poisson count component, on the log-rate scale.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "count", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- SeqDesignInferenceCountMultiHurdlePoissonRegr$new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceCountMultiHurdlePoissonRegr = R6::R6Class("SeqDesignInferenceCountMultiHurdlePoissonRegr",
	inherit = SeqDesignInferenceCountUnivHurdlePoissonRegr,
	private = list(
		predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		}
	)
)
