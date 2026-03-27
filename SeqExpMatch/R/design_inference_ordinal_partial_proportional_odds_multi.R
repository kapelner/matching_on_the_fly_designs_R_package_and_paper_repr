#' Partial Proportional-Odds Inference for Ordinal Responses with Covariates
#'
#' @description
#' Fits a partial proportional-odds model for ordinal responses with treatment
#' and observed baseline covariates. Treatment remains parallel across
#' thresholds, while user-selected covariates may have threshold-specific
#' effects. When all covariates are parallel, the fit uses the package's fast
#' Rcpp proportional-odds solver.
#'
#' @param des_obj A completed \code{SeqDesign} object whose response type is
#'   \code{"ordinal"}.
#' @param nonparallel Character vector of covariate names from \code{des_obj$
#'   get_X()} that should vary by threshold. Treatment is always kept parallel.
#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
#' @param verbose Whether to print progress messages.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- DesignInferenceOrdinalMultiPartialProportionalOddsRegr$
#'   new(
#'   seq_des,
#'   nonparallel = "x1",
#'   verbose = FALSE
#' )
#' infer
DesignInferenceOrdinalMultiPartialProportionalOddsRegr = R6::R6Class(
	"DesignInferenceOrdinalMultiPartialProportionalOddsRegr",
	inherit = DesignInferenceOrdinalPartialProportionalOddsAbstract,
	public = list(
		#' @description
		#' Initialize a multivariate partial proportional-odds inference object.
		#' @param des_obj A completed \code{SeqDesign} object with an ordinal
		#'   response.
		#' @param nonparallel Covariate names that should vary across thresholds.
		#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj,
						nonparallel = character(0),
						num_cores = 1,
						verbose = FALSE){
			super$initialize(
				des_obj = des_obj,
				nonparallel = nonparallel,
				num_cores = num_cores,
				verbose = verbose
			)
		}
	),
	private = list(
		ppo_covariate_matrix = function(){
			X_raw = private$get_X()
			if (ncol(X_raw) == 0){
				return(matrix(0, nrow = private$n, ncol = 0))
			}
			X_full = cbind(treatment = private$w, as.matrix(X_raw))
			qr_X = qr(X_full)
			if (qr_X$rank < ncol(X_full)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(1L %in% keep)) keep[qr_X$rank] = 1L
				keep = sort(unique(keep))
				X_full = X_full[, keep, drop = FALSE]
			}
			X_cov = X_full[, setdiff(colnames(X_full), "treatment"), drop = FALSE]
			if (ncol(X_cov) == 0){
				return(matrix(0, nrow = private$n, ncol = 0))
			}
			X_cov
		}
	)
)

#' Partial Proportional-Odds Inference for Ordinal Responses
#'
#' @description
#' Backward-compatible alias for the multivariate partial proportional-odds
#' estimator.
#'
#' @param des_obj A completed \code{SeqDesign} object whose response type is
#'   \code{"ordinal"}.
#' @param nonparallel Character vector of covariate names from \code{des_obj$
#'   get_X()} that should vary by threshold. Treatment is always kept parallel.
#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
#' @param verbose Whether to print progress messages.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- DesignInferenceOrdinalPartialProportionalOdds$
#'   new(
#'   seq_des,
#'   nonparallel = "x1",
#'   verbose = FALSE
#' )
#' infer
DesignInferenceOrdinalPartialProportionalOdds = R6::R6Class(
	"DesignInferenceOrdinalPartialProportionalOdds",
	inherit = DesignInferenceOrdinalMultiPartialProportionalOddsRegr,
	public = list(
		#' @description
		#' Initialize the backward-compatible multivariate PPO estimator.
		#' @param des_obj A completed \code{SeqDesign} object with an ordinal
		#'   response.
		#' @param nonparallel Covariate names that should vary across thresholds.
		#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj,
						nonparallel = character(0),
						num_cores = 1,
						verbose = FALSE){
			super$initialize(
				des_obj = des_obj,
				nonparallel = nonparallel,
				num_cores = num_cores,
				verbose = verbose
			)
		}
	)
)
