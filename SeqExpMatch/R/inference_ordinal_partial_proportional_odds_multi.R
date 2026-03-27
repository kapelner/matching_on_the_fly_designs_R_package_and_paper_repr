#' Partial Proportional-Odds Inference for Ordinal Responses with Covariates
#'
#' @description
#' Fits a partial proportional-odds model for ordinal responses with treatment
#' and observed baseline covariates. Treatment remains parallel across
#' thresholds, while user-selected covariates may have threshold-specific
#' effects. When all covariates are parallel, the fit uses the package's fast
#' Rcpp proportional-odds solver.
#'
#' @param des_obj A completed \code{DesignSeqOneByOne} object whose response type is
#'   \code{"ordinal"}.
#' @param nonparallel Character vector of covariate names from \code{des_obj$
#'   get_X()} that should vary by threshold. Treatment is always kept parallel.
#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
#' @param verbose Whether to print progress messages.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalMultiPartialProportionalOddsRegr$
#'   new(
#'   seq_des,
#'   nonparallel = "x1",
#'   verbose = FALSE
#' )
#' infer
InferenceOrdinalMultiPartialProportionalOddsRegr = R6::R6Class(
	"InferenceOrdinalMultiPartialProportionalOddsRegr",
	inherit = InferenceOrdinalPartialProportionalOddsAbstract,
	public = list(
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
#' @param des_obj A completed \code{DesignSeqOneByOne} object whose response type is
#'   \code{"ordinal"}.
#' @param nonparallel Character vector of covariate names from \code{des_obj$
#'   get_X()} that should vary by threshold. Treatment is always kept parallel.
#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
#' @param verbose Whether to print progress messages.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalPartialProportionalOdds$
#'   new(
#'   seq_des,
#'   nonparallel = "x1",
#'   verbose = FALSE
#' )
#' infer
InferenceOrdinalPartialProportionalOdds = R6::R6Class(
	"InferenceOrdinalPartialProportionalOdds",
	inherit = InferenceOrdinalMultiPartialProportionalOddsRegr,
	public = list(
	)
)
