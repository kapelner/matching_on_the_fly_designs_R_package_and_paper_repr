#' Multivariate Survival Transformation Regression with Dependent Censoring
#'
#' @description
#' Fits a lognormal transformation model for survival responses that jointly models
#' the event and censoring times using the treatment indicator and all recorded
#' covariates. Dependence between the two transformed times is represented with a
#' Gaussian correlation parameter, allowing the censoring process to be informative
#' rather than independent. The treatment effect is reported on the transformed
#' log-time scale.
#'
#' @inherit InferenceRand methods
#' @inherit InferenceBoot methods
#' @inherit InferenceAsymp methods
#' @inherit InferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
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
#' infer <- InferenceSurvivalMultiDepCensTransformRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalMultiDepCensTransformRegr = R6::R6Class("InferenceSurvivalMultiDepCensTransformRegr",
	inherit = InferenceSurvivalUniDepCensTransformRegr,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj A DesignSeqOneByOne object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{FALSE}.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
		}
	),

	private = list(
		build_design_matrix = function(){
			X_cov_orig = private$get_X()
			if (ncol(X_cov_orig) == 0L){
				X = matrix(private$w, ncol = 1)
				colnames(X) = "treatment"
				return(X)
			}

			thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			best_X = NULL
			for (thresh in thresholds){
				X_cov = if (is.finite(thresh)) drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M else X_cov_orig
				X = cbind(X_cov, treatment = private$w)
				qr_X = qr(X)
				if (qr_X$rank < ncol(X)){
					keep = qr_X$pivot[seq_len(qr_X$rank)]
					treat_col = which(colnames(X) == "treatment")
					if (!(treat_col %in% keep)) keep = c(keep, treat_col)
					keep = sort(unique(keep))
					X = X[, keep, drop = FALSE]
				}
				if ("treatment" %in% colnames(X)){
					best_X = X
					break
				}
			}
			if (is.null(best_X)){
				best_X = matrix(private$w, ncol = 1)
				colnames(best_X) = "treatment"
			}
			best_X
		}
	)
)
