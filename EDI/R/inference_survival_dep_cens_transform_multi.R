#' Multivariate Survival Transformation Regression with Dependent Censoring
#'
#' Fits a lognormal transformation model for survival responses that jointly models
#' the event and censoring times using the treatment indicator and all recorded
#' covariates. Dependence between the two transformed times is represented with a
#' Gaussian correlation parameter, allowing the censoring process to be informative
#' rather than independent. The treatment effect is reported on the transformed
#' log-time scale.
#'
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
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
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
	lock_objects = FALSE,
	inherit = InferenceSurvivalUniDepCensTransformRegr,
	public = list(
		#' @description
		#' Initialize the Inference object.
		#'
		#' @param des_obj The design object.
		#' @param verbose If TRUE, print additional information.
		initialize = function(des_obj, verbose = FALSE) {
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose)
		}

	),

	private = list(
		build_design_matrix_candidates = function(){
			X_cov_orig = private$get_X()
			if (is.null(X_cov_orig) || ncol(X_cov_orig) == 0L){
				X = matrix(private$w, ncol = 1)
				colnames(X) = "treatment"
				return(list(X))
			}

			if (!private$harden) {
				X_full = cbind(treatment = private$w, X_cov_orig)
				return(list(X_full))
			}

			thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			candidates = list()
			keys = character()
			base_colnames = colnames(X_cov_orig)
			if (is.null(base_colnames)) base_colnames = rep("x", ncol(X_cov_orig))

			# This model has 2*(p+1) + 3 parameters. To keep optimization feasible, 
			# we cap p such that parameters < n / 4 (approx).
			max_allowed_p = max(1L, floor((private$n - 5) / 4))

			for (thresh in thresholds){
				X_cov = if (is.finite(thresh)) drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M else X_cov_orig
				
				# If we still have too many columns, we must drop more based on correlation
				# but here we just skip this threshold if it's over the cap.
				if (ncol(X_cov) > max_allowed_p) next

				X_full = cbind(treatment = private$w, X_cov)
				reduced = qr_reduce_preserve_cols_cpp(as.matrix(X_full), required_cols = 1L)
				X_red = reduced$X_reduced
				keep_idx = as.integer(reduced$keep)
				if (length(keep_idx) == 0L) next
				keep_idx0 = keep_idx
				colnames(X_red) = colnames(X_full)[keep_idx0]
				key = paste(colnames(X_red), collapse = "|")
				if (!(key %in% keys)){
					candidates[[length(candidates) + 1L]] = X_red
					keys = c(keys, key)
				}
			}
			if (length(candidates) == 0L){
				X = matrix(private$w, ncol = 1)
				colnames(X) = "treatment"
				candidates[[1L]] = X
			}
			candidates
		},

		build_design_matrix = function(){
			private$build_design_matrix_candidates()[[1L]]
		}
	)
	)
