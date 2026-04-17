#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
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
#'   response_type = "proportion",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0.10, 0.25, 0.20, 0.40, 0.35, 0.55, 0.60, 0.75))
#' infer <- InferencePropUniBetaRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferencePropUniBetaRegr = R6::R6Class("InferencePropUniBetaRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
				assertNumeric(private$y, any.missing = FALSE)
			}
		}
		),

	private = list(
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			
			Xmm = if (length(Xmm_cols) == 0L){
				# Univariate case
				cbind(1, private$w)
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				cbind(1, treatment = private$w, X_cov)
			}

			y_beta = private$sanitize_beta_response(private$y)
			res = fast_beta_regression(Xmm = Xmm, y = y_beta)
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			as.numeric(res$b[2])
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		sanitize_beta_response = function(y){
			if (should_run_asserts()) {
				assertNumeric(y, any.missing = FALSE)
			}
			n = length(y)
			if (n == 0L) return(y)
			y = pmin(1, pmax(0, as.numeric(y)))
			(y * (n - 1) + 0.5) / n
		},

		generate_mod = function(estimate_only = FALSE){
			Xmm = cbind(1, private$w)
			full_names = c("(Intercept)", "treatment")
			colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
			y_beta = private$sanitize_beta_response(private$y)
			
			if (estimate_only) {
				res = fast_beta_regression(Xmm = Xmm, y = y_beta)
				# Ensure names are set for shared()
				names(res$b) = colnames(Xmm)
				return(list(
					b = res$b,
					ssq_b_2 = NA_real_
				))
			} else {
				res = fast_beta_regression_with_var(Xmm = Xmm, y = y_beta)
				# Ensure names are set for shared()
				names(res$b) = colnames(Xmm)
				return(list(
					b = res$b,
					ssq_b_2 = res$ssq_b_2
				))
			}
		}
	)
) # End of R6::R6Class
