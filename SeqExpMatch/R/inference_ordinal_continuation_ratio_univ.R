#' Continuation-Ratio Logit Inference based on Maximum Likelihood
#'
#' @description
#' Continuation-ratio logit model inference for ordinal responses. This model
#' conditions on "having reached category k" before comparing category k vs >k.
#' It is particularly useful when categories represent a progression.
#'
#' @export
SeqDesignInferenceOrdinalContRatioRegr = R6::R6Class("SeqDesignInferenceOrdinalContRatioRegr",
	inherit = SeqDesignInferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param	verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},


		#' @description
		#' Computes the appropriate estimate
		#'
		#' @return	The setting-appropriate (see description) numeric estimate of the treatment effect
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		}
		),

	private = list(
		generate_mod = function(){
			# Continuation-ratio model logit(P(Y = k | Y >= k)) = alpha_k + beta * w
			# We use fast_continuation_ratio_regression_with_var_cpp
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = c("treatment")
			res = fast_continuation_ratio_regression_with_var_cpp(X = Xmm, y = as.numeric(private$y))

			# Return in expected format
			list(
				b = c(NA, res$b[1]), # Match the [2] indexing in shared()
				ssq_b_2 = res$ssq_b_2
			)
		}
	)
) # End of R6::R6Class
