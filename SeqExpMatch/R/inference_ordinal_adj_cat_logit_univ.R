#' Adjacent-Category Logit Inference for Ordinal Responses
#'
#' @description
#' Adjacent-category logit inference for ordinal responses. The model assumes a
#' common treatment effect across the logits of each category versus the next
#' category.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalUniAdjCatLogitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceOrdinalUniAdjCatLogitRegr = R6::R6Class(
	"SeqDesignInferenceOrdinalUniAdjCatLogitRegr",
	inherit = SeqDesignInferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the adjacent-category treatment effect estimate.
		#'
		#' @return	The estimated treatment log-odds ratio shared across adjacent categories.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		adjacent_category_design_matrix = function(){
			Xmm = matrix(private$w, ncol = 1)
			colnames(Xmm) = "treatment"
			Xmm
		},

		generate_mod = function(){
			res = fast_adjacent_category_logit_with_var_cpp(
				X = private$adjacent_category_design_matrix(),
				y = as.numeric(private$y)
			)
			list(
				b = c(NA, res$b[1]),
				ssq_b_2 = res$ssq_b_1
			)
		}
	)
)
