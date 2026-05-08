#' Adjacent Category Logit Inference for KK Matching-on-the-fly Designs
#'
#' Fits a conditional adjacent category logit model for ordinal responses
#' under a KK matching-on-the-fly design. This model is implemented via data
#' expansion and conditional logistic regression.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKCondAdjCatLogitRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKCondAdjCatLogitRegr = R6::R6Class("InferenceOrdinalKKCondAdjCatLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK \code{DesignSeqOneByOne} object.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose Flag for progress messages.
		#' @param smart_default Whether to use smart optimizer start values by default.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_default = TRUE){
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_default = smart_default)
		},

		#' @description
		#' Returns the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Returns the asymptotic confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared()
			ordinal_cond_clogit_assert_finite_se(private, class(self)[1])
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Returns the asymptotic p-value.
		#' @param delta Null hypothesis treatment effect.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared()
			ordinal_cond_clogit_assert_finite_se(private, class(self)[1])
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		shared = function(estimate_only = FALSE){
			ordinal_cond_clogit_shared_multi(private, expand_adjacent_category_data_cpp, function(y_i, n_alpha) {
				trials = integer(0)
				if (y_i <= n_alpha) trials = c(trials, y_i)
				if (y_i > 1) trials = c(trials, y_i - 1L)
				sort(unique(trials))
			})
		}
	)
)
