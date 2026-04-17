#' Adjacent-Category Logit Inference for Ordinal Responses
#'
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
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = nrow(x_dat), response_type = "ordinal",
#'   verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUniAdjCatLogitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalUniAdjCatLogitRegr = R6::R6Class(
	"InferenceOrdinalUniAdjCatLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose                 A flag indicating whether messages should be displayed
		#'   to the user. Default is \code{TRUE}.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj,  verbose = FALSE, harden = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose, harden)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
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

			y_ord = as.integer(factor(private$y, ordered = TRUE))
			K = max(y_ord)
			if (K < 2L) return(NA_real_)

			expanded = expand_adjacent_category_data_cpp(as.integer(y_ord), as.integer(private$w), as.integer(seq_len(private$n)), as.integer(K))
			
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			
			if (length(Xmm_cols) == 0L){
				# Univariate case: treatment only
				mod = clogit_helper(expanded$y, data.frame(), expanded$w, expanded$strata)
			} else {
				# Multivariate case
				X_stack_list = list()
				for (i in 1:private$n) {
					yi = y_ord[i]
					Xi = X_data[i, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
					trials_i = integer()
					for (j in 1 : (K-1L)) {
						if (yi == j || yi == j + 1) trials_i = c(trials_i, j)
					}
					if (length(trials_i) > 0) {
						X_stack_list[[i]] = matrix(rep(Xi, each = length(trials_i)), nrow = length(trials_i))
					}
				}
				X_stack = do.call(rbind, X_stack_list)
				mod = clogit_helper(expanded$y, as.data.frame(X_stack), expanded$w, expanded$strata)
			}

			if (is.null(mod) || !is.finite(mod$b[1])) return(NA_real_)
			as.numeric(mod$b[1])
		},

		adjacent_category_design_matrix = function(){
			Xmm = matrix(private$w, ncol = 1)
			full_names = "treatment"
			colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
			Xmm
		},

		generate_mod = function(estimate_only = FALSE){
			if (estimate_only) {
				res = fast_adjacent_category_logit_cpp(
					X = private$adjacent_category_design_matrix(),
					y = as.numeric(private$y)
				)
				return(list(
					b = c(NA, res$b[1]),
					ssq_b_2 = NA_real_
				))
			}

			res = fast_adjacent_category_logit_with_var_cpp(
				X = private$adjacent_category_design_matrix(),
				y = as.numeric(private$y)
			)
			list(
				b = c(NA, res$b[1]),
				ssq_b_2 = res$ssq_b_1
			)
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			if (is.null(w_mat)) return(NULL)
			X_covars = private$get_X()
			compute_adj_cat_logit_distr_parallel_cpp(as.numeric(y), X_covars, w_mat, as.numeric(delta), private$n_cpp_threads(ncol(w_mat)))
		}
	)
)
