#' Stereotype Logistic Inference for Ordinal Responses
#'
#' Stereotype logistic inference for ordinal responses. The model allows
#' category-specific treatment effects through latent score weights while keeping
#' a common treatment coefficient.
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
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUniStereotypeLogitRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalUniStereotypeLogitRegr = R6::R6Class("InferenceOrdinalUniStereotypeLogitRegr",
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
			assertResponseType(des_obj$get_response_type(), "ordinal")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		stereotype_design_matrix = function(){
			Xmm = matrix(private$w, ncol = 1)
			full_names = "treatment"
			colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
			Xmm
		},

		generate_mod = function(estimate_only = FALSE){
			X = private$stereotype_design_matrix()
			y = as.numeric(private$y)
			res = fast_stereotype_logit_with_var_cpp(X = X, y = y)
			if (isTRUE(res$converged) && is.finite(res$b[1])){
				if (estimate_only) return(list(b = c(NA, res$b[1]), ssq_b_2 = NA_real_))
				if (is.finite(res$ssq_b_1) && res$ssq_b_1 > 0) return(list(b = c(NA, res$b[1]), ssq_b_2 = res$ssq_b_1))
			}
			fallback = private$stereotype_polr_fallback()
			if (!is.null(fallback)){
				if (estimate_only) return(list(b = fallback$b, ssq_b_2 = NA_real_))
				return(fallback)
			}
			list(b = c(NA, res$b[1]), ssq_b_2 = if (estimate_only) NA_real_ else res$ssq_b_1)
		},

		profile_loglik = function(beta){
			fast_stereotype_profile_loglik_cpp(
				X = private$stereotype_design_matrix(),
				y = as.numeric(private$y),
				beta_fixed = beta
			)
		},

		stereotype_polr_fallback = function(link = "logistic"){
			if (!check_package_installed("MASS")) return(NULL)
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)
			mod = tryCatch(
				MASS::polr(y ~ w, data = dat, method = link, Hess = TRUE),
				error = function(e) NULL
			)
			if (is.null(mod) || !"w" %in% names(stats::coef(mod))) return(NULL)
			coef_w = as.numeric(stats::coef(mod)["w"])
			var_w = tryCatch(vcov(mod)["w", "w"], error = function(e) NA_real_)
			ssq = if (is.finite(var_w) && var_w > 0) var_w else NA_real_
			list(b = c(NA, coef_w), ssq_b_2 = ssq)
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			if (!is.null(private[["custom_randomization_statistic_function"]])) return(NULL)
			w_mat = permutations$w_mat
			if (is.null(w_mat)) return(NULL)
			X_covars = private$get_X()
			compute_stereotype_logit_distr_parallel_cpp(as.numeric(y), X_covars, w_mat, as.numeric(delta), private$n_cpp_threads(ncol(w_mat)))
		}
	)
)
