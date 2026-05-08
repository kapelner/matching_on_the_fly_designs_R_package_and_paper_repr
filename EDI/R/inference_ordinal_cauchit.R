#' Cumulative Cauchit Inference for Ordinal Responses
#'
#' Cumulative Cauchit model inference for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalCauchitRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceOrdinalCauchitRegr = R6::R6Class("InferenceOrdinalCauchitRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
						#' @description
		#' Initialize a cumulative-cauchit inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_default Whether to use smart optimizer start values by default.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_default = smart_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),

	private = list(
		best_X_colnames = NULL,

		cauchit_polr_fallback = function(){
			if (!check_package_installed("MASS")) return(NULL)
			y_fac = factor(private$y, levels = sort(unique(private$y)))
			if (length(levels(y_fac)) < 2) return(NULL)
			dat = data.frame(y = y_fac, w = private$w)
			mod = tryCatch(
				MASS::polr(y ~ w, data = dat, method = "cauchit", Hess = TRUE),
				error = function(e) NULL
			)
			if (is.null(mod) || !"w" %in% names(stats::coef(mod))) return(NULL)
			coef_w = as.numeric(stats::coef(mod)["w"])
			var_w = tryCatch(vcov(mod)["w", "w"], error = function(e) NA_real_)
			ssq = if (is.finite(var_w) && var_w > 0) var_w else NA_real_
			list(b = c(NA, coef_w), ssq_b_2 = ssq)
		},

		supports_likelihood_tests = function(){
			TRUE
		},

		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			full_fit = list(params = ctx$full_params, neg_loglik = ctx$full_neg_loglik)
			list(
				X = X_fit, y = y, j = j_treat,
				full_fit = full_fit,
				fit_null = function(delta, start = NULL){
					res = tryCatch(
						fast_ordinal_cauchit_regression_cpp(
							X_fit, y,
							start_params = start %||% private$get_fit_warm_start_for_length("params", length(ctx$full_params)),
							fixed_idx = j_treat, fixed_values = delta,
							smart_start = private$smart_default
						),
						error = function(e) NULL
					)
					if (is.null(res) || length(res) == 0) return(NULL)
					list(params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
				},
				extract_start = function(fit){
					as.numeric(fit$params)
				},
				score = function(fit){
					get_ordinal_cauchit_regression_score_cpp(X_fit, y, as.numeric(fit$params))
				},
				observed_information = function(fit){
					-get_ordinal_cauchit_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				fisher_information = function(fit){
					-get_ordinal_cauchit_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				information = function(fit){
					-get_ordinal_cauchit_regression_hessian_cpp(X_fit, y, as.numeric(fit$params))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},

		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit, keep){
					start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + length(sort(unique(private$y))) - 1L)
					if (estimate_only) {
						res = fast_ordinal_cauchit_regression_cpp(
							X = X_fit, y = as.numeric(private$y),
							start_params = start_params,
							smart_start = private$smart_default
						)
						if (is.null(res) || length(res) == 0) return(NULL)
						list(b = c(NA_real_, as.numeric(res$b)), ssq_b_2 = NA_real_,
						     params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
					} else {
						res = fast_ordinal_cauchit_regression_with_var_cpp(
							X = X_fit, y = as.numeric(private$y),
							start_params = start_params,
							smart_start = private$smart_default
						)
						if (is.null(res) || length(res$b) == 0 || is.na(res$b[1])) return(NULL)
						list(b = c(NA_real_, as.numeric(res$b)), ssq_b_2 = as.numeric(res$ssq_b_2),
						     params = as.numeric(res$params), neg_loglik = as.numeric(res$neg_loglik))
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$params, "params")
				private$best_X_colnames = setdiff(colnames(attempt$X), "treatment")
				n_alpha = length(attempt$fit$params) - ncol(attempt$X)
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					j_treat = as.integer(n_alpha + 1L),
					full_params = attempt$fit$params,
					full_neg_loglik = attempt$fit$neg_loglik
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		},

		build_design_matrix = function(){
			X_cov = private$X
			if (is.null(X_cov) || ncol(X_cov) == 0) {
				X = matrix(private$w, ncol = 1L)
				colnames(X) = "treatment"
			} else {
				X = cbind(treatment = private$w, X_cov)
			}
			X
		}
	)
)
