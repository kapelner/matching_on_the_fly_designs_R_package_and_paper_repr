#' Zero/One-Inflated Beta Inference for Proportion Responses
#'
#' Internal class for non-KK zero/one-inflated beta regression models. The
#' response is modeled as a three-component mixture with point masses at 0 and 1
#' plus a beta-distributed interior component on \eqn{(0, 1)}. The reported
#' treatment effect is the treatment coefficient from the beta mean submodel, on
#' the logit scale.
#'
#' @details
#' The beta mean submodel uses treatment alone in the univariate class and
#' treatment plus covariates in the multivariate class. The zero/one inflation
#' submodels use \code{model_formula_zero_one}, which defaults to \code{~ .}
#' so that treatment plus all available covariates enter those auxiliary pieces.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'proportion')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferencePropZeroOneInflatedBetaRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferencePropZeroOneInflatedBetaRegr = R6::R6Class("InferencePropZeroOneInflatedBetaRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a zero-one-inflated beta regression inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL}
		#'   (default), the formula from the design object is used and its pre-computed
		#'   design matrix is reused. If a formula is provided, a new design matrix is
		#'   constructed from the design's imputed covariates.
		#' @param model_formula_zero_one Formula for the zero/one inflation submodels.
		#'   Defaults to \code{~ .}, meaning treatment plus all available covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, model_formula_zero_one = ~ ., verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
				assertFormula(model_formula_zero_one, null.ok = FALSE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$model_formula_zero_one = model_formula_zero_one
		}
	),
	private = list(
		best_X_colnames = NULL,
		best_X_zero_one_colnames = NULL,
		model_formula_zero_one = NULL,
		build_component_matrix = function(model_formula, selected_colnames = NULL, treatment_name = "treatment"){
			if (is.null(selected_colnames)) {
				if (identical(model_formula, ~ .)) {
					X_cov = private$get_X()
				} else {
					X_imp = private$des_obj$get_X_imp()
					X_cov = if (is.null(X_imp)) matrix(NA_real_, nrow = private$n, ncol = 0) else create_model_matrix_from_features(model_formula, X_imp)
				}
			} else {
				if (identical(model_formula, ~ .)) {
					X_cov_all = private$get_X()
				} else {
					X_imp = private$des_obj$get_X_imp()
					X_cov_all = if (is.null(X_imp)) matrix(NA_real_, nrow = private$n, ncol = 0) else create_model_matrix_from_features(model_formula, X_imp)
				}
				X_cov = if (is.null(X_cov_all) || !length(selected_colnames)) matrix(NA_real_, nrow = private$n, ncol = 0) else as.matrix(X_cov_all[, intersect(selected_colnames, colnames(X_cov_all)), drop = FALSE])
			}
			if (is.null(X_cov) || ncol(as.matrix(X_cov)) == 0L) {
				X_fit = cbind(`(Intercept)` = 1, treatment = private$w)
				colnames(X_fit)[2] = treatment_name
				return(X_fit)
			}
			X_fit = cbind(`(Intercept)` = 1, treatment = private$w, as.matrix(X_cov))
			colnames(X_fit)[2] = treatment_name
			X_fit
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X = private$build_component_matrix(private$model_formula, private$best_X_colnames)
			X_zero_one = private$build_component_matrix(private$model_formula_zero_one, private$best_X_zero_one_colnames)
			res = tryCatch(
				fast_zero_one_inflated_beta_cpp(X, X_zero_one, private$y, init = rep(0, ncol(X) + 1L + 2L * ncol(X_zero_one))),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 2 || !is.finite(res$b[2])) return(NA_real_)
			as.numeric(res$b[2])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		supports_likelihood_tests = function(){
			TRUE
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			X_zero_one = ctx$X_zero_one
			y = as.numeric(private$y)
			j_treat = as.integer(ctx$j_treat)
			full_params = ctx$full_params
			list(
				X = X_fit, X_zero_one = X_zero_one, y = y, j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta){
					init = if (!is.null(full_params)) as.numeric(full_params) else rep(0, ncol(X_fit) + 1L + 2L * ncol(X_zero_one))
					res = tryCatch(
						fast_zero_one_inflated_beta_cpp(
							X_fit,
							X_zero_one,
							y,
							init = init,
							fixed_idx = j_treat,
							fixed_values = delta
						),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					res
				},
				score = function(fit){
					get_zero_one_inflated_beta_score_cpp(X_fit, X_zero_one, y, as.numeric(fit$params))
				},
				observed_information = function(fit){
					-get_zero_one_inflated_beta_hessian_cpp(X_fit, X_zero_one, y, as.numeric(fit$params))
				},
				fisher_information = function(fit){
					-get_zero_one_inflated_beta_hessian_cpp(X_fit, X_zero_one, y, as.numeric(fit$params))
				},
				information = function(fit){
					-get_zero_one_inflated_beta_hessian_cpp(X_fit, X_zero_one, y, as.numeric(fit$params))
				},
				neg_loglik = function(fit){ as.numeric(fit$neg_loglik) }
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_component_matrix(private$model_formula)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 2L,
				fit_fun = function(X_fit, keep){
					j_treat = which(keep == 2L)
					if (is.null(private$best_X_zero_one_colnames)) {
						X_zero_one_full = private$build_component_matrix(private$model_formula_zero_one)
						red_zo = private$reduce_design_matrix_preserving_treatment(X_zero_one_full)
						X_zero_one = red_zo$X
						if (is.null(X_zero_one)) return(NULL)
						colnames(X_zero_one) = colnames(X_zero_one_full)[red_zo$keep]
					} else {
						X_zero_one = private$build_component_matrix(private$model_formula_zero_one, private$best_X_zero_one_colnames)
					}
					init = rep(0, ncol(X_fit) + 1L + 2L * ncol(X_zero_one))
					res = tryCatch(
						fast_zero_one_inflated_beta_cpp(X_fit, X_zero_one, private$y, init = init),
						error = function(e) NULL
					)
					if (is.null(res)) return(NULL)
					ssq_b_j = if (!is.null(res$vcov) && nrow(res$vcov) >= j_treat) {
						as.numeric(res$vcov[j_treat, j_treat])
					} else {
						NA_real_
					}
					list(
						b = as.numeric(res$b),
						ssq_b_j = ssq_b_j,
						j_treat = j_treat,
						params = as.numeric(res$params),
						neg_loglik = as.numeric(res$neg_loglik),
						X_zero_one = X_zero_one
					)
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = mod$j_treat
					if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$best_X_colnames = setdiff(colnames(attempt$X), c("(Intercept)", "treatment"))
				private$best_X_zero_one_colnames = setdiff(colnames(attempt$fit$X_zero_one), c("(Intercept)", "treatment"))
				private$cached_values$likelihood_test_context = list(
					X = attempt$X,
					X_zero_one = attempt$fit$X_zero_one,
					j_treat = attempt$fit$j_treat,
					full_params = attempt$fit$params
				)
			} else {
				private$cached_values$likelihood_test_context = NULL
			}
			attempt$fit
		}
	)
)
