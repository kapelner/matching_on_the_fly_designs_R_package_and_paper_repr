#' Stereotype Logit Regression Inference for Ordinal Responses
#'
#' Fits a stereotype logit regression for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @name InferenceOrdinalStereotypeLogitRegr
#' @description Internal base class for stereotype logit regression.
#' @keywords internal
InferenceOrdinalStereotypeLogitRegr = R6::R6Class("InferenceOrdinalStereotypeLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a stereotype logit inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		#' @param smart_cold_start_default Whether to use smart cold starts.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),
	private = list(
		best_Xmm_colnames = NULL,
		get_complexity_tier = function() "heavy",
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			if (length(X_cols) == 0L){
				X = as.matrix(private$w)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}

			n_params = (length(sort(unique(private$y))) - 1L) + ncol(X)
			if (length(sort(unique(private$y))) >= 3) n_params = n_params + (length(sort(unique(private$y))) - 2L)
			
			ws_args = private$get_backend_warm_start_args(n_params)
			
			res = tryCatch(
				fast_stereotype_logit_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_params = ws_args$start_params,
					warm_start_fisher_info = ws_args$warm_start_fisher_info
				),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[1])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$params, "params", fisher = res$fisher_information)
			as.numeric(res$b[1])
		},
		supports_reusable_bootstrap_worker = function(){
			TRUE
		},
		get_bootstrap_worker_spec = function(){
			private$shared(estimate_only = FALSE)
			list(
				X_full = private$build_design_matrix(),
				best_X_cols = private$best_Xmm_colnames,
				fit_fun = function(X_fit, keep){
					K = length(sort(unique(private$y)))
					n_params = (K - 1L) + ncol(X_fit)
					if (K >= 3) n_params = n_params + (K - 2L)
					ws_args = private$get_backend_warm_start_args(n_params)
					res = fast_stereotype_logit_cpp(
						X_fit, private$y,
						warm_start_params = ws_args$start_params,
						warm_start_fisher_info = ws_args$warm_start_fisher_info
					)
					list(b = res$b, ssq_b_j = NA_real_, params = res$params, fisher_information = res$fisher_information)
				}
			)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					K = length(sort(unique(private$y)))
					n_params = (K - 1L) + ncol(X_fit)
					if (K >= 3) n_params = n_params + (K - 2L)
					ws_args = private$get_backend_warm_start_args(n_params)
					if (estimate_only) {
						res = fast_stereotype_logit_cpp(
							X_fit, private$y,
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(b = res$b, ssq_b_j = NA_real_, params = res$params, fisher_information = res$fisher_information)
					} else {
						res = fast_stereotype_logit_with_var_cpp(
							X_fit, private$y,
							warm_start_params = ws_args$start_params,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(b = res$b, ssq_b_j = res$ssq_b_1, params = res$params, fisher_information = res$fisher_information)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 1L || !is.finite(mod$b[1])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$params, "params", fisher = attempt$fit$fisher_information)
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), "treatment")
				list(b = c(0, attempt$fit$b[1]), ssq_b_2 = attempt$fit$ssq_b_j)
			} else {
				NULL
			}
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

#' Univariate Stereotype Logit Regression Inference for Ordinal Responses
#' @export
InferenceOrdinalUniStereotypeLogitRegr = R6::R6Class("InferenceOrdinalUniStereotypeLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalStereotypeLogitRegr,
	private = list(
		ppo_covariate_matrix = function(){
			matrix(0, nrow = private$n, ncol = 0)
		}
	)
)

#' Multivariate Stereotype Logit Regression Inference for Ordinal Responses
#' @export
InferenceOrdinalMultiStereotypeLogitRegr = R6::R6Class("InferenceOrdinalMultiStereotypeLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalStereotypeLogitRegr,
	private = list(
		ppo_covariate_matrix = function(){
			private$get_X()
		}
	)
)

#' Continuation Ratio Regression Inference for Ordinal Responses
#'
#' Fits a continuation ratio regression for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @name InferenceOrdinalContRatioRegr
#' @description Internal base class for continuation ratio regression.
#' @keywords internal
InferenceOrdinalContRatioRegr = R6::R6Class("InferenceOrdinalContRatioRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize a continuation ratio inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param verbose Whether to print progress messages.
		#' @param harden Whether to apply robustness measures.
		#' @param smart_cold_start_default Whether to use smart cold starts.
		initialize = function(des_obj, verbose = FALSE, harden = TRUE, model_formula = NULL, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),
	private = list(
		best_Xmm_colnames = NULL,
		get_complexity_tier = function() "heavy",
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}
			X_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			if (length(X_cols) == 0L){
				X = as.matrix(private$w)
				colnames(X) = "treatment"
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				X = cbind(treatment = private$w, X_cov)
			}

			n_params = (length(sort(unique(private$y))) - 1L) + ncol(X)
			ws_args = private$get_backend_warm_start_args(n_params)
			res = tryCatch(
				fast_continuation_ratio_regression_cpp(
					X = X, y = as.numeric(private$y),
					warm_start_beta = ws_args$warm_start_beta,
					warm_start_fisher_info = ws_args$warm_start_fisher_info
				),
				error = function(e) NULL
			)
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[length(res$b)])){
				return(NA_real_)
			}
			private$set_fit_warm_start(res$b, "beta", fisher = res$fisher_information)
			as.numeric(res$b[length(res$b)])
		},
		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					n_params = ncol(X_fit) + length(sort(unique(private$y))) - 1L
					ws_args = private$get_backend_warm_start_args(n_params)
					if (estimate_only) {
						res = fast_continuation_ratio_regression_cpp(
							X_fit, private$y,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(b = res$b, ssq_b_j = NA_real_, fisher_information = res$fisher_information)
					} else {
						res = fast_continuation_ratio_regression_with_var_cpp(
							X_fit, private$y,
							warm_start_beta = ws_args$warm_start_beta,
							warm_start_fisher_info = ws_args$warm_start_fisher_info
						)
						list(b = res$b, ssq_b_j = res$ssq_b_1, fisher_information = res$fisher_information)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || length(mod$b) < 1L || !is.finite(mod$b[1])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)
			if (!is.null(attempt$fit)){
				private$set_fit_warm_start(attempt$fit$b, "beta", fisher = attempt$fit$fisher_information)
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), "treatment")
				list(b = c(0, attempt$fit$b[1]), ssq_b_2 = attempt$fit$ssq_b_j)
			} else {
				NULL
			}
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

#' Univariate Continuation Ratio Regression Inference
#' @export
InferenceOrdinalUniContRatioRegr = R6::R6Class("InferenceOrdinalUniContRatioRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalContRatioRegr
)

#' Multivariate Continuation Ratio Regression Inference
#' @export
InferenceOrdinalMultiContRatioRegr = R6::R6Class("InferenceOrdinalMultiContRatioRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalContRatioRegr
)

#' Univariate Cumulative Probit Regression Inference
#' @export
InferenceOrdinalUniCumulProbitRegr = R6::R6Class("InferenceOrdinalUniCumulProbitRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalOrderedProbitRegr
)

#' Multivariate Cumulative Probit Regression Inference
#' @export
InferenceOrdinalMultiCumulProbitRegr = R6::R6Class("InferenceOrdinalMultiCumulProbitRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalOrderedProbitRegr
)
