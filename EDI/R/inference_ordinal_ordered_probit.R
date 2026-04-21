#' Ordered Probit Regression Inference for Ordinal Responses
#'
#' Fits an ordered probit regression for ordinal responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceOrdinalOrderedProbitRegr = R6::R6Class("InferenceOrdinalOrderedProbitRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(
		#' @description
		#' Initialize an ordered-probit inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
				assertFormula(model_formula, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		}
	),

	private = list(
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()

			if (length(Xmm_cols) == 0L){
				# Univariate case
				Xmm = as.matrix(private$w)
				colnames(Xmm) = "treatment"
			} else {
				# Multivariate case
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				Xmm = cbind(treatment = private$w, X_cov)
			}

			res = fast_ordinal_probit_regression_cpp(X = Xmm, y = as.numeric(private$y))
			if (is.null(res) || length(res$b) < 1L || !is.finite(res$b[length(res$b)])){
				return(NA_real_)
			}
			as.numeric(res$b[length(res$b)])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		generate_mod = function(estimate_only = FALSE){
			X_full = private$build_design_matrix()
			
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L, # treatment only (no intercept in these ordinal models)
				fit_fun = function(X_fit){
					if (estimate_only) {
						res = fast_ordinal_probit_regression_cpp(X_fit, private$y)
						list(b = res$b, ssq_b_j = NA_real_)
					} else {
						res = fast_ordinal_probit_regression_with_var_cpp(X_fit, private$y)
						list(b = res$b, ssq_b_j = res$ssq_b_j)
					}
				},
				fit_ok = function(mod, X_fit, keep){
					j_treat = length(mod$b) # treatment is last
					if (is.null(mod) || j_treat < 1L || !is.finite(mod$b[j_treat])) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq_b_j) && mod$ssq_b_j > 0
				}
			)

			if (!is.null(attempt$fit)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X), "treatment")
				# Format for shared(): b[2] is treatment
				list(b = c(0, attempt$fit$b[length(attempt$fit$b)]), ssq_b_2 = attempt$fit$ssq_b_j)
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
