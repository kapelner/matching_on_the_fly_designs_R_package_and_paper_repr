#' Survival Transformation Regression with Dependent Censoring
#'
#' Fits a lognormal transformation model for survival responses that jointly models
#' the event and censoring times.
#'
#' @export
InferenceSurvivalDepCensTransformRegr = R6::R6Class("InferenceSurvivalDepCensTransformRegr",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMSummaryTable,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose)
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
		}
	),

	private = list(
		include_covariates = NULL,

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},

		build_design_matrix = function(){
			if (private$include_covariates) {
				private$create_design_matrix()
			} else {
				X = matrix(private$w, ncol = 1)
				colnames(X) = "treatment"
				X
			}
		},

		build_design_matrix_candidates = function(){
			list(private$build_design_matrix())
		},

		generate_mod = function(estimate_only = FALSE){
			for (Xmm in private$build_design_matrix_candidates()){
				treatment_col = match("treatment", colnames(Xmm))
				if (!is.finite(treatment_col)) treatment_col = 1L
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = Xmm,
					required_cols = treatment_col,
					fit_fun = function(X_fit){
						.fit_dep_cens_transform_model(
							y = private$y,
							dead = private$dead,
							Xmm = X_fit,
							estimate_only = estimate_only
						)
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod) || is.null(mod$coefficients)) return(FALSE)
						if (!"treatment" %in% names(mod$coefficients)) return(FALSE)
						if (!is.finite(as.numeric(mod$coefficients["treatment"]))) return(FALSE)
						if (estimate_only) return(TRUE)
						if (is.null(mod$vcov) || !"treatment" %in% rownames(mod$vcov)) return(FALSE)
						treatment_var = as.numeric(mod$vcov["treatment", "treatment"])
						is.finite(treatment_var) && treatment_var > 0
					}
				)
				if (!is.null(attempt$fit)) return(attempt$fit)
			}
			warning("Dependent-censoring transformation model failed to converge; returning NA.")
			private$failed_dep_cens_transform_mod()
		},

		failed_dep_cens_transform_mod = function(){
			coef_names = "treatment"
			coef = stats::setNames(NA_real_, coef_names)
			vcov_mat = matrix(NA_real_, nrow = 1L, ncol = 1L,
				dimnames = list(coef_names, coef_names))
			list(
				coefficients = coef,
				vcov = vcov_mat
			)
		}
	)
)
