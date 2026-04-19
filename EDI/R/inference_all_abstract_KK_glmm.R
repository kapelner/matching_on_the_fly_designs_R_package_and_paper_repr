#' Abstract class for GLMM-based Inference
#'
#' @keywords internal
InferenceAbstractKKGLMM = R6::R6Class("InferenceAbstractKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object (must be a KK design).
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as fixed-effect predictors. If \code{FALSE}, only the treatment
		#'   indicator is used. If \code{NULL} (default), it is set to \code{TRUE} if the
		#'   design contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), private$glmm_response_type())
				assertFlag(include_covariates, null.ok = TRUE)
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or FixedDesignBinaryMatch.")
				}
			}
			super$initialize(des_obj, verbose)
			if (is(des_obj, "FixedDesignBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_bms_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			if (identical(private$glmm_response_type(), "proportion")) {
				private$y = .sanitize_proportion_response(private$y, interior = FALSE)
			}
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts()) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
				}
			}
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
				private$shared(estimate_only = TRUE)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		m = NULL,
		include_covariates = NULL,

		# Abstract: subclasses must return the expected response type string.
		glmm_response_type = function() stop(class(self)[1], " must implement glmm_response_type()"),

		# Abstract: subclasses must return the glm family object for glmmTMB.
		glmm_family = function() stop(class(self)[1], " must implement glmm_family()"),

		# Default (multivariate): all covariates + treatment.
		# Univariate subclasses override this to return data.frame(w = private$w).
		glmm_predictors_df = function(){
			if (private$include_covariates) {
				as.data.frame(private$create_design_matrix()[, -1, drop = FALSE])
			} else {
				data.frame(w = private$w)
			}
		},

		glmm_predictors_df_candidates = function(){
			predictors_df = private$glmm_predictors_df()
			if (!private$harden || is.null(predictors_df) || ncol(predictors_df) <= 1L){
				return(list(predictors_df))
			}

			X_full = as.matrix(predictors_df)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = match("w", colnames(X_full)),
				fit_fun = function(X_fit){
					private$fit_glmm_on_data(as.data.frame(X_fit), se = TRUE)
				},
				fit_ok = function(mod, X_fit, keep){
					private$.is_usable_glmm_fit(mod, se = TRUE)
				}
			)
			
			if (is.null(attempt$fit)) return(list(data.frame(w = predictors_df$w)))
			
			candidates = list(as.data.frame(attempt$X_fit))
			if (ncol(attempt$X_fit) > 1L){
				if (!("w" %in% unlist(lapply(candidates, colnames), use.names = FALSE))){
					candidates[[length(candidates) + 1L]] = data.frame(w = predictors_df$w)
				}
			}
			candidates
		},

		max_abs_reasonable_coef = 1e4,

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

			mod = private$fit_glmm(se = !estimate_only)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_glmm_fit_failed")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			# glmmTMB fixed effects for the conditional model
			beta = glmmTMB::fixef(mod)$cond
			if ("w" %in% names(beta)){
				private$cached_values$beta_hat_T = as.numeric(beta["w"])
			} else {
				private$cached_values$beta_hat_T = NA_real_
			}

			if (estimate_only) return(invisible(NULL))

			coef_table = summary(mod)$coefficients$cond
			if ("w" %in% rownames(coef_table) && "Std. Error" %in% colnames(coef_table)){
				private$cached_values$s_beta_hat_T = as.numeric(coef_table["w", "Std. Error"])
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}

			private$cached_values$is_z = TRUE
			private$cached_values$df = Inf
			private$cached_values$summary_table = coef_table
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				return(invisible(NULL))
		},

		fit_glmm_on_data = function(predictors_df, se = TRUE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)
			dat = data.frame(y = private$y, predictors_df, group_id = factor(group_id))
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			glmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))
			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					glmmTMB::glmmTMB(
						glmm_formula,
						family  = private$glmm_family(),
						data    = dat,
						control = glmm_control,
						se      = se
					)
				)))
				mod
			}, error = function(e) {
				message(paste("GLMM FIT ERROR:", e$message))
				NULL
			})
		},

		fit_glmm = function(se = TRUE){
			for (predictors_df in private$glmm_predictors_df_candidates()){
				mod = private$fit_glmm_on_data(predictors_df, se = se)
				if (private$.is_usable_glmm_fit(mod, se)) return(mod)
			}
			NULL
		},

		.is_usable_glmm_fit = function(mod, se){
			if (is.null(mod)) return(FALSE)
			beta = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
			if (is.null(beta) || !("w" %in% names(beta)) || !is.finite(beta["w"])) return(FALSE)
			if (!se) return(TRUE)
			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table)) || !("Std. Error" %in% colnames(coef_table))) return(FALSE)
			se_w = suppressWarnings(as.numeric(coef_table["w", "Std. Error"]))
			is.finite(se_w) && se_w > 0
		}
	)
)
