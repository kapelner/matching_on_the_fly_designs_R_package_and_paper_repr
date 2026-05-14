#' Abstract class for GLMM-based Inference
#'
#' @keywords internal
InferenceAbstractKKGLMM = R6::R6Class("InferenceAbstractKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), private$glmm_response_type())
				assertFormula(model_formula, null.ok = TRUE)
			}
			if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or DesignFixedBinaryMatch.")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (inherits(des_obj, "DesignFixedBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			if (identical(private$glmm_response_type(), "proportion")) {
				private$y = .sanitize_proportion_response(private$y, interior = FALSE)
			}
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !isTRUE(private$skip_glmm_pkg_check)) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
				}
			}
		},
		#' @description Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the asymptotic confidence interval.
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		m = NULL,
		optimization_alg = "lbfgs",
		# Subclasses that provide their own fitter (e.g. Rcpp) set this to TRUE
		# before calling super$initialize() to suppress the glmmTMB package check.
			skip_glmm_pkg_check = FALSE,
		# Abstract: subclasses must return the expected response type string.
		glmm_response_type = function() stop(class(self)[1], " must implement glmm_response_type()"),
		# Abstract: subclasses must return the glm family object for glmmTMB.
		glmm_family = function() stop(class(self)[1], " must implement glmm_family()"),
		# Default (multivariate): all covariates + treatment.
		# Univariate subclasses override this to return data.frame(w = private$w).
		glmm_predictors_df = function(){
			df = as.data.frame(private$create_design_matrix()[, -1, drop = FALSE])
			# create_design_matrix uses "treatment"; glmmTMB path expects "w"
			if ("treatment" %in% colnames(df))
				colnames(df)[colnames(df) == "treatment"] = "w"
			df
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
			
			candidates = list(as.data.frame(attempt$X))
			if (ncol(attempt$X) > 1L){
				if (!("w" %in% unlist(lapply(candidates, colnames), use.names = FALSE))){
					candidates[[length(candidates) + 1L]] = data.frame(w = predictors_df$w)
				}
			}
			candidates
		},
			max_abs_reasonable_coef = 1e4,
			get_standard_error = function(){
				private$shared(estimate_only = FALSE)
				se = private$compute_standard_error_from_information_matrix()
				if (is.finite(se)) return(se)
				private$cached_values$s_beta_hat_T
			},
			get_degrees_of_freedom = function(){
				private$shared(estimate_only = FALSE)
				private$cached_values$df
			},
			supports_likelihood_tests = function(){
				isTRUE(private$use_rcpp)
			},
			get_likelihood_test_spec = function(){
				private$shared(estimate_only = FALSE)
				ctx = private$cached_values$likelihood_test_context
				if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
				X_fit = ctx$X
				y = as.numeric(ctx$y)
				group_id = as.integer(ctx$group_id)
				j_treat = as.integer(ctx$j_treat)
				j_T = as.integer(ctx$j_T)
				n_gh = as.integer(ctx$n_gh %||% 20L)
				start_len = length(ctx$start)
				list(
					X = X_fit,
					y = y,
					group_id = group_id,
					j = j_treat,
					full_fit = private$cached_mod,
					fit_null = function(delta, start = NULL){
						start_par = start %||% private$get_fit_warm_start_for_length("params", start_len) %||% ctx$start
						warm_fisher = private$get_fit_warm_start_fisher(start_len)
						fast_logistic_glmm_cpp(
							X = X_fit,
							y = y,
							group_id = group_id,
							j_T = j_T,
							start_par = start_par,
							warm_start_fisher_info = warm_fisher,
							estimate_only = FALSE,
							n_gh = n_gh,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					},
					extract_start = function(fit){
						as.numeric(fit$params %||% c(fit$b, fit$log_sigma))
					},
					score = function(fit){
						as.numeric(fit$score %||% get_logistic_glmm_score_cpp(X_fit, y, group_id, as.numeric(fit$params %||% c(fit$b, fit$log_sigma)), n_gh = n_gh))
					},
					observed_information = function(fit){
						as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% get_logistic_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(fit$params %||% c(fit$b, fit$log_sigma)), n_gh = n_gh))
					},
					fisher_information = function(fit){
						as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% get_logistic_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(fit$params %||% c(fit$b, fit$log_sigma)), n_gh = n_gh))
					},
					information = function(fit){
						as.matrix(fit$fisher_information %||% fit$information %||% fit$observed_information %||% get_logistic_glmm_hessian_cpp(X_fit, y, group_id, as.numeric(fit$params %||% c(fit$b, fit$log_sigma)), n_gh = n_gh))
					},
					neg_loglik = function(fit){
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_logistic_glmm_neg_loglik_cpp(X_fit, y, group_id, as.numeric(fit$params %||% c(fit$b, fit$log_sigma)), n_gh = n_gh))
					}
				)
			},
			shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			private$clear_nonestimable_state()
			mod = private$fit_glmm(se = !estimate_only)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_glmm_fit_failed")
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
