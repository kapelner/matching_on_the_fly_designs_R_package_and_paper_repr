#' KK Hurdle Poisson IVWC Inference for Count Responses
#'
#' Internal base class for KK hurdle-Poisson inverse-variance weighted combined
#' inference. The matched-pair component is fit with a hurdle-Poisson mixed model
#' using pair random intercepts, and the reservoir component is fit with an
#' ordinary Poisson log-link regression. The reported treatment effect is on the
#' log-rate scale.
#'
#' @keywords internal
#' @noRd
InferenceAbstractKKHurdlePoissonIVWC = R6::R6Class("InferenceAbstractKKHurdlePoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB} for
		#'   the matched-pair component.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(use_rcpp)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or DesignFixedBinaryMatch.")
				}
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (inherits(des_obj, "DesignFixedBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_matching_structure_computed()
			}
			private$m = des_obj$.__enclos_env__$private$m
			private$use_rcpp = use_rcpp
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !private$use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
		},
		#' @description Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute asymp confidence interval
		#' @param alpha The significance level (default 0.05).
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
		#' @description Compute asymp two sided pval for treatment effect
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				if (should_run_asserts()) {
					stop("Testing non-zero delta is not yet implemented for this class.")
				}
				NA_real_
			}
		},
		#' @description Computes the treatment effect estimate for a bootstrap sample.
		#' @param subject_or_block_weights Row weights for the bootstrap sample.
		#' @param estimate_only If TRUE, skip variance calculations.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			private$shared_combined_bootstrap(subject_or_block_weights, estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		use_rcpp = TRUE,
		max_abs_reasonable_coef = 1e4,
		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the fixed-effect coefficient for "w" directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			X = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(NA_integer_, nrow(X))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)
			beta_m = NA_real_
			ssq_m = NA_real_
			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(X, matched_idx, m_vec, se = FALSE)
				beta_m = res_m$beta_hat
				ssq_m = res_m$se^2
			}
			m_ok = !is.na(beta_m) && is.finite(beta_m) && !is.na(ssq_m) && is.finite(ssq_m) && ssq_m > 0
			beta_r = NA_real_
			ssq_r = NA_real_
			if (length(reservoir_idx) > 1L && length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(X, reservoir_idx, estimate_only = estimate_only)
				beta_r = res_r$beta_hat
				ssq_r = res_r$ssq_hat
			}
			r_ok = !is.na(beta_r) && is.finite(beta_r) && !is.na(ssq_r) && is.finite(ssq_r) && ssq_r > 0
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		},
		compute_basic_match_data = function(){
			private$cached_values$KKstats = .compute_kk_basic_match_data_cached(
				private_env = private,
				des_priv     = private$des_obj_priv_int,
				X = private$get_X(),
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		build_model_matrix = function(){
			if (ncol(as.matrix(private$X)) > 0){
				X = private$create_design_matrix()
				full_names = c("(Intercept)", "w", if (ncol(X) > 2L) paste0("x", seq_len(ncol(X) - 2L)) else NULL)
				colnames(X) = full_names[seq_len(ncol(X))]
			} else {
				X = cbind(1, private$w)
				colnames(X) = c("(Intercept)", "w")
			}
			X
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			X = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(NA_integer_, nrow(X))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)
			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(X, matched_idx, m_vec, se = !estimate_only)
				private$cached_values$beta_T_matched = res_m$beta_hat
				if (!estimate_only) private$cached_values$ssq_beta_T_matched = res_m$se^2 else private$cached_values$ssq_beta_T_matched = 1.0
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
				!is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0
			if (length(reservoir_idx) > 1L &&
				length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(X, reservoir_idx, estimate_only = estimate_only)
				private$cached_values$beta_T_reservoir = res_r$beta_hat
				if (!estimate_only) private$cached_values$ssq_beta_T_reservoir = res_r$ssq_hat else private$cached_values$ssq_beta_T_reservoir = 1.0
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) &&
				!is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cache_nonestimable_estimate("kk_hurdle_poisson_ivwc_no_usable_component")
			}
			invisible(NULL)
		},
		build_glmm_formula = function(dat){
			fixed_terms = setdiff(colnames(dat), c("y", "pair_group"))
			rhs = paste(c(fixed_terms, "(1 | pair_group)"), collapse = " + ")
			stats::as.formula(paste("y ~", rhs))
		},
		fit_hurdle_for_matched_pairs = function(X, matched_idx, m_vec, se = TRUE){
			X_matched = X[matched_idx, , drop = FALSE]
			if (is.null(dim(X_matched)) || ncol(X_matched) < 2L) {
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			reduced = private$reduce_design_matrix_preserving_treatment(X_matched)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			if (private$use_rcpp) {
				res = private$fit_hurdle_for_matched_pairs_rcpp(
					X_fit = X_fit,
					y_fit = private$y[matched_idx],
					group_id = m_vec[matched_idx],
					j_treat = reduced$j_treat,
					se = se
				)
				if (is.finite(res$beta_hat) && (!se || (is.finite(res$se) && res$se > 0))) {
					return(res)
				}
			}
			private$fit_hurdle_for_matched_pairs_glmm_tmb(
				X_fit = X_fit,
				y_fit = private$y[matched_idx],
				group_id = m_vec[matched_idx],
				se = se
			)
		},
		fit_hurdle_for_matched_pairs_rcpp = function(X_fit, y_fit, group_id, j_treat, se = TRUE){
			fit = tryCatch(
				fast_hurdle_poisson_glmm_cpp(
					X = as.matrix(X_fit),
					y = as.numeric(y_fit),
					group_id = as.integer(group_id),
					j_T = as.integer(j_treat - 1L),
					estimate_only = !se,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) return(list(beta_hat = NA_real_, se = NA_real_))
			beta_hat = as.numeric(fit$b[j_treat])
			if (!is.finite(beta_hat) || abs(beta_hat) > private$max_abs_reasonable_coef) {
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			if (!se) return(list(beta_hat = beta_hat, se = 1.0))
			ssq = as.numeric(fit$ssq_b_T)
			if (!is.finite(ssq) || ssq <= 0) return(list(beta_hat = NA_real_, se = NA_real_))
			se_val = sqrt(ssq)
			if (!is.finite(se_val) || se_val <= 0 || se_val > private$max_abs_reasonable_coef) {
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			list(beta_hat = beta_hat, se = se_val)
		},
		fit_hurdle_for_matched_pairs_glmm_tmb = function(X_fit, y_fit, group_id, se = TRUE){
			if (!check_package_installed("glmmTMB")) return(list(beta_hat = NA_real_, se = NA_real_))
			pred_df = as.data.frame(X_fit[, -1, drop = FALSE])
			colnames(pred_df)[1] = "w"
			dat = data.frame(
				y = y_fit,
				pred_df,
				pair_group = factor(group_id)
			)
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)
			formula_cond = private$build_glmm_formula(dat)
			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = stats::as.formula(sub("^y ~ ", "~ ", deparse(formula_cond))),
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat,
						control = glmm_control,
						se = se
					)
				)),
				error = function(e) NULL
			)
			if (is.null(mod) && ncol(dat) > 3L){
				dat = dat[, c("y", "w", "pair_group"), drop = FALSE]
				mod = tryCatch(
					suppressWarnings(suppressMessages(
						glmmTMB::glmmTMB(
							y ~ w + (1 | pair_group),
							ziformula = ~ w + (1 | pair_group),
							family = glmmTMB::truncated_poisson(link = "log"),
							data = dat,
							control = glmm_control,
							se = se
						)
					)),
					error = function(e) NULL
				)
			}
			if (is.null(mod)) return(list(beta_hat = NA_real_, se = NA_real_))
			if (!se){
				beta = glmmTMB::fixef(mod)$cond
				if ("w" %in% names(beta)){
					return(list(beta_hat = as.numeric(beta["w"]), se = 1.0)) # Return dummy SE > 0
				}
				return(list(beta_hat = NA_real_, se = NA_real_))
			}
			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table))) return(list(beta_hat = NA_real_, se = NA_real_))
			beta_hat = as.numeric(coef_table["w", "Estimate"])
			se_val = as.numeric(coef_table["w", "Std. Error"])
			if (!is.finite(beta_hat) || !is.finite(se_val) || se_val <= 0) return(list(beta_hat = NA_real_, se = NA_real_))
			list(beta_hat = beta_hat, se = se_val)
		},
		fit_poisson_for_reservoir = function(X, reservoir_idx, estimate_only = FALSE){
			X_res = X[reservoir_idx, , drop = FALSE]
			if (is.null(dim(X_res)) || ncol(X_res) < 2L) {
				return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			}
			reduced = private$reduce_design_matrix_preserving_treatment(X_res)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			}
			mod = tryCatch({
				if (estimate_only) {
					fast_poisson_regression_cpp(X_fit, private$y[reservoir_idx])
				} else {
					fast_poisson_regression_with_var_cpp(X_fit, private$y[reservoir_idx], j = reduced$j_treat)
				}
			}, error = function(e) NULL)
			if (is.null(mod) || !isTRUE(mod$converged)) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			beta_hat = as.numeric(mod$b[reduced$j_treat])
			if (estimate_only) {
				if (!is.finite(beta_hat)) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
				return(list(beta_hat = beta_hat, ssq_hat = 1))
			}
			ssq_hat = as.numeric(mod$ssq_b_j)
			if (!is.finite(beta_hat) || !is.finite(ssq_hat) || ssq_hat <= 0) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			list(beta_hat = beta_hat, ssq_hat = ssq_hat)
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)))
)
#' Abstract class for Conditional Poisson Combined-Likelihood Compound Inference
#'
#' Internal base class for KK hurdle-Poisson combined-likelihood models.
#'
#' @keywords internal
#' @noRd
InferenceAbstractKKHurdlePoissonOneLik = R6::R6Class("InferenceAbstractKKHurdlePoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB}.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param optimization_alg Optimization algorithm.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(use_rcpp)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$use_rcpp = use_rcpp
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			if (should_run_asserts() && !private$use_rcpp) {
				if (!check_package_installed("glmmTMB")){
					stop("Package 'glmmTMB' is required for ", class(self)[1], " when use_rcpp = FALSE. Please install it.")
				}
			}
			private$init_kk_passthrough(des_obj)
		},
		#' @description Compute treatment estimate
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_hurdle(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute asymp confidence interval
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute asymp two sided pval
		compute_asymp_two_sided_pval = function(delta = 0){
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		use_rcpp = TRUE,
		shared_combined_hurdle = function(estimate_only = FALSE){
			# Placeholder for combined hurdle logic
		}
	)))
)
#' Abstract class for Conditional Poisson / Negative Binomial Compound Inference
#'
#' @keywords internal
InferenceAbstractKKPoissonCPoissonIVWC = R6::R6Class("InferenceAbstractKKPoissonCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = TRUE){
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		shared = function(estimate_only = FALSE){
			# Placeholder for IVWC logic
		}
	)))
)
#' Conditional-Poisson Inference for KK Designs with IVWC
#' @export
InferenceCountKKCPoissonIVWC = R6::R6Class("InferenceCountKKCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonIVWC,
	public = list()
)
#' Conditional-Poisson Inference for KK Designs with Combined Likelihood
#' @export
InferenceCountKKCPoissonOneLik = R6::R6Class("InferenceCountKKCPoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonOneLik,
	public = list()
)
#' KK Hurdle Poisson IVWC Inference for Count Responses
#' @export
InferenceCountKKHurdlePoissonIVWC = R6::R6Class("InferenceCountKKHurdlePoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKHurdlePoissonIVWC,
	public = list()
)
