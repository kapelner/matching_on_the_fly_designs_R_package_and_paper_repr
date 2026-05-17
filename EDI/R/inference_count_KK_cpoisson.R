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
	public = list(
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
		}
	),
	private = list(
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
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
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
	)
)
#' KK Hurdle Poisson Combined-Likelihood Inference for Count Responses
#'
#' Internal base class for KK hurdle-Poisson combined-likelihood models. The
#' count component is a zero-truncated Poisson GLMM fitted over all subjects,
#' with pair-specific random intercepts for matched subjects and singleton
#' groups for reservoir subjects. The reported treatment effect is the
#' treatment coefficient on the log-rate scale.
#'
#' When \code{use_rcpp = TRUE} (default) the likelihood is maximised by an
#' internal Rcpp/L-BFGS routine. Set \code{use_rcpp = FALSE} to fall back to
#' \pkg{glmmTMB}.
#'
#' @keywords internal
#' @noRd
InferenceAbstractKKHurdlePoissonOneLik = R6::R6Class("InferenceAbstractKKHurdlePoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = utils::modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB}.
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
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_hurdle(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Compute asymp confidence interval
		#' @param alpha The significance level (default 0.05).
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(private$fallback_bootstrap_ci(alpha = alpha))
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute asymp two sided pval for treatment effect
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(private$fallback_bootstrap_pval(delta = delta))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
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
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		use_rcpp = TRUE,
		max_abs_reasonable_coef = 1e4,
		warn_bootstrap_fallback_once = function(){
				if (!isTRUE(private$cached_values$warned_bootstrap_se_unavailable)) {
					private$cached_values$warned_bootstrap_se_unavailable = TRUE
					warning("KK hurdle-Poisson combined-likelihood: falling back to bootstrap because standard error is unavailable.")
				}
			},
			fallback_bootstrap_ci = function(alpha){
				private$warn_bootstrap_fallback_once()
				self$compute_bootstrap_confidence_interval(alpha = alpha)
			},
			fallback_bootstrap_pval = function(delta){
				private$warn_bootstrap_fallback_once()
				self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE)
			},
			get_standard_error = function(){
				private$shared_combined_hurdle(estimate_only = FALSE)
				se = private$compute_standard_error_from_information_matrix()
				if (is.finite(se)) return(se)
				private$cached_values$s_beta_hat_T
			},
			get_degrees_of_freedom = function(){
				private$shared_combined_hurdle(estimate_only = FALSE)
				private$cached_values$df
			},
			supports_likelihood_tests = function(){
				isTRUE(private$use_rcpp)
			},
			get_likelihood_test_spec = function(){
				private$shared_combined_hurdle(estimate_only = FALSE)
				ctx = private$cached_values$likelihood_test_context
				if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
				X_fit = ctx$X
				y = as.numeric(ctx$y)
				group_id = as.integer(ctx$group_id)
				j_treat = as.integer(ctx$j_treat)
				j_T = as.integer(ctx$j_T)
				n_gh = as.integer(ctx$n_gh %||% 7L)
				list(
					X = X_fit,
					y = y,
					group_id = group_id,
					j = j_treat,
					full_fit = private$cached_mod,
					fit_null = function(delta, start = NULL){
						warm_start_params = start %||% private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L)
						warm_fisher = private$get_fit_warm_start_fisher(ncol(X_fit) + 1L)
						fast_hurdle_poisson_glmm_cpp(
							X = X_fit,
							y = y,
							group_id = group_id,
							j_T = j_T,
							warm_start_params = warm_start_params,
							warm_start_fisher_info = warm_fisher,
							smart_cold_start = private$smart_cold_start_default,
							estimate_only = FALSE,
							n_gh = n_gh,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					},
					extract_start = function(fit){
						as.numeric(fit$params)
					},
					score = function(fit){
						as.numeric(fit$score %||% get_hurdle_poisson_glmm_score_cpp(X_fit, y, group_id, as.numeric(fit$params), n_gh = n_gh))
					},
					observed_information = function(fit){
						as.matrix(fit$information)
					},
					information = function(fit){
						as.matrix(fit$information)
					},
					neg_loglik = function(fit){
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_hurdle_poisson_glmm_neg_loglik_cpp(X_fit, y, group_id, as.numeric(fit$params), n_gh = n_gh))
					}
				)
			},
			compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
				private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		# ── Group ID construction (matched pairs + singletons for reservoir) ──────
		build_group_ids = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			as.integer(group_id)
		},
		# ── Rcpp path ────────────────────────────────────────────────────────────
			shared_rcpp = function(estimate_only = FALSE){
				if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
				if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
				private$clear_nonestimable_state()
				private$cached_values$likelihood_test_context = NULL
				group_id = private$build_group_ids()
			if (ncol(as.matrix(private$X)) > 0){
				X_fit = private$create_design_matrix()
				if ("treatment" %in% colnames(X_fit))
					colnames(X_fit)[colnames(X_fit) == "treatment"] = "w"
			} else {
				X_fit = cbind(`(Intercept)` = 1, w = private$w)
			}
			X_fit = as.matrix(X_fit)
			j_T = 1L  # 0-based index of treatment column (w is always column 2 = index 1)
			fit = tryCatch(
				fast_hurdle_poisson_glmm_cpp(
					X             = X_fit,
					y             = as.numeric(private$y),
					group_id      = as.integer(group_id),
					j_T           = j_T,
					warm_start_params = private$get_fit_warm_start_for_length("params", ncol(X_fit) + 1L),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(X_fit) + 1L),
					estimate_only = estimate_only,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				return(private$shared_glmm_tmb(estimate_only = estimate_only))
			}
			beta_hat_T = as.numeric(fit$b[j_T + 1L])
				if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
					return(private$shared_glmm_tmb(estimate_only = estimate_only))
				}
				private$cached_mod = fit
				private$set_fit_warm_start(as.numeric(fit$params), "params", fisher = fit$fisher_information)
				private$cached_values$likelihood_test_context = list(
					X = X_fit,
					y = as.numeric(private$y),
					group_id = as.integer(group_id),
					j_T = j_T,
					j_treat = j_T + 1L,
					n_gh = 7L
				)
				private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		},
		# ── glmmTMB fallback path ─────────────────────────────────────────────────
		shared_glmm_tmb = function(estimate_only = FALSE){
			if (!check_package_installed("glmmTMB")){
				private$cache_nonestimable_estimate("kk_hurdle_poisson_combined_no_glmmTMB")
				return(invisible(NULL))
			}
			model_data = private$build_model_data()
			if (is.null(model_data)){
				private$cache_nonestimable_estimate("kk_hurdle_poisson_combined_no_model_data")
				return(invisible(NULL))
			}
			mod = private$fit_hurdle_model_glmm(model_data)
			if (is.null(mod)){
				private$cache_nonestimable_estimate("kk_hurdle_poisson_combined_fit_unavailable")
				return(invisible(NULL))
			}
			cond_fixef = tryCatch(glmmTMB::fixef(mod)$cond, error = function(e) NULL)
			if (is.null(cond_fixef) || !("w" %in% names(cond_fixef)) || !is.finite(cond_fixef[["w"]])){
				private$cache_nonestimable_estimate("kk_hurdle_poisson_combined_treatment_missing")
				return(invisible(NULL))
			}
			private$cached_values$beta_hat_T = as.numeric(cond_fixef[["w"]])
			if (!estimate_only) {
				coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
				se = if (!is.null(coef_table) && ("w" %in% rownames(coef_table))) as.numeric(coef_table["w", "Std. Error"]) else NA_real_
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			}
		},
		shared_combined_hurdle = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				private$shared_glmm_tmb(estimate_only)
			}
		},
		# ── glmmTMB helpers ───────────────────────────────────────────────────────
		build_model_data = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			X_cov = if (!is.null(private$X) && ncol(as.matrix(private$X)) > 0) as.data.frame(private$X) else NULL
			dat = data.frame(
				y          = private$y,
				w          = private$w,
				pair_group = factor(ifelse(m_vec > 0L, m_vec, private$n + seq_along(m_vec)))
			)
			if (!is.null(X_cov)) dat = cbind(dat, X_cov)
			has_pairs = any(m_vec > 0L)
			list(dat = dat, has_pairs = has_pairs)
		},
		build_cond_formula = function(dat, has_pairs = TRUE){
			fixed_terms = setdiff(colnames(dat), c("y", "pair_group"))
			if (has_pairs){
				rhs = paste(c(fixed_terms, "(1 | pair_group)"), collapse = " + ")
			} else {
				rhs = paste(fixed_terms, collapse = " + ")
			}
			stats::as.formula(paste("y ~", rhs))
		},
		build_zi_formula = function(dat, has_pairs = TRUE){
			fixed_terms = setdiff(colnames(dat), c("y", "pair_group"))
			if (has_pairs){
				rhs = paste(c(fixed_terms, "(1 | pair_group)"), collapse = " + ")
			} else {
				rhs = paste(fixed_terms, collapse = " + ")
			}
			stats::as.formula(paste("~", rhs))
		},
		fit_hurdle_model_glmm = function(model_data){
			dat       = model_data$dat
			has_pairs = isTRUE(model_data$has_pairs)
			glmm_control = glmmTMB::glmmTMBControl(parallel = self$num_cores)
			formula_cond = private$build_cond_formula(dat, has_pairs = has_pairs)
			formula_zi   = private$build_zi_formula(dat, has_pairs = has_pairs)
			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = formula_zi,
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
			if (!is.null(mod)) return(mod)
			# Fallback: treatment + pair_group only
			dat_fallback = dat[, intersect(c("y", "w", "pair_group"), colnames(dat)), drop = FALSE]
			formula_cond = private$build_cond_formula(dat_fallback, has_pairs = has_pairs)
			formula_zi   = private$build_zi_formula(dat_fallback, has_pairs = has_pairs)
			tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = formula_zi,
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat_fallback,
						control = glmm_control
					)
				)),
				error = function(e) NULL
			)
		}
	))
)
#' Abstract class for Conditional Poisson / Negative Binomial Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' count responses. For matched pairs, it uses conditional Poisson regression
#' (equivalent to weighted logistic regression) and for the reservoir, it uses
#' standard negative binomial regression. The two estimates are combined via
#' inverse-variance weighted combination (IVWC).
#'
#' @keywords internal
InferenceAbstractKKPoissonCPoissonIVWC = R6::R6Class("InferenceAbstractKKPoissonCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = utils::modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the estimated treatment effect (log-rate ratio).
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
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. Default is 0.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		},
		#' @description Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			inf_obj = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			inf_obj
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Use the same joint-likelihood logic for the point estimate
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_treatment_estimate_during_power_simulation = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = 1e-4){
			private$compute_treatment_estimate_during_power_simulation_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (m > 0){
				private$cpoisson_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) && 
			       (!estimate_only && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0 || estimate_only)
			if (nRT > 0 && nRC > 0){
				private$negbin_for_reservoir(estimate_only = estimate_only)
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) && 
			       (!estimate_only && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0 || estimate_only)
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
				if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
		},
		cpoisson_for_matched_pairs = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			yT = KKstats$yTs_matched
			yC = KKstats$yCs_matched
			# Filter pairs where total count is zero (provide no information for the conditional likelihood)
			y_total = yT + yC
			valid_idx = which(y_total > 0)
			if (length(valid_idx) == 0) return(invisible(NULL))
			y_prop = yT[valid_idx] / y_total[valid_idx]
			weights = y_total[valid_idx]
			# If no covariates, use weighted intercept-only logistic
			X = matrix(1, nrow = length(valid_idx), ncol = 1)
			colnames(X) = "(Intercept)"
			if (ncol(as.matrix(private$X)) > 0 && !is.null(private$best_X_colnames_matched)){
				X = cbind(1, KKstats$X_matched_diffs[valid_idx, intersect(private$best_X_colnames_matched, colnames(KKstats$X_matched_diffs)), drop = FALSE])
			}
			mod = tryCatch({
				res = fast_logistic_regression_weighted_cpp(X = X, y = y_prop, weights = weights)
				ssq_b_1 = if (estimate_only) NA_real_ else {
					vcov_mat = tryCatch(solve(res$XtWX), error = function(e) NULL)
					if (!is.null(vcov_mat)) vcov_mat[1L, 1L] else NA_real_
				}
				list(b = res$b, ssq_b_1 = ssq_b_1, X_fit = X)
			}, error = function(e) NULL)
			if (!is.null(mod) && is.finite(mod$b[1])){
				private$cached_values$beta_T_matched = as.numeric(mod$b[1])
				private$best_X_colnames_matched = setdiff(colnames(mod$X_fit), "(Intercept)")
				if (!estimate_only) private$cached_values$ssq_beta_T_matched = as.numeric(mod$ssq_b_1)
			}
		},
		negbin_for_reservoir = function(estimate_only = FALSE){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L
			if (ncol(as.matrix(private$X)) > 0){
				X_full = cbind(1, w_r, X_r)
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = 2L, # intercept and treatment
					fit_fun = function(X_fit){
						fast_neg_bin_with_var_cpp(X = X_fit, y = as.integer(y_r))
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod) || !is.finite(mod$b[2L])) return(FALSE)
						if (estimate_only) return(TRUE)
						j_col = which(colnames(X_fit) == "w_r")
						if (length(j_col) == 0L) j_col = 2L
						ssq = mod$vcov[j_col, j_col]
						is.finite(ssq) && ssq > 0
					}
				)
				mod = attempt$fit
				if (!is.null(mod)){
					j_treat = which(colnames(attempt$X) == "w_r")
					if (length(j_treat) == 0) j_treat = 2L
				}
			} else {
				X = cbind(1, w_r)
				mod = tryCatch(fast_neg_bin_with_var_cpp(X = X, y = as.integer(y_r)), error = function(e) NULL)
			}
			if (!is.null(mod) && is.finite(mod$b[j_treat])){
				private$cached_values$beta_T_reservoir = as.numeric(mod$b[j_treat])
				if (ncol(as.matrix(private$X)) > 0 && !is.null(attempt$fit)){
					private$best_X_colnames_reservoir = setdiff(colnames(attempt$X), c("(Intercept)", "w_r"))
					private$best_X_j_treat_reservoir = j_treat
				}
			}
			if (!estimate_only && !is.null(mod)) private$cached_values$ssq_beta_T_reservoir = as.numeric(mod$vcov[j_treat, j_treat])
		},
		best_X_colnames_matched = NULL,
		best_X_colnames_reservoir = NULL,
		best_X_j_treat_reservoir = 2L
	))
)
#' Abstract class for Conditional Poisson Combined-Likelihood Compound Inference
#'
#' Fits a single joint likelihood over all KK design data for count responses.
#' The matched-pair component uses the conditional Poisson likelihood
#' (which simplifies to a binomial contribution per pair), and the reservoir
#' uses an ordinary Poisson likelihood.
#'
#' @keywords internal
InferenceAbstractKKPoissonCPoissonOneLik = R6::R6Class("InferenceAbstractKKPoissonCPoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = utils::modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
		},
		#' @description Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an asymptotic confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an asymptotic p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
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
		},
		#' @description Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			inf_obj = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			inf_obj
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Re-read design variables which might have been transformed during randomization
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y

			# Recompute basic match data for the new w/y
			private$compute_basic_match_data()

			# Clear combined_cov_keep to allow re-reduction if needed
			private$cached_values$combined_cov_keep = NULL

			# Use the same joint-likelihood logic for the point estimate
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		compute_treatment_estimate_during_power_simulation = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = 1e-4){
			private$compute_treatment_estimate_during_power_simulation_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
			assert_finite_se = function(){
				if (!is.finite(private$cached_values$s_beta_hat_T)){
					return(invisible(NULL))
				}
			},
			get_standard_error = function(){
				private$shared_combined_likelihood(estimate_only = FALSE)
				se = private$compute_standard_error_from_information_matrix()
				if (is.finite(se)) return(se)
				private$cached_values$s_beta_hat_T
			},
			get_degrees_of_freedom = function(){
				private$cached_values$df %||% NA_real_
			},
			supports_likelihood_tests = function(){
				!is.null(private$cached_values$likelihood_test_context)
			},
			supports_fisher_information = function(){
				TRUE
			},
			get_likelihood_test_spec = function(){
				private$shared_combined_likelihood(estimate_only = FALSE)
				ctx = private$cached_values$likelihood_test_context
				if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
				yT_v = as.numeric(ctx$yT_v)
				n_k_v = as.numeric(ctx$n_k_v)
				X_diff_v = as.matrix(ctx$X_diff_v)
				y_r_v = as.numeric(ctx$y_r_v)
				w_r_v = as.numeric(ctx$w_r_v)
				X_r_v = as.matrix(ctx$X_r_v)
				list(
					j = 2L,
					full_fit = private$cached_mod,
					fit_null = function(delta, start = NULL){
						n_params = ncol(X_diff_v) + 2L
						fast_cpoisson_combined_with_var_cpp(
							yT_v = yT_v,
							n_k_v = n_k_v,
							X_diff_v = X_diff_v,
							y_r = y_r_v,
							w_r = w_r_v,
							X_r = X_r_v,
							fixed_idx = 2L,
							fixed_values = delta,
							warm_start_params = start %||% private$get_fit_warm_start_for_length("params", n_params),
							warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
						)
					},
					score = function(fit){
						as.numeric(fit$score %||% get_cpoisson_combined_score_cpp(yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v, as.numeric(fit$params %||% fit$b)))
					},
					observed_information = function(fit){
						as.matrix(fit$observed_information %||% fit$information)
					},
					fisher_information = function(fit){
						as.matrix(fit$fisher_information %||% fit$information)
					},
					information = function(fit){
						as.matrix(fit$information %||% -get_cpoisson_combined_hessian_cpp(yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v, as.numeric(fit$params %||% fit$b)))
					},
					neg_loglik = function(fit){
						as.numeric(fit$neg_loglik %||% fit$neg_ll)
					}
				)
			},
		set_failed_combined_cache = function(){
			private$cached_values$beta_hat_T   = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cache_nonestimable_estimate("kk_cpoisson_combined_fit_failed")
		},
		reduce_combined_covariates = function(X_diff, X_r, w_r){
			# Use a heuristic to ensure X_diff and X_r are jointly full-rank when combined
			# with treatment and reservoir intercept.
			n_p = nrow(as.matrix(X_diff))
			n_r = length(w_r)
			p   = if (is.null(private$X)) 0L else ncol(as.matrix(private$X))
			
			# Layout: [beta_0, beta_T, beta_xs]
			# Pairs:     [0, 1, X_diff]
			# Reservoir: [1, w_r, X_r]
			
			# Handle empty components to avoid cbind recycling warnings
			pairs_part = if (n_p > 0) {
				cbind(Intercept = rep(0, n_p), treatment = rep(1, n_p), as.matrix(X_diff))
			} else {
				matrix(0, nrow = 0, ncol = p + 2)
			}
			
			reservoir_part = if (n_r > 0) {
				cbind(Intercept = rep(1, n_r), treatment = w_r, as.matrix(X_r))
			} else {
				matrix(0, nrow = 0, ncol = p + 2)
			}
			
			X_stack = rbind(pairs_part, reservoir_part)
			
			if (nrow(X_stack) == 0) return(integer(0))
			
			qr_res = qr(X_stack)
			if (is.finite(qr_res$rank) && qr_res$rank < ncol(X_stack)){
				keep = qr_res$pivot[seq_len(qr_res$rank)]
				# Always keep intercept(1) and treatment(2)
				if (!(1L %in% keep)) keep = c(1L, keep)
				if (!(2L %in% keep)) keep = c(2L, keep)
				keep = sort(unique(keep))
				# Extract original covariate indices (offset by 2)
				return(keep[keep > 2L] - 2L)
			}
			seq_len(p)
		},
		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
			try_combined_fit = function(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v){
				n_params = ncol(X_diff_v) + 2L
				mod = tryCatch(
					fast_cpoisson_combined_with_var_cpp(
						yT_v = as.numeric(yT_v),
						n_k_v = as.numeric(n_k_v),
						X_diff_v = as.matrix(X_diff_v),
						y_r = as.numeric(y_r_v),
						w_r = as.numeric(w_r_v),
						X_r = as.matrix(X_r_v),
						warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
						warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params)
					),
					error = function(e) NULL
				)
				if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
				if (!estimate_only){
					ssq = mod$ssq_b_j
					if (is.null(ssq) || !is.finite(ssq) || ssq < 0) return(FALSE)
				}
				private$cached_mod = mod
				private$set_fit_warm_start(as.numeric(mod$params), "params", fisher = mod$fisher_information)
				
				private$cached_values$likelihood_test_context = list(
					yT_v = as.numeric(yT_v),
					n_k_v = as.numeric(n_k_v),
					X_diff_v = as.matrix(X_diff_v),
					y_r_v = as.numeric(y_r_v),
					w_r_v = as.numeric(w_r_v),
					X_r_v = as.matrix(X_r_v)
				)
				private$cached_values$beta_hat_T = as.numeric(mod$b[2])
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
				private$cached_values$df = NA_real_
				TRUE
			},
		try_pairs_only = function(estimate_only, yT_v, n_k_v, X_diff_v){
			if (length(yT_v) == 0L) return(FALSE)
			y_prop = yT_v / n_k_v
			X = if (ncol(X_diff_v) > 0L) cbind(1, X_diff_v) else matrix(1, nrow = length(yT_v), ncol = 1L)
			mod = tryCatch({
				res = fast_logistic_regression_weighted_cpp(X = X, y = y_prop, weights = n_k_v)
				list(b = res$b, ssq_b_1 = NA_real_, X_fit = X) # weights version doesn't return var yet
			}, error = function(e) NULL)
			if (is.null(mod) || length(mod$b) < 1L || !is.finite(mod$b[1])) return(FALSE)
			private$cached_values$beta_hat_T = as.numeric(mod$b[1])
			TRUE
		},
		try_reservoir_only = function(estimate_only, y_r_v, w_r_v, X_r_v){
			if (length(y_r_v) == 0L) return(FALSE)
			X_full = if (ncol(X_r_v) > 0L) cbind(1, w_r_v, X_r_v) else cbind(1, w_r_v)
			X_fit = X_full
			j_treat = 2L
			if (ncol(X_full) > 2L){
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (is.finite(r_full) && r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(2L %in% keep)) {
						# Try to force treatment column in
						keep = c(2L, keep)
					}
					keep = sort(unique(keep))
					X_fit = X_full[, keep, drop = FALSE]
					j_treat = which(colnames(X_fit) == "w_r_v")
					if (length(j_treat) == 0) j_treat = 2L
				}
			}
			mod = tryCatch(
				fast_poisson_regression_with_var_cpp(X = X_fit, y = as.numeric(y_r_v), j = as.integer(j_treat[1])),
				error = function(e) NULL
			)
			if (is.null(mod) || length(j_treat) == 0 || length(mod$b) < j_treat[1] || !is.finite(mod$b[j_treat[1]])) return(FALSE)
			if (!estimate_only){
				ssq = mod$ssq_b_j
				if (is.null(ssq) || !is.finite(ssq) || ssq < 0) return(FALSE)
			}
			private$cached_values$beta_hat_T = as.numeric(mod$b[j_treat[1]])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			TRUE
		},
		# The combined case is handled by fast_cpoisson_combined_with_var_cpp
		# (Newton's method with analytic Fisher-information Hessian).
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$cached_values$likelihood_test_context = NULL
			private$cached_mod = NULL
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
				private$cached_values$combined_cov_keep = NULL
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			p             = if (is.null(private$X)) 0L else ncol(as.matrix(private$X))
			has_reservoir = nRT > 0 && nRC > 0
			# ---- Pair data (conditional Poisson, zero-total pairs discarded) ----
			yT_v     = numeric(0)
			n_k_v    = numeric(0)
			X_diff_v = matrix(nrow = 0L, ncol = p)
			if (m > 0){
				yT    = KKstats$yTs_matched
				yC    = KKstats$yCs_matched
				n_k   = yT + yC
				valid = which(n_k > 0)
				if (length(valid) > 0){
					yT_v  = yT[valid]
					n_k_v = n_k[valid]
					if (p > 0L) {
						# Use the full-width pair-difference matrix here so any later
						# covariate reduction stays aligned with the reservoir matrix.
						X_diff_v = as.matrix(KKstats$X_matched_diffs_full[valid, , drop = FALSE])
					}
				}
			}
			has_pairs = length(yT_v) > 0
			# ---- Reservoir data (marginal Poisson) ----
			y_r_v = numeric(0)
			w_r_v = numeric(0)
			X_r_v = matrix(nrow = 0L, ncol = p)
			if (has_reservoir){
				y_r_v = KKstats$y_reservoir
				w_r_v = KKstats$w_reservoir
				if (p > 0L) X_r_v = as.matrix(KKstats$X_reservoir)
			}
			if (!has_pairs && !has_reservoir){
				private$set_failed_combined_cache()
				return(invisible(NULL))
			}
			# ---- Covariate Reduction (Aligned) -----------------------------------
			if (p > 0L && is.null(private$cached_values$combined_cov_keep)){
				private$cached_values$combined_cov_keep = private$reduce_combined_covariates(X_diff_v, X_r_v, w_r_v)
			}
			if (p > 0L){
				keep = private$cached_values$combined_cov_keep
				if (length(keep) > 0L){
					X_diff_v = X_diff_v[, keep, drop = FALSE]
					X_r_v    = X_r_v[, keep, drop = FALSE]
				} else {
					X_diff_v = matrix(nrow = nrow(X_diff_v), ncol = 0L)
					X_r_v    = matrix(nrow = nrow(X_r_v), ncol = 0L)
				}
			}
			# ---- Joint Likelihood Fit --------------------------------------------
			success = FALSE
			if (has_pairs && has_reservoir){
				success = private$try_combined_fit(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v)
			}
			
			if (!success){
				# Fallback: one or the other
				fallback_success = FALSE
				if (has_pairs){
					fallback_success = private$try_pairs_only(estimate_only, yT_v, n_k_v, X_diff_v)
				}
				if (!fallback_success && has_reservoir){
					fallback_success = private$try_reservoir_only(estimate_only, y_r_v, w_r_v, X_r_v)
				}
				if (!fallback_success){
					private$set_failed_combined_cache()
				}
			}
			invisible(NULL)
		}
	))
)
#' KK Hurdle Poisson IVWC Inference for Count Responses
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with count
#' responses using a hurdle-Poisson mixed model for matched pairs and an ordinary
#' Poisson regression for reservoir subjects. Estimates are combined via
#' inverse-variance weighting.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountKKHurdlePoissonIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountKKHurdlePoissonIVWC = R6::R6Class("InferenceCountKKHurdlePoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKHurdlePoissonIVWC,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB} for
		#'   the matched-pair component.
		#' @param optimization_alg Optimization algorithm. Default is dispatched via policy.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			super$initialize(
				des_obj,
				model_formula = model_formula,
				use_rcpp = use_rcpp,
				optimization_alg = optimization_alg,
				verbose = verbose,
				smart_cold_start_default = smart_cold_start_default
			)
		}
	)
)
#' Conditional-Poisson Inference for KK Designs with IVWC
#'
#' Fits a conditional-Poisson regression for count responses under a KK design
#' using the Independent-Variables-as-Working-Covariates (IVWC) approach.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountKKCPoissonIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountKKCPoissonIVWC = R6::R6Class("InferenceCountKKCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonIVWC,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
		}
	)
)
#' Conditional-Poisson Inference for KK Designs with Combined Likelihood
#'
#' Fits a conditional-Poisson regression for count responses under a KK design
#' using the combined-likelihood approach.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'count')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rpois(10, 2))
#' inf = InferenceCountKKCPoissonOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceCountKKCPoissonOneLik = R6::R6Class("InferenceCountKKCPoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKPoissonCPoissonOneLik,
	public = list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed \code{DesignSeqOneByOneKK14} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, smart_cold_start_default = TRUE){
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose, smart_cold_start_default = smart_cold_start_default)
		}
	)
)
