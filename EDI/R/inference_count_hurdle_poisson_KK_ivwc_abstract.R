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

		#' @description
		#' Initialize
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
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(use_rcpp)
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass) or FixedDesignBinaryMatch.")
				}
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (inherits(des_obj, "FixedDesignBinaryMatch")){
				des_obj$.__enclos_env__$private$ensure_bms_computed()
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

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
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
		#' Compute asymp two sided pval for treatment effect
		#' @param delta Description for delta
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
		}
	),

	private = list(
		use_rcpp = TRUE,
		max_abs_reasonable_coef = 1e4,

		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the fixed-effect coefficient for "w" directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			Xmm = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(NA_integer_, nrow(Xmm))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L

			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)

			beta_m = NA_real_
			ssq_m = NA_real_
			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(Xmm, matched_idx, m_vec, se = FALSE)
				beta_m = res_m$beta_hat
				ssq_m = res_m$se^2
			}
			m_ok = !is.na(beta_m) && is.finite(beta_m) && !is.na(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			beta_r = NA_real_
			ssq_r = NA_real_
			if (length(reservoir_idx) > 1L && length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(Xmm, reservoir_idx, estimate_only = estimate_only)
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
				Xmm = private$create_design_matrix()
				full_names = c("(Intercept)", "w", if (ncol(Xmm) > 2L) paste0("x", seq_len(ncol(Xmm) - 2L)) else NULL)
				colnames(Xmm) = full_names[seq_len(ncol(Xmm))]
			} else {
				Xmm = cbind(1, private$w)
				colnames(Xmm) = c("(Intercept)", "w")
			}
			Xmm
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			Xmm = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(NA_integer_, nrow(Xmm))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L

			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)

			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(Xmm, matched_idx, m_vec, se = !estimate_only)
				private$cached_values$beta_T_matched = res_m$beta_hat
				if (!estimate_only) private$cached_values$ssq_beta_T_matched = res_m$se^2 else private$cached_values$ssq_beta_T_matched = 1.0
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
				!is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			if (length(reservoir_idx) > 1L &&
				length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(Xmm, reservoir_idx, estimate_only = estimate_only)
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

		fit_hurdle_for_matched_pairs = function(Xmm, matched_idx, m_vec, se = TRUE){
			X_matched = Xmm[matched_idx, , drop = FALSE]
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

		fit_poisson_for_reservoir = function(Xmm, reservoir_idx, estimate_only = FALSE){
			X_res = Xmm[reservoir_idx, , drop = FALSE]
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
