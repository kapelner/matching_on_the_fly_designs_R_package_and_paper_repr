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
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use our internal Rcpp
		#'   implementations where available. If \code{FALSE}, use \pkg{glmmTMB}.
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
			private$shared_combined_hurdle(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("KK hurdle-Poisson combined-likelihood: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_confidence_interval(alpha = alpha))
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
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared_combined_hurdle()
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				warning("KK hurdle-Poisson combined-likelihood: falling back to bootstrap because standard error is unavailable.")
				return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		use_rcpp = TRUE,
			max_abs_reasonable_coef = 1e4,

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
					fit_null = function(delta){
						fast_hurdle_poisson_glmm_cpp(
							X = X_fit,
							y = y,
							group_id = group_id,
							j_T = j_T,
							estimate_only = FALSE,
							n_gh = n_gh,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
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
	)
)
