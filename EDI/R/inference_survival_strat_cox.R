#' Stratified Cox PH Inference for Survival Responses
#'
#' Fits an auto-stratified Cox PH regression. Stratification variables are chosen
#' automatically from the recorded low-cardinality covariates. If no suitable
#' stratification covariates are found, the fit falls back to the corresponding
#' standard Cox PH model.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalStratCoxPHRegr$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceSurvivalStratCoxPHRegr = R6::R6Class("InferenceSurvivalStratCoxPHRegr",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with a survival response.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   covariates from the design object are included. Use \code{~ 1} for univariate.
		#' @param use_rcpp Logical. If \code{TRUE} (default), enable internal Rcpp score/information helpers for likelihood inference.
		#'   Cox optimization uses \pkg{survival::coxph.fit}.
		#' @param optimization_alg Optimization algorithm: \code{"newton_raphson"} (default) or \code{"lbfgs"}.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = "lbfgs", verbose = FALSE, smart_cold_start_default = NULL) {
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
				assertFlag(use_rcpp)
			}
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$use_rcpp = use_rcpp
		},
		#' @description Compute the treatment effect estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the stratified Cox PH treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) return(beta_hat_T)
			}
			X_cov = private$get_X()
			X_cov = private$reduce_covariates_preserving_treatment(X_cov)
			X_fit = if (!is.null(X_cov) && ncol(X_cov) > 0) cbind(treatment = private$w, X_cov) else matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			strata_info = private$compute_strata_info(X_cov)
			use_strata = isTRUE(strata_info$num_strata > 1L)
			fit = weighted_cox_bootstrap_surrogate_fit(private$y, private$dead, X_fit, row_weights, strata = if (use_strata) strata_info$strata_id else NULL)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Computes an approximate confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Compute confidence interval rand
		#' @param alpha The significance level (default 0.05).
		#' @param r Number of vectors to draw.
		#' @param pval_epsilon The bisection convergence tolerance.
		#' @param show_progress Whether to show a progress bar.
		#' @param ci_search_control Unused.
		compute_rand_confidence_interval = function(alpha = 0.05, r = 501, pval_epsilon = 0.005, show_progress = TRUE, ci_search_control = NULL){
			stop("Randomization confidence intervals are not supported for stratified Cox PH models because the estimator units (Log-Hazard Ratio) are inconsistent with the randomization test's required transformed scale (Log-Time Ratio / AFT effect).")
		}
	),
		private = list(
		use_rcpp = TRUE,
		cached_mod = NULL,
		coxph_control = NULL,
		get_complexity_tier = function() "light",
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			mod = private$generate_mod(estimate_only = estimate_only)
			private$cached_values$beta_hat_T = as.numeric(mod$b[2])
			if (estimate_only) return(invisible(NULL))
			se = if (is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0) sqrt(mod$ssq_b_2) else NA_real_
			private$cached_values$s_beta_hat_T = se
			private$cached_values$df = NA_real_
		},
		supports_likelihood_tests = function(){
			isTRUE(private$use_rcpp)
		},
		supports_lik_ratio_param_bootstrap = function() isTRUE(private$use_rcpp),
		compute_lik_ratio_confidence_interval_impl = function(alpha){
			est = self$compute_estimate()
			if (!is.finite(est)) return(c(NA_real_, NA_real_))
			private$shared(estimate_only = FALSE)
			se = private$cached_values$s_beta_hat_T
			step = if (is.finite(se) && se > 0) se else max(abs(est), 1)
			step = max(step, 1e-4)
			wald_ci = private$compute_wald_confidence_interval_impl(alpha)
			seeds = c(
				if (length(wald_ci) >= 1L && is.finite(wald_ci[[1L]])) wald_ci[[1L]] else est - step,
				if (length(wald_ci) >= 2L && is.finite(wald_ci[[2L]])) wald_ci[[2L]] else est + step
			)
			pval_fn = function(delta){
				tryCatch(self$compute_lik_ratio_two_sided_pval(delta), error = function(e) NA_real_)
			}
			find_bound = function(direction, seed){
				f_est = pval_fn(est) - alpha
				if (!is.finite(f_est)) return(NA_real_)
				if (f_est <= 0) return(est)
				outer = NA_real_
				if (is.finite(seed) && ((direction < 0 && seed < est) || (direction > 0 && seed > est))) {
					f_seed = pval_fn(seed) - alpha
					if (is.finite(f_seed) && f_seed <= 0) outer = seed
				}
				if (!is.finite(outer)) {
					for (i in 0:59) {
						d = est + direction * step * 2^i
						f_d = pval_fn(d) - alpha
						if (is.finite(f_d) && f_d <= 0) {
							outer = d
							break
						}
					}
				}
				if (!is.finite(outer)) return(NA_real_)
				root = tryCatch(
					stats::uniroot(function(d) pval_fn(d) - alpha, lower = min(est, outer), upper = max(est, outer), tol = 1e-6)$root,
					error = function(e) NA_real_
				)
				as.numeric(root)
			}
			ci = c(find_bound(-1, seeds[1]), find_bound(1, seeds[2]))
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$coefficients %||% null_fit$b)
			if (!all(is.finite(b_null))) return(NULL)
			X_fit = spec$X
			y_obs = as.numeric(spec$y)
			dead_obs = as.numeric(spec$dead)
			strata = as.integer(spec$strata)
			stratified = isTRUE(spec$stratified)
			j = spec$j
			if (stratified && !is.null(strata)){
				sim = .cox_simulate_stratified(y_obs, dead_obs, X_fit, b_null, strata)
			} else {
				breslow = .breslow_hazard(y_obs, dead_obs, X_fit, b_null)
				if (length(breslow$times) == 0L) return(NULL)
				sim = .cox_simulate_from_breslow(breslow, y_obs, dead_obs, X_fit, b_null)
			}
			y_sim = sim$y_sim; dead_sim = sim$dead_sim
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			strata_arg = if (stratified && !is.null(strata)) strata else NULL
			full_res = .fit_survival_coxph_kernel(X_fit, y_sim, dead_sim, strata = strata_arg)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$coefficients[j])) return(NULL)
			full_fit_boot = list(b = as.numeric(full_res$coefficients), neg_loglik = as.numeric(full_res$neg_ll))
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = .fit_survival_coxph_fixed_kernel(X_fit, y_sim, dead_sim, strata = strata_arg, fixed_idx = j, fixed_value = d)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), neg_loglik = as.numeric(res$neg_ll))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik)
			)
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			dead = as.numeric(ctx$dead)
			strata = as.integer(ctx$strata)
			stratified = isTRUE(ctx$stratified)
			j_treat = as.integer(ctx$j_treat %||% 1L)
			list(
				X = X_fit,
				y = y,
				dead = dead,
				strata = strata,
				stratified = stratified,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					strata_arg = if (stratified) strata else NULL
					.fit_survival_coxph_fixed_kernel(X_fit, y, dead, strata = strata_arg, fixed_idx = j_treat, fixed_value = delta)
				},
				extract_start = function(fit){
					as.numeric(fit$coefficients %||% fit$b)
				},
				score = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						.cox_score_breslow_fd_r(X_fit, y, dead, beta, strata = strata)
					} else {
						get_coxph_score_cpp(X_fit, y, dead, beta)
					}
				},
				observed_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						.cox_information_breslow_fd_r(X_fit, y, dead, beta, strata = strata)
					} else {
						-get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				fisher_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					fit$fisher_information %||% if (stratified) {
						.cox_information_breslow_fd_r(X_fit, y, dead, beta, strata = strata)
					} else {
						-get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					fit$information %||% fit$fisher_information %||% if (stratified) {
						.cox_information_breslow_fd_r(X_fit, y, dead, beta, strata = strata)
					} else {
						-get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_ll %||% fit$neg_loglik)
				}
			)
		},
		compute_strata_info = function(X_full) {
			n = length(private$y)
			if (is.null(X_full) || ncol(X_full) == 0){
				return(list(strata_id = rep.int(1L, n), selected_cols = integer(0), num_strata = 1L))
			}
			info = tryCatch(
				compute_survival_strata_ids_cpp(as.matrix(X_full)),
				error = function(e) NULL
			)
			if (is.null(info)){
				return(list(strata_id = rep.int(1L, n), selected_cols = integer(0), num_strata = 1L))
			}
			list(
				strata_id = as.integer(info$strata_id),
				selected_cols = as.integer(info$selected_cols),
				num_strata = as.integer(info$num_strata)
			)
		},
		reduce_covariates_preserving_treatment = function(X_covars){
			if (is.null(X_covars) || ncol(X_covars) == 0){
				return(matrix(numeric(0), nrow = length(private$y), ncol = 0))
			}
			X_covars = as.matrix(X_covars)
			if (ncol(X_covars) == 0){
				return(matrix(numeric(0), nrow = nrow(X_covars), ncol = 0))
			}
			full_design = cbind(w = private$w, X_covars)
			reduced = drop_linearly_dependent_cols(full_design)
			X_keep = reduced$M
			if (ncol(X_keep) == 0){
				return(matrix(numeric(0), nrow = nrow(full_design), ncol = 0))
			}
			if (!("w" %in% colnames(X_keep))){
				return(matrix(numeric(0), nrow = nrow(full_design), ncol = 0))
			}
			X_keep[, colnames(X_keep) != "w", drop = FALSE]
		},
		get_informative_rows = function(strata_id){
			if (length(strata_id) != length(private$y)) return(integer(0))
			good = rep(FALSE, length(strata_id))
			for (s in unique(strata_id)){
				i_s = which(strata_id == s)
				if (length(i_s) < 2) next
				if (length(unique(private$w[i_s])) < 2) next
				if (!any(private$dead[i_s] == 1, na.rm = TRUE)) next
				good[i_s] = TRUE
			}
			which(good)
		},
		fit_cox_with_formula = function(dat, formula_str){
			tryCatch(
				suppressWarnings(survival::coxph(stats::as.formula(formula_str), data = dat)),
				error = function(e) NULL
			)
		},
		format_mod_output = function(mod){
			if (is.null(mod)){
				return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_))
			}
			coef_w = tryCatch(as.numeric(stats::coef(mod)["w"]), error = function(e) NA_real_)
			ssq_w = tryCatch(as.numeric(stats::vcov(mod)["w", "w"]), error = function(e) NA_real_)
			list(
				b = c(0, coef_w),
				ssq_b_2 = if (is.finite(ssq_w) && ssq_w > 0) ssq_w else NA_real_
			)
		},
		# Build (y, dead, X_mat) for a given row set.  X_mat = cbind(w, X_linear).
		build_rcpp_inputs = function(rows, X_linear){
			y_r    = as.numeric(private$y[rows])
			dead_r = as.numeric(private$dead[rows])
			w_r    = private$w[rows]
			if (ncol(X_linear) > 0){
				X_mat = cbind(w = w_r, X_linear[rows, , drop = FALSE])
			} else {
				X_mat = matrix(w_r, ncol = 1L, dimnames = list(NULL, "w"))
			}
			list(
				y = y_r,
				dead = dead_r,
				X = as.matrix(X_mat),
				surv_y = survival::Surv(y_r, dead_r),
				rownames = as.character(seq_along(y_r))
			)
		},
		fit_coxph_estimate_only_fast = function(X, surv_y, strata = NULL, rownames = NULL){
			if (is.null(private$coxph_control)) {
				private$coxph_control = survival::coxph.control()
			}
			fit = tryCatch(
				survival::coxph.fit(
					x = X,
					y = surv_y,
					strata = if (is.null(strata)) NULL else as.integer(strata),
					offset = NULL,
					init = NULL,
					control = private$coxph_control,
					weights = NULL,
					method = "breslow",
					rownames = rownames %||% as.character(seq_len(nrow(X))),
					resid = FALSE
				),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NULL)
			b = as.numeric(fit$coefficients %||% numeric(0))
			if (length(b) != ncol(X) || !all(is.finite(b))) return(NULL)
			names(b) = colnames(X)
			list(b = b, coefficients = b, vcov = NULL, var = NULL,
				neg_ll = NA_real_, neg_loglik = NA_real_, neg_log_lik = NA_real_,
				fisher_information = NULL, converged = TRUE)
		},
		fit_rcpp_stratified = function(rows, X_linear, strata_id, estimate_only = FALSE){
			inp  = private$build_rcpp_inputs(rows, X_linear)
			strata_sub = as.integer(strata_id[rows])
			fit = if (estimate_only) {
				private$fit_coxph_estimate_only_fast(inp$X, inp$surv_y, strata = strata_sub, rownames = inp$rownames)
			} else {
				.fit_survival_coxph_kernel(inp$X, inp$y, inp$dead, strata = strata_sub, estimate_only = FALSE)
			}
			if (is.null(fit)) return(NULL)
			if (estimate_only) return(list(fit = fit, stratified = TRUE))
			list(fit = fit, X = inp$X, y = inp$y, dead = inp$dead, strata = strata_sub, stratified = TRUE)
		},

		fit_rcpp_unstratified = function(rows, X_linear, estimate_only = FALSE){
			inp = private$build_rcpp_inputs(rows, X_linear)
			fit = if (estimate_only) {
				private$fit_coxph_estimate_only_fast(inp$X, inp$surv_y, rownames = inp$rownames)
			} else {
				.fit_survival_coxph_kernel(inp$X, inp$y, inp$dead, estimate_only = FALSE)
			}
			if (is.null(fit)) return(NULL)
			if (estimate_only) return(list(fit = fit, stratified = FALSE))
			list(fit = fit, X = inp$X, y = inp$y, dead = inp$dead, strata = NULL, stratified = FALSE)
		},

		format_rcpp_output = function(fit){
			if (is.null(fit) || !isTRUE(fit$converged)) return(list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_))
			beta_w  = as.numeric(fit$coefficients[1])
			ssq_w   = if (!is.null(fit$vcov) && nrow(fit$vcov) >= 1L && is.finite(fit$vcov[1, 1]) && fit$vcov[1, 1] > 0)
				fit$vcov[1, 1] else NA_real_
			list(b = c(0, beta_w), ssq_b_2 = ssq_w, fisher_information = fit$fisher_information)
		},
		generate_mod = function(estimate_only = FALSE){
			X_full      = private$X
			strata_info = private$compute_strata_info(X_full)
			X_linear = matrix(numeric(0), nrow = length(private$y), ncol = 0)
			if (ncol(as.matrix(private$X)) > 0 && !is.null(X_full) && ncol(X_full) > 0){
				keep_cols = setdiff(seq_len(ncol(X_full)), strata_info$selected_cols)
				if (length(keep_cols) > 0){
					X_linear = private$reduce_covariates_preserving_treatment(X_full[, keep_cols, drop = FALSE])
				}
			}
			informative_rows = integer(0)
			if (!is.null(strata_info$strata_id) && isTRUE(strata_info$num_strata > 1L)){
				informative_rows = private$get_informative_rows(strata_info$strata_id)
			}
			if (length(informative_rows) >= 4){
				if (private$use_rcpp){
					res = private$fit_rcpp_stratified(informative_rows, X_linear, strata_info$strata_id, estimate_only = estimate_only)
					if (!is.null(res) && isTRUE(res$fit$converged)){
						private$cached_mod = res$fit
						if (!estimate_only) {
							private$set_fit_warm_start(as.numeric(res$fit$coefficients), "params", fisher = res$fit$fisher_information)
							private$cached_values$likelihood_test_context = list(
								X = res$X, y = res$y, dead = res$dead, strata = res$strata,
								stratified = TRUE, j_treat = 1L
							)
						}
						return(private$format_rcpp_output(res$fit))
					}
					res = private$fit_rcpp_unstratified(informative_rows, X_linear, estimate_only = estimate_only)
					if (!is.null(res) && isTRUE(res$fit$converged)){
						private$cached_mod = res$fit
						if (!estimate_only) {
							private$set_fit_warm_start(as.numeric(res$fit$coefficients), "params", fisher = res$fit$fisher_information)
							private$cached_values$likelihood_test_context = list(
								X = res$X, y = res$y, dead = res$dead, strata = res$strata,
								stratified = FALSE, j_treat = 1L
							)
						}
						return(private$format_rcpp_output(res$fit))
					}
				} else {
					colnames(X_linear) = paste0("x", seq_len(ncol(X_linear)))
					dat_full = data.frame(y = private$y, dead = private$dead, w = private$w)
					if (ncol(X_linear) > 0) dat_full = cbind(dat_full, as.data.frame(X_linear))
					base_terms   = c("w", setdiff(colnames(dat_full), c("y", "dead", "w")))
					base_formula = paste("survival::Surv(y, dead) ~", paste(base_terms, collapse = " + "))
					dat_strat    = dat_full[informative_rows, , drop = FALSE]
					dat_strat$strata_id = factor(strata_info$strata_id[informative_rows])
					mod = private$fit_cox_with_formula(dat_strat, paste(base_formula, "+ strata(strata_id)"))
					if (!is.null(mod)){
						private$cached_mod = mod
						private$cached_values$likelihood_test_context = list(
							X = as.matrix(cbind(w = private$w[informative_rows], X_linear[informative_rows, , drop = FALSE])),
							y = private$y[informative_rows],
							dead = private$dead[informative_rows],
							strata = strata_info$strata_id[informative_rows],
							stratified = TRUE,
							j_treat = 1L
						)
						return(private$format_mod_output(mod))
					}
					if (ncol(X_linear) > 0){
						mod = private$fit_cox_with_formula(
							dat_strat[, c("y", "dead", "w", "strata_id"), drop = FALSE],
							"survival::Surv(y, dead) ~ w + strata(strata_id)"
						)
						if (!is.null(mod)){
							private$cached_mod = mod
							private$cached_values$likelihood_test_context = list(
								X = as.matrix(cbind(w = private$w[informative_rows])),
								y = private$y[informative_rows],
								dead = private$dead[informative_rows],
								strata = strata_info$strata_id[informative_rows],
								stratified = TRUE,
								j_treat = 1L
							)
							return(private$format_mod_output(mod))
						}
					}
				}
			}
			all_rows = seq_len(length(private$y))
			if (private$use_rcpp){
				res = private$fit_rcpp_unstratified(all_rows, X_linear, estimate_only = estimate_only)
				if (!is.null(res) && isTRUE(res$fit$converged)){
					private$cached_mod = res$fit
					if (!estimate_only) {
						private$cached_values$likelihood_test_context = list(
							X = res$X, y = res$y, dead = res$dead, strata = NULL,
							stratified = FALSE, j_treat = 1L
						)
					}
					return(private$format_rcpp_output(res$fit))
				}
			} else {
				colnames(X_linear) = paste0("x", seq_len(ncol(X_linear)))
				dat_full = data.frame(y = private$y, dead = private$dead, w = private$w)
				if (ncol(X_linear) > 0) dat_full = cbind(dat_full, as.data.frame(X_linear))
				base_terms   = c("w", setdiff(colnames(dat_full), c("y", "dead", "w")))
				base_formula = paste("survival::Surv(y, dead) ~", paste(base_terms, collapse = " + "))
				mod = private$fit_cox_with_formula(dat_full, base_formula)
				if (is.null(mod) && ncol(X_linear) > 0){
					mod = private$fit_cox_with_formula(
						dat_full[, c("y", "dead", "w"), drop = FALSE],
						"survival::Surv(y, dead) ~ w"
					)
				}
				if (!is.null(mod)){
					private$cached_mod = mod
					private$cached_values$likelihood_test_context = list(
						X = as.matrix(cbind(w = private$w, X_linear)),
						y = private$y,
						dead = private$dead,
						strata = NULL,
						stratified = FALSE,
						j_treat = 1L
					)
				}
				return(private$format_mod_output(mod))
			}
			list(b = c(NA_real_, NA_real_), ssq_b_2 = NA_real_)
		}
	)
)
