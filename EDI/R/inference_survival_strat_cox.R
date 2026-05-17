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
		#' @param use_rcpp Logical. If \code{TRUE} (default), use internal Rcpp Cox PH optimiser.
		#'   If \code{FALSE}, use \pkg{survival::coxph}.
		#' @param optimization_alg Optimization algorithm: \code{"newton_raphson"} (default) or \code{"lbfgs"}.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, optimization_alg = "newton_raphson", verbose = FALSE, smart_cold_start_default = TRUE) {
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
		get_complexity_tier = function() "light",
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			mod = private$generate_mod()
			private$cached_values$beta_hat_T = as.numeric(mod$b[2])
			if (estimate_only) return(invisible(NULL))
			se = if (is.finite(mod$ssq_b_2) && mod$ssq_b_2 > 0) sqrt(mod$ssq_b_2) else NA_real_
			private$cached_values$s_beta_hat_T = se
			private$cached_values$df = NA_real_
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
			dead = as.numeric(ctx$dead)
			strata = as.integer(ctx$strata)
			stratified = isTRUE(ctx$stratified)
			j_treat = as.integer(ctx$j_treat %||% 1L)
			list(
				X = X_fit,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					ws_args = private$get_backend_warm_start_args(ncol(X_fit))
					warm_start_beta = start %||% ws_args$warm_start_beta
					if (stratified) {
						fast_stratified_coxph_regression_cpp(
							X = X_fit,
							y = y,
							dead = dead,
							strata = strata,
							warm_start_beta = warm_start_beta,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					} else {
						fast_coxph_regression_cpp(
							X = X_fit,
							y = y,
							dead = dead,
							warm_start_beta = warm_start_beta,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg,
							fixed_idx = j_treat,
							fixed_values = delta
						)
					}
				},
				extract_start = function(fit){
					as.numeric(fit$coefficients %||% fit$b)
				},
				score = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						get_stratified_coxph_score_cpp(X_fit, y, dead, strata, beta)
					} else {
						get_coxph_score_cpp(X_fit, y, dead, beta)
					}
				},
				observed_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						-get_stratified_coxph_hessian_cpp(X_fit, y, dead, strata, beta)
					} else {
						-get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				fisher_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						-get_stratified_coxph_hessian_cpp(X_fit, y, dead, strata, beta)
					} else {
						-get_coxph_hessian_cpp(X_fit, y, dead, beta)
					}
				},
				information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					if (stratified) {
						-get_stratified_coxph_hessian_cpp(X_fit, y, dead, strata, beta)
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
			list(y = y_r, dead = dead_r, X = as.matrix(X_mat))
		},
		fit_rcpp_stratified = function(rows, X_linear, strata_id){
			inp  = private$build_rcpp_inputs(rows, X_linear)
			strata_sub = as.integer(strata_id[rows])
			fit = tryCatch(
				fast_stratified_coxph_regression_cpp(
					X             = inp$X,
					y             = inp$y,
					dead          = inp$dead,
					strata        = strata_sub,
					warm_start_beta    = private$get_fit_warm_start_for_length("params", ncol(inp$X)),
					smart_cold_start   = private$smart_cold_start_default,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(inp$X)),
					estimate_only = FALSE,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NULL)
			list(fit = fit, X = inp$X, y = inp$y, dead = inp$dead, strata = strata_sub, stratified = TRUE)
		},

		fit_rcpp_unstratified = function(rows, X_linear){
			inp = private$build_rcpp_inputs(rows, X_linear)
			fit = tryCatch(
				fast_coxph_regression_cpp(
					X             = inp$X,
					y             = inp$y,
					dead          = inp$dead,
					warm_start_beta    = private$get_fit_warm_start_for_length("params", ncol(inp$X)),
					smart_cold_start   = private$smart_cold_start_default,
					warm_start_fisher_info = private$get_fit_warm_start_fisher(ncol(inp$X)),
					estimate_only = FALSE,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NULL)
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
					res = private$fit_rcpp_stratified(informative_rows, X_linear, strata_info$strata_id)
					if (!is.null(res) && isTRUE(res$fit$converged)){
						private$cached_mod = res$fit
						private$set_fit_warm_start(as.numeric(res$fit$coefficients), "params", fisher = res$fit$fisher_information)
						private$cached_values$likelihood_test_context = list(
							X = res$X,
							y = res$y,
							dead = res$dead,
							strata = res$strata,
							stratified = TRUE,
							j_treat = 1L
						)
						return(private$format_rcpp_output(res$fit))
					}
					res = private$fit_rcpp_unstratified(informative_rows, X_linear)
					if (!is.null(res) && isTRUE(res$fit$converged)){
						private$cached_mod = res$fit
						private$set_fit_warm_start(as.numeric(res$fit$coefficients), "params", fisher = res$fit$fisher_information)
						private$cached_values$likelihood_test_context = list(
							X = res$X,
							y = res$y,
							dead = res$dead,
							strata = res$strata,
							stratified = FALSE,
							j_treat = 1L
						)
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
				res = private$fit_rcpp_unstratified(all_rows, X_linear)
				if (!is.null(res) && isTRUE(res$fit$converged)){
					private$cached_mod = res$fit
					private$cached_values$likelihood_test_context = list(
						X = res$X,
						y = res$y,
						dead = res$dead,
						strata = NULL,
						stratified = FALSE,
						j_treat = 1L
					)
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
