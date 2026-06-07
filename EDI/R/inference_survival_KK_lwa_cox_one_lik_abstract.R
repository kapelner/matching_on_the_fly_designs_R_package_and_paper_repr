#' Abstract class for LWA-style Marginal Cox Combined-Likelihood Inference
#'
#' Fits a single joint marginal Cox model over all KK design data for survival
#' responses. Matched subjects share their pair ID as a cluster, and reservoir
#' subjects are treated as independent (unique clusters). Standard errors are
#' obtained via the Huber-White cluster-robust sandwich estimator (LWA style).
#'
#' @keywords internal
InferenceAbstractKKLWACoxOneLik = R6::R6Class("InferenceAbstractKKLWACoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceParamBootstrap,
	public = utils::modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Recomputes the LWA one-likelihood treatment estimate under
		#'   Bayesian-bootstrap weights.
		#' @param subject_or_block_weights Subject-, block-, cluster-, or matched-set
		#'   bootstrap weights.
		#' @param estimate_only If \code{TRUE}, compute only the weighted point
		#'   estimate.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE){
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(self$compute_estimate(estimate_only = TRUE))[1L]
				if (is.finite(beta_hat_T)) {
					private$cached_values$beta_hat_T = beta_hat_T
					private$cached_values$s_beta_hat_T = NA_real_
					return(private$cached_values$beta_hat_T)
				}
			}
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			cluster_ids = m_vec
			res_idx = which(cluster_ids == 0L)
			if (length(res_idx) > 0L){
				max_m = max(cluster_ids)
				cluster_ids[res_idx] = max_m + seq_along(res_idx)
			}
			X_fit = private$design_matrix_candidates()
			fit = weighted_cox_bootstrap_surrogate_fit(
				private$y, private$dead, X_fit, row_weights,
				cluster = cluster_ids,
				warm_start_beta = private$get_fit_warm_start_for_length("params", ncol(X_fit)) %||% private$get_fit_warm_start_for_length("beta", ncol(X_fit))
			)
			private$cached_values$beta_hat_T = if (is.null(fit)) NA_real_ else as.numeric(fit$beta_hat)
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$beta_hat_T
		},
		#' @description Compute an asymptotic confidence interval.
		#' @param alpha Significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Compute an asymptotic two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
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
			eval(body(InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T))
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		max_abs_reasonable_coef = 1e4,
		optimization_alg = "lbfgs",
		get_standard_error = function(){
			private$shared_combined_likelihood()
			as.numeric(private$cached_values$s_beta_hat_T)
		},
		get_degrees_of_freedom = function() Inf,
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
		supports_likelihood_tests = function() TRUE,
		supports_lik_ratio_param_bootstrap = function() TRUE,
		simulate_under_lik_null = function(spec, delta, null_fit){
			b_null = as.numeric(null_fit$coefficients %||% null_fit$b)
			if (!all(is.finite(b_null))) return(NULL)
			X_fit = spec$X
			y_obs = as.numeric(spec$y)
			dead_obs = as.numeric(spec$dead)
			group_id = as.integer(spec$group_id)
			j = spec$j
			breslow = .breslow_hazard(y_obs, dead_obs, X_fit, b_null)
			if (length(breslow$times) == 0L) return(NULL)
			sim = .cox_simulate_from_breslow(breslow, y_obs, dead_obs, X_fit, b_null)
			y_sim = sim$y_sim; dead_sim = sim$dead_sim
			if (!all(is.finite(y_sim)) || any(y_sim <= 0)) return(NULL)
			full_res = tryCatch(
				fast_coxph_regression_cpp(
					X = X_fit, y = y_sim, dead = dead_sim,
					cluster = group_id,
					estimate_only = FALSE,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(full_res) || !isTRUE(full_res$converged) || !is.finite(full_res$coefficients[j])) return(NULL)
			full_fit_boot = list(b = as.numeric(full_res$coefficients), neg_loglik = as.numeric(full_res$neg_ll))
			list(
				full_fit = full_fit_boot,
				fit_null = function(d, start = NULL){
					res = tryCatch(
						fast_coxph_regression_cpp(
							X = X_fit, y = y_sim, dead = dead_sim,
							cluster = group_id,
							warm_start_beta = start %||% as.numeric(full_res$coefficients),
							fixed_idx = j, fixed_values = d,
							estimate_only = FALSE,
							optimization_alg = private$optimization_alg
						),
						error = function(e) NULL
					)
					if (is.null(res) || !isTRUE(res$converged)) return(NULL)
					list(b = as.numeric(res$coefficients), neg_loglik = as.numeric(res$neg_ll))
				},
				neg_loglik = function(fit) as.numeric(fit$neg_loglik %||% fit$neg_ll)
			)
		},
		get_likelihood_test_spec = function(){
			private$shared_combined_likelihood(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X_fit = ctx$X
			y = as.numeric(ctx$y)
			dead = as.numeric(ctx$dead)
			group_id = as.integer(ctx$group_id)
			j_treat = 1L
			list(
				X = X_fit,
				y = y,
				dead = dead,
				group_id = group_id,
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					warm_start_beta = start %||% private$get_fit_warm_start_for_length("params", ncol(X_fit))
					fast_coxph_regression_cpp(
						X = X_fit,
						y = y,
						dead = dead,
						cluster = group_id,
						warm_start_beta = warm_start_beta,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg,
						fixed_idx = j_treat,
						fixed_values = delta
				)
			},
				extract_start = function(fit){
					as.numeric(fit$coefficients %||% fit$b)
				},
				score = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					get_coxph_score_cpp(X_fit, y, dead, beta)
				},
				observed_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					-get_coxph_hessian_cpp(X_fit, y, dead, beta)
				},
				fisher_information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					-get_coxph_hessian_cpp(X_fit, y, dead, beta)
				},
				information = function(fit){
					beta = as.numeric(fit$coefficients %||% fit$b)
					-get_coxph_hessian_cpp(X_fit, y, dead, beta)
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_ll %||% fit$neg_loglik)
				}
			)
		},
		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Re-read w, y, dead because they might have been transformed for randomization
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			
			# Recompute basic match data for the new w/y/dead
			private$compute_basic_match_data()
			
			# Ensure we have the best design from the original data
			if (is.null(private$cached_values$best_X_colnames)){
				private$shared_combined_likelihood(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$cached_values$best_X_colnames)){
				return(NA_real_)
			}
			X_data = private$get_X()
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			X_covs_filtered = X_data[, intersect(private$cached_values$best_X_colnames, colnames(X_data)), drop = FALSE]
			if (ncol(X_covs_filtered) > 0){
				X_full = cbind(X_full, X_covs_filtered)
			}
			
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			# Clusters: pairs + reservoir
			cluster_ids = m_vec
			res_idx = which(cluster_ids == 0L)
			if (length(res_idx) > 0L){
				max_m = if (any(cluster_ids > 0)) max(cluster_ids) else 0L
				cluster_ids[res_idx] = max_m + seq_along(res_idx)
			}
			fit = tryCatch(
				fast_coxph_regression_cpp(
					X = as.matrix(X_full),
					y = private$y,
					dead = private$dead,
					cluster = as.integer(cluster_ids),
					estimate_only = estimate_only,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || length(fit$coefficients) < 1 || !is.finite(fit$coefficients[1])) return(NA_real_)
			as.numeric(fit$coefficients[1])
		},
		design_matrix_candidates = function(){
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			X_covs = private$get_X()
			if (ncol(as.matrix(X_covs)) > 0){
				X_full = cbind(X_full, as.matrix(X_covs))
			}
			X_full
		},
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			# Clusters: pairs + reservoir
			cluster_ids = m_vec
			res_idx = which(cluster_ids == 0L)
			if (length(res_idx) > 0L){
				max_m = max(cluster_ids)
				cluster_ids[res_idx] = max_m + seq_along(res_idx)
			}
			X_full = private$design_matrix_candidates()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L, # treatment
				fit_fun = function(X_fit){
					fast_coxph_regression_cpp(
						X = as.matrix(X_fit),
						y = private$y,
						dead = private$dead,
						cluster = as.integer(cluster_ids),
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					)
				},
				fit_ok = function(res, X_fit, keep){
					if (is.null(res) || !isTRUE(res$converged)) return(FALSE)
					beta = res$coefficients[1L]
					if (!is.finite(beta) || abs(beta) > private$max_abs_reasonable_coef) return(FALSE)
					if (estimate_only) return(TRUE)
					se   = tryCatch(sqrt(res$vcov[1L, 1L]), error = function(e) NA_real_)
					is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
				}
			)
				res = attempt$fit
				if (!is.null(res)){
					private$cached_values$best_X_colnames = setdiff(colnames(attempt$X), "w")
					private$cached_mod = res
					private$cached_values$likelihood_test_context = list(
						X = attempt$X,
						y = private$y,
						dead = private$dead,
						group_id = cluster_ids,
						j_treat = 1L
					)
					private$cached_values$beta_hat_T = as.numeric(res$coefficients[1])
					if (!estimate_only) {
						se = sqrt(res$vcov[1L, 1L])
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
				return(invisible(NULL))
			}
			private$cache_nonestimable_estimate("kk_lwa_cox_combined_fit_failed")
			invisible(NULL)
		}
	))
)
