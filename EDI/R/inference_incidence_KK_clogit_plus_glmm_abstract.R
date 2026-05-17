#' Abstract Conditional Logistic Plus GLMM Inference
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs and reservoir subjects.
#'
#' @keywords internal
InferenceAbstractKKClogitPlusGLMM = R6::R6Class("InferenceAbstractKKClogitPlusGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = utils::modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object with an incidence response.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param max_abs_reasonable_coef Cap for reasonable coefficient estimates.
		#' @param max_abs_log_sigma Cap for reasonable log random effect variance.
		#' @param verbose Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart optimizer start values.
		initialize = function(des_obj, model_formula = NULL, max_abs_reasonable_coef = 1e4, max_abs_log_sigma = 8, verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$max_abs_reasonable_coef = max_abs_reasonable_coef
			private$max_abs_log_sigma = max_abs_log_sigma
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$init_kk_passthrough(des_obj)
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
			private$shared(estimate_only = FALSE)
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		#' @description Computes an approximate two-sided p-value.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval = function(delta = 0){
			private$shared(estimate_only = FALSE)
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
		}
	)),
	private = utils::modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		max_abs_reasonable_coef = 1e4,
		max_abs_log_sigma = 8,
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && isTRUE(private$cached_values$s_beta_hat_T > 0)) return(invisible(NULL))
			private$clear_nonestimable_state()
			private$cached_mod = NULL
			private$cached_values$likelihood_test_context = NULL
			
			d = private$prepare_clogit_plus_glmm_data()
			if (!d$has_discordant && !d$has_concordant) {
				private$cache_nonestimable_estimate("no_data_for_clogit_plus_glmm")
				return(invisible(NULL))
			}
			
			n_params = ncol(d$X_conc) + 1L # betas + log_sigma
			fit = tryCatch(
				fast_clogit_plus_glmm_cpp(
					X_disc = d$X_disc, y_disc = d$y_disc,
					X_conc = d$X_conc, y_conc = d$y_conc,
					group_conc = d$group_conc,
					has_discordant = d$has_discordant,
					has_concordant = d$has_concordant,
					warm_start_params = private$get_fit_warm_start_for_length("params", n_params),
					warm_start_fisher_info = private$get_fit_warm_start_fisher(n_params),
					estimate_only = estimate_only,
					max_abs_log_sigma = private$max_abs_log_sigma,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) {
				private$cache_nonestimable_estimate("joint_likelihood_failed_to_converge")
				return(invisible(NULL))
			}
			
			# Treatment effect index:
			# If has_concordant is true, shared beta starts at par[0] (intercept), so treatment is par[1].
			# If has_concordant is false, shared beta starts at par[0], so treatment is par[0].
			j_T = if (d$has_concordant) 2L else 1L
			
			beta_hat_T = as.numeric(fit$params[j_T])
			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				private$cache_nonestimable_estimate("joint_likelihood_nonestimable")
				return(invisible(NULL))
			}
			
			private$cached_mod = fit
			private$set_fit_warm_start(as.numeric(fit$params), "params", fisher = fit$fisher_information)
			private$cached_values$likelihood_test_context = list(
				d = d,
				j_T = j_T,
				start = as.numeric(fit$params)
			)
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df   = Inf
			if (estimate_only) return(invisible(NULL))
			ssq = fit$ssq_b_j 
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		},
		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			d = ctx$d
			j_treat = as.integer(ctx$j_T)
			list(
				X = d$X_conc, 
				j = j_treat,
				full_fit = private$cached_mod,
				fit_null = function(delta, start = NULL){
					fast_clogit_plus_glmm_cpp(
						X_disc = d$X_disc, y_disc = d$y_disc,
						X_conc = d$X_conc, y_conc = d$y_conc,
						group_conc = d$group_conc,
						has_discordant = d$has_discordant,
						has_concordant = d$has_concordant,
						warm_start_params = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% ctx$start,
						warm_start_fisher_info = private$get_fit_warm_start_fisher(length(ctx$start)),
						estimate_only = FALSE,
						max_abs_log_sigma = private$max_abs_log_sigma,
						fixed_idx = j_treat, fixed_values = delta,
						optimization_alg = private$optimization_alg
					)
				},
				extract_start = function(fit){ as.numeric(fit$params) },
				score = function(fit){
					as.numeric(get_clogit_plus_glmm_score_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				observed_information = function(fit){
					as.matrix(get_clogit_plus_glmm_hessian_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				fisher_information = function(fit){
					as.matrix(get_clogit_plus_glmm_hessian_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				information = function(fit){
					as.matrix(get_clogit_plus_glmm_hessian_cpp(
						d$X_disc, d$y_disc, d$X_conc, d$y_conc, d$group_conc,
						as.numeric(fit$params), d$has_discordant, d$has_concordant,
						private$max_abs_log_sigma
					))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll)
				}
			)
		},
		prepare_clogit_plus_glmm_data = function(){
			private$compute_basic_match_data()
			KKstats = private$cached_values$KKstats
			
			# Discordant pairs
			yTs = as.numeric(KKstats$yTs_matched)
			yCs = as.numeric(KKstats$yCs_matched)
			y_m = yTs - yCs # 1 if (1,0), -1 if (0,1)
			
			i_m_disc = which(abs(y_m) == 1)
			X_disc_cov = KKstats$X_matched_diffs_full[i_m_disc, , drop = FALSE]
			# Prepend 1 for the treatment effect in clogit part
			X_disc = cbind(treatment = 1, X_disc_cov)
			# y for clogit should be 1 if (1,0) and 0 if (0,1)
			y_disc = (y_m[i_m_disc] + 1) / 2
			
			# Concordant pairs and reservoir for GLMM part
			m_vec = private$m
			i_m_conc = which(y_m == 0)
			
			i_conc_pair = which(m_vec %in% i_m_conc)
			i_res = which(is.na(m_vec) | m_vec == 0L)
			
			i_glmm = if (private$combine_reservoir_into_glmm()) c(i_conc_pair, i_res) else i_conc_pair
			
			X_glmm_cov = private$get_X()[i_glmm, , drop = FALSE]
			y_glmm = private$y[i_glmm]
			w_glmm = private$w[i_glmm]
			group_glmm = m_vec[i_glmm]
			i_glmm_res = which(is.na(group_glmm) | group_glmm == 0L)
			if (length(i_glmm_res) > 0L) {
				group_glmm[i_glmm_res] = max(c(0L, m_vec), na.rm = TRUE) + seq_along(i_glmm_res)
			}
			
			X_glmm = cbind(`(Intercept)` = 1, treatment = w_glmm, X_glmm_cov)
			
			list(
				X_disc = as.matrix(X_disc),
				y_disc = as.numeric(y_disc),
				X_conc = as.matrix(X_glmm),
				y_conc = as.numeric(y_glmm),
				group_conc = as.integer(group_glmm),
				has_discordant = length(i_m_disc) > 0L,
				has_concordant = length(i_glmm) > 0L
			)
		},
		combine_reservoir_into_glmm = function() stop("must implement combine_reservoir_into_glmm()"),
		log_sum_exp = function(x){
			m = max(x)
			if (!is.finite(m)) return(m)
			m + log(sum(exp(x - m)))
		},
		log1pexp = function(x){
			ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
		}
	))
)
