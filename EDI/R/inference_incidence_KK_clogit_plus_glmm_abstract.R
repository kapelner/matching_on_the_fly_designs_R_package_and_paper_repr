#'  internalAbstract Conditional Logistic Plus GLMM Inference
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The fixed effects are shared by both likelihood components. The
#' reservoir component is also included as additional independent observations
#' in the GLMM part.
#'
#' @keywords internal
InferenceAbstractKKClogitPlusGLMM = R6::R6Class("InferenceAbstractKKClogitPlusGLMM",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param max_abs_reasonable_coef Cap for coefficient magnitudes.
		#' @param max_abs_log_sigma Cap for log-random-effect-SD.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, max_abs_reasonable_coef = 1e4, max_abs_log_sigma = 8, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			private$max_abs_reasonable_coef = max_abs_reasonable_coef
			private$max_abs_log_sigma = max_abs_log_sigma
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha Confidence level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_confidence_interval(alpha = alpha))
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta Null difference.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			if (!identical(self$get_testing_type(), "wald")) {
				return(super$compute_asymp_two_sided_pval(delta = delta))
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
			max_abs_reasonable_coef = 1e4,
			max_abs_log_sigma = 8,

			get_standard_error = function(){
				private$shared(estimate_only = FALSE)
				se = private$compute_standard_error_from_information_matrix()
				if (is.finite(se)) return(se)
				private$cached_values$s_beta_hat_T
			},

			get_degrees_of_freedom = function(){
				private$cached_values$df %||% NA_real_
			},

			supports_likelihood_tests = function(){
				isTRUE(private$combine_reservoir_into_glmm())
			},

			get_likelihood_test_spec = function(){
				private$shared(estimate_only = FALSE)
				ctx = private$cached_values$likelihood_test_context
				if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
				j_treat = as.integer(ctx$j_treat)
				list(
					j = j_treat,
					full_fit = private$cached_mod,
					fit_null = function(delta){
						fast_clogit_plus_glmm_cpp(
							X_disc = ctx$X_disc,
							y_disc = ctx$y_disc,
							X_conc = ctx$X_conc,
							y_conc = ctx$y_conc,
							group_conc = ctx$group_conc,
							start = as.numeric(ctx$start),
							has_discordant = ctx$has_discordant,
							has_concordant = ctx$has_concordant,
							estimate_only = FALSE,
							max_abs_log_sigma = private$max_abs_log_sigma,
							fixed_idx = j_treat,
							fixed_values = delta,
							optimization_alg = private$optimization_alg
						)
					},
					score = function(fit){
						get_clogit_plus_glmm_score_cpp(
							ctx$X_disc, ctx$y_disc, ctx$X_conc, ctx$y_conc, ctx$group_conc,
							as.numeric(fit$params %||% fit$b),
							ctx$has_discordant, ctx$has_concordant,
							private$max_abs_log_sigma
						)
					},
					observed_information = function(fit){
						-get_clogit_plus_glmm_hessian_cpp(
							ctx$X_disc, ctx$y_disc, ctx$X_conc, ctx$y_conc, ctx$group_conc,
							as.numeric(fit$params %||% fit$b),
							ctx$has_discordant, ctx$has_concordant,
							private$max_abs_log_sigma
						)
					},
					information = function(fit){
						-get_clogit_plus_glmm_hessian_cpp(
							ctx$X_disc, ctx$y_disc, ctx$X_conc, ctx$y_conc, ctx$group_conc,
							as.numeric(fit$params %||% fit$b),
							ctx$has_discordant, ctx$has_concordant,
							private$max_abs_log_sigma
						)
					},
					neg_loglik = function(fit){
						as.numeric(fit$neg_loglik %||% fit$neg_ll)
					}
				)
			},

			shared = function(estimate_only = FALSE){
				if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
				if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
				private$clear_nonestimable_state()
				private$cached_values$likelihood_test_context = NULL
				private$cached_mod = NULL

			data_parts = private$build_clogit_plus_glmm_data()
			if (is.null(data_parts)){
				private$cache_nonestimable_estimate("clogit_plus_glmm_no_data")
				return(invisible(NULL))
			}

			fit = if (private$combine_reservoir_into_glmm()){
				private$fit_clogit_plus_glmm(data_parts, estimate_only = estimate_only)
			} else {
				private$fit_clogit_plus_glmm_ivwc(data_parts, estimate_only = estimate_only)
			}

				if (is.null(fit)){
				private$cache_nonestimable_estimate("clogit_plus_glmm_fit_failed")
				return(invisible(NULL))
				}

				if (private$combine_reservoir_into_glmm() && !is.null(fit$mod)) {
					j_treat = if (isTRUE(data_parts$has_concordant)) 2L else 1L
					private$cached_mod = fit$mod
					private$cached_values$likelihood_test_context = list(
						X_disc = as.matrix(data_parts$X_disc),
						y_disc = as.numeric(data_parts$y_disc),
						X_conc = as.matrix(data_parts$X_conc),
						y_conc = as.numeric(data_parts$y_conc),
						group_conc = as.integer(data_parts$group_conc),
						start = as.numeric(fit$mod$b),
						has_discordant = data_parts$has_discordant,
						has_concordant = data_parts$has_concordant,
						j_treat = j_treat
					)
				}
				private$cached_values$beta_hat_T = as.numeric(fit$beta_T)
			if (!estimate_only) private$cached_values$s_beta_hat_T = as.numeric(fit$se_beta_T)
			private$cached_values$df = NA_real_
			invisible(NULL)
		},

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Re-read design variables which might have been transformed during randomization
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			
			# Recompute basic match data for the new w/y
			private$compute_basic_match_data()
			
			# Use the same joint-likelihood logic for the point estimate
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		build_clogit_plus_glmm_data = function(include_reservoir = private$combine_reservoir_into_glmm()){
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			i_matched = which(m_vec > 0L)

			p = ncol(as.matrix(private$X))
			X_m = if (p > 0L) as.matrix(private$get_X()[i_matched, drop = FALSE]) else matrix(nrow = length(i_matched), ncol = 0L)
			y_m = as.numeric(private$y[i_matched])
			w_m = as.numeric(private$w[i_matched])
			strata_m = as.integer(m_vec[i_matched])

			discordant_rows = integer(0)
			concordant_rows = integer(0)
			if (length(i_matched) > 0L){
				pair_rows = split(seq_along(i_matched), strata_m)
				for (rows in pair_rows){
					if (length(unique(y_m[rows])) > 1L){
						discordant_rows = c(discordant_rows, rows)
					} else {
						concordant_rows = c(concordant_rows, rows)
					}
				}
			}

			has_discordant = length(discordant_rows) > 0L
			has_concordant = length(concordant_rows) > 0L

			# Discordant: Clogit (diff-matrix)
			X_disc = matrix(0, 0, 0)
			y_disc = numeric(0)
			if (has_discordant){
				disc = collect_discordant_pairs_cpp(
					as.double(y_m[discordant_rows]),
					as.double(w_m[discordant_rows]),
					X_m[discordant_rows, drop = FALSE],
					as.integer(strata_m[discordant_rows])
				)
				if (disc$nd > 0L){
					X_disc = if (p > 0L) cbind(disc$t_diffs, disc$X_diffs) else matrix(disc$t_diffs, ncol = 1L)
					y_disc = disc$y_01
				} else {
					has_discordant = FALSE
				}
			}

			# Concordant (plus optionally Reservoir): GLMM
			X_conc = matrix(0, 0, 0)
			y_conc = numeric(0)
			group_conc = integer(0)
			
			curr_y = numeric(0)
			curr_w = numeric(0)
			curr_X = matrix(nrow = 0, ncol = p)
			curr_g = integer(0)

			if (has_concordant){
				curr_y = y_m[concordant_rows]
				curr_w = w_m[concordant_rows]
				if (p > 0L) curr_X = X_m[concordant_rows, , drop = FALSE]
				curr_g = strata_m[concordant_rows]
			}

			if (include_reservoir && KKstats$nRT > 0L && KKstats$nRC > 0L){
				curr_y = c(curr_y, as.numeric(KKstats$y_reservoir))
				curr_w = c(curr_w, as.numeric(KKstats$w_reservoir))
				if (p > 0L) curr_X = rbind(curr_X, as.matrix(KKstats$X_reservoir))
				max_g = if (length(curr_g) > 0) max(curr_g) else 0L
				curr_g = c(curr_g, max_g + seq_len(KKstats$nRT + KKstats$nRC))
			}

			if (length(curr_y) > 0L){
				# Ensure consistent n_rows for Intercept and treatment
				n_conc = length(curr_y)
				X_conc = if (p > 0L) {
					cbind(Intercept = rep(1, n_conc), treatment = curr_w, curr_X)
				} else {
					cbind(Intercept = rep(1, n_conc), treatment = curr_w)
				}
				y_conc = curr_y
				group_conc = curr_g
			}

			if (!has_discordant && length(y_conc) == 0L) return(NULL)

			list(
				has_discordant = has_discordant,
				X_disc         = X_disc,
				y_disc         = y_disc,
				has_concordant = length(y_conc) > 0L,
				X_conc         = X_conc,
				y_conc         = y_conc,
				group_conc     = group_conc
			)
		},

		fit_clogit_plus_glmm = function(data_parts, estimate_only = FALSE){
			start = private$get_clogit_plus_glmm_start(data_parts)
			
			res = tryCatch(
				fast_clogit_plus_glmm_cpp(
					X_disc          = as.matrix(data_parts$X_disc),
					y_disc          = as.numeric(data_parts$y_disc),
					X_conc          = as.matrix(data_parts$X_conc),
					y_conc          = as.numeric(data_parts$y_conc),
					group_conc      = as.integer(data_parts$group_conc),
					start           = as.numeric(start),
					has_discordant  = data_parts$has_discordant,
					has_concordant  = data_parts$has_concordant,
					estimate_only   = estimate_only,
					max_abs_log_sigma = private$max_abs_log_sigma,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)
				if (is.null(res) || !is.finite(res$beta_T)) return(NULL)
				if (!estimate_only && (!is.finite(res$se_beta_T) || res$se_beta_T <= 0)) return(NULL)
				
				res$mod = res
				res
			},

		fit_clogit_plus_glmm_ivwc = function(data_parts, estimate_only = FALSE){
			# Fit discordant part (Clogit)
			beta_m = NA_real_
			se_m = NA_real_
			if (data_parts$has_discordant){
				res_m = tryCatch(
					fast_logistic_regression_with_var(data_parts$X_disc, data_parts$y_disc, j = 1L),
					error = function(e) NULL
				)
				if (!is.null(res_m) && is.finite(res_m$b[1])){
					beta_m = res_m$b[1]
					se_m = sqrt(res_m$ssq_b_j)
				}
			}

			# Fit concordant/reservoir part (GLMM)
			beta_r = NA_real_
			se_r = NA_real_
			if (data_parts$has_concordant){
				# Parameters: [Intercept, beta_T, beta_xs, log_sigma]
				# We want beta_T (index 2)
				start = private$get_clogit_plus_glmm_start(data_parts)
				res_r = tryCatch(
					fast_logistic_glmm_cpp(
						X = as.matrix(data_parts$X_conc),
						y = as.numeric(data_parts$y_conc),
						group_id = as.integer(data_parts$group_conc),
						j_T = 2L,
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					),
					error = function(e) NULL
				)
				if (!is.null(res_r) && is.finite(res_r$beta_T)){
					beta_r = res_r$beta_T
					se_r = res_r$se_beta_T
				}
			}

			m_ok = is.finite(beta_m) && (estimate_only || (is.finite(se_m) && se_m > 0))
			r_ok = is.finite(beta_r) && (estimate_only || (is.finite(se_r) && se_r > 0))

			if (m_ok && r_ok){
				pooled = private$weighted_average(beta_m, se_m %||% 1, beta_r, se_r %||% 1)
				list(beta_T = pooled$beta, se_beta_T = if (estimate_only) NA_real_ else pooled$se)
			} else if (m_ok){
				list(beta_T = beta_m, se_beta_T = se_m)
			} else if (r_ok){
				list(beta_T = beta_r, se_beta_T = se_r)
			} else {
				NULL
			}
		},

		get_clogit_plus_glmm_start = function(data_parts){
			if (data_parts$has_concordant){
				X = data_parts$X_conc
				y = data_parts$y_conc
				fit = tryCatch(fast_logistic_regression_cpp(X, y), error = function(e) NULL)
				if (!is.null(fit) && all(is.finite(fit$b))) return(c(fit$b, 0.0))
			}
			p = if (is.null(data_parts$X_disc)) 0L else max(ncol(data_parts$X_disc) - 1L, 0L)
			return(rep(0, p + 3)) # [Intercept, beta_T, beta_xs, log_sigma]
		},

		# Overridden by concrete classes
		combine_reservoir_into_glmm = function() TRUE,

		compute_basic_match_data = function(){
			private$cached_values$KKstats = compute_zhang_match_data_cpp(private$w, private$m, private$y, private$get_X())
		},

		clear_nonestimable_state = function(){
			private$cached_values$beta_hat_T = NULL
			private$cached_values$s_beta_hat_T = NULL
		},

		# --- Utils ---
		weighted_average = function(beta1, se1, beta2, se2){
			w1 = 1 / se1^2
			w2 = 1 / se2^2
			list(
				beta = (w1 * beta1 + w2 * beta2) / (w1 + w2),
				se   = sqrt(1 / (w1 + w2))
			)
		},

		log_sum_exp = function(x){
			m = max(x)
			if (!is.finite(m)) return(m)
			m + log(sum(exp(x - m)))
		},

		log1pexp = function(x){
			ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
		}
	)
)
