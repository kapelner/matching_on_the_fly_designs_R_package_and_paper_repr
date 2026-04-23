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
	inherit = InferenceAsymp,
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
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
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

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

			data_parts = private$build_clogit_plus_glmm_data()
			if (is.null(data_parts)){
				private$cache_nonestimable_estimate("clogit_plus_glmm_no_data")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			fit = if (private$combine_reservoir_into_glmm()){
				private$fit_clogit_plus_glmm(data_parts, estimate_only = estimate_only)
			} else {
				private$fit_clogit_plus_glmm_ivwc(data_parts, estimate_only = estimate_only)
			}

			if (is.null(fit)){
				private$cache_nonestimable_estimate("clogit_plus_glmm_fit_failed")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(fit$beta_T)
			if (!estimate_only) private$cached_values$s_beta_hat_T = as.numeric(fit$se_beta_T)
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
			invisible(NULL)
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
			if (length(i_matched) == 0L) return(NULL)

			p = ncol(as.matrix(private$X))
			X_m = if (p > 0L) as.matrix(private$get_X()[i_matched, drop = FALSE]) else matrix(nrow = length(i_matched), ncol = 0L)
			y_m = as.numeric(private$y[i_matched])
			w_m = as.numeric(private$w[i_matched])
			strata_m = as.integer(m_vec[i_matched])

			pair_rows = split(seq_along(i_matched), strata_m)
			discordant_rows = integer(0)
			concordant_rows = integer(0)
			for (rows in pair_rows){
				if (length(unique(y_m[rows])) > 1L){
					discordant_rows = c(discordant_rows, rows)
				} else {
					concordant_rows = c(concordant_rows, rows)
				}
			}

			has_discordant = length(discordant_rows) > 0L
			has_concordant = length(concordant_rows) > 0L

			# Discordant: Clogit (diff-matrix)
			X_disc = NULL
			y_disc = NULL
			if (has_discordant){
				# Clogit expects diff-matrix and 0/1 pair outcomes. 
				# The C++ helper collect_discordant_pairs_cpp already does this.
				disc = collect_discordant_pairs_cpp(
					as.double(y_m[discordant_rows]),
					as.double(w_m[discordant_rows]),
					X_m[discordant_rows, drop = FALSE],
					as.integer(strata_m[discordant_rows])
				)
				if (disc$nd > 0L){
					X_disc = if (p > 0L) cbind(res$t_diffs, res$X_diffs) else matrix(res$t_diffs, ncol = 1L)
					y_disc = res$y_01
				} else {
					has_discordant = FALSE
				}
			}

			# Concordant (plus optionally Reservoir): GLMM
			X_conc = NULL
			y_conc = NULL
			group_conc = NULL
			if (has_concordant || (include_reservoir && KKstats$nRT > 0L && KKstats$nRC > 0L)){
				X_conc = X_m[concordant_rows, , drop = FALSE]
				y_conc = y_m[concordant_rows]
				w_conc = w_m[concordant_rows]
				group_conc = strata_m[concordant_rows]

				if (include_reservoir){
					y_conc = c(y_conc, as.numeric(KKstats$y_reservoir))
					w_conc = c(w_conc, as.numeric(KKstats$w_reservoir))
					if (p > 0L) X_conc = rbind(X_conc, as.matrix(KKstats$X_reservoir))
					# Reservoir subjects are independent → unique group IDs
					if (KKstats$nRT + KKstats$nRC > 0L){
						max_g = if (length(group_conc) > 0) max(group_conc) else 0L
						group_conc = c(group_conc, max_g + seq_len(KKstats$nRT + KKstats$nRC))
					}
				}
				X_conc = cbind(1, w_conc, X_conc)
			}

			if (!has_discordant && is.null(X_conc)) return(NULL)

			list(
				has_discordant = has_discordant,
				X_disc         = X_disc,
				y_disc         = y_disc,
				has_concordant = !is.null(X_conc),
				X_conc         = X_conc,
				y_conc         = y_conc,
				group_conc     = group_conc
			)
		},

		fit_clogit_plus_glmm = function(data_parts, estimate_only = FALSE){
			# Call the joint likelihood optimizer in C++
			# Parameters: [beta_T, beta_xs, log_sigma] (discordant part has no intercept)
			# BUT GLMM part has intercept.  This specific Rcpp backend handles the logic.
			start = private$get_clogit_plus_glmm_start(data_parts)
			
			res = tryCatch(
				fast_clogit_plus_glmm_cpp(
					X_disc          = as.matrix(data_parts$X_disc %||% matrix(0, 0, 0)),
					y_disc          = as.numeric(data_parts$y_disc %||% numeric(0)),
					X_conc          = as.matrix(data_parts$X_conc %||% matrix(0, 0, 0)),
					y_conc          = as.numeric(data_parts$y_conc %||% numeric(0)),
					group_conc      = as.integer(data_parts$group_conc %||% integer(0)),
					start           = as.numeric(start),
					has_discordant  = data_parts$has_discordant,
					has_concordant  = data_parts$has_concordant,
					estimate_only   = estimate_only,
					max_abs_log_sigma = private$max_abs_log_sigma,
					optimization_alg = private$optimization_alg
				),
				error = function(e) {
					print(paste("DEBUG: fast_clogit_plus_glmm_cpp error:", e$message))
					NULL
				}
			)
			if (is.null(res) || !is.finite(res$beta_T)) return(NULL)
			if (!estimate_only && (!is.finite(res$se_beta_T) || res$se_beta_T <= 0)) return(NULL)
			
			res
		},

		get_clogit_plus_glmm_start = function(data_parts){
			if (data_parts$has_concordant){
				X = data_parts$X_conc
				y = data_parts$y_conc
				# Logistic regression to get starting values
				# Param layout: [Intercept, beta_T, beta_xs]
				fit = fast_logistic_regression_cpp(X, y)
				# clogit_plus_glmm_cpp expects [Intercept, beta_T, beta_xs, log_sigma]
				return(c(fit$b, 0.0))
			}
			# fallback
			p = if (is.null(data_parts$X_disc)) 0L else ncol(data_parts$X_disc)
			return(rep(0, p + 2))
		},

		# Overridden by concrete classes
		combine_reservoir_into_glmm = function() TRUE,

		compute_basic_match_data = function(){
			# Standard Zhang-style summary
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
