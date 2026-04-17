#'  internalAbstract Conditional Logistic Plus GLMM Inference
#'
#' Fits one likelihood with a conditional-logistic contribution from discordant
#' matched pairs and a random-intercept logistic GLMM contribution from concordant
#' matched pairs. The fixed effects are shared by both likelihood components. The
#' asymptotic variance is computed from the inverse of the summed observed Fisher
#' information matrices for the two components.
#'
#' @export
InferenceAbstractKKClogitPlusGLMM = R6::R6Class("InferenceAbstractKKClogitPlusGLMM",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose			Whether to print progress messages.
		#' @param max_abs_reasonable_coef Maximum absolute value allowed for the
		#'   treatment estimate and its standard error before the fit is treated as
		#'   non-estimable. Defaults to \code{1e4}.
		#' @param max_abs_log_sigma Maximum absolute value allowed for the log
		#'   random-intercept standard deviation in the concordant-pair GLMM
		#'   likelihood. Defaults to \code{8}.
		initialize = function(
				des_obj,
				verbose = FALSE,
				max_abs_reasonable_coef = 1e4,
				max_abs_log_sigma = 8
			){
			if (should_run_asserts()) {
				assertNumeric(max_abs_reasonable_coef, lower = .Machine$double.xmin, len = 1, any.missing = FALSE)
				assertNumeric(max_abs_log_sigma, lower = .Machine$double.xmin, len = 1, any.missing = FALSE)
			}
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose)
			private$max_abs_reasonable_coef = max_abs_reasonable_coef
			private$max_abs_log_sigma = max_abs_log_sigma
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_clogit_plus_glmm(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		max_abs_reasonable_coef = 1e4,
		max_abs_log_sigma = 8,

		include_covariates = function() TRUE,
		combine_reservoir_into_glmm = function() stop(class(self)[1], " must implement combine_reservoir_into_glmm()"),

		get_standard_error = function(){
			private$shared_clogit_plus_glmm(estimate_only = FALSE)
			private$cached_values$s_beta_hat_T
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			self$compute_treatment_estimate(estimate_only = TRUE)
		},

		shared_clogit_plus_glmm = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

			data_parts = private$build_clogit_plus_glmm_data()
			if (is.null(data_parts) || !data_parts$has_discordant){
				private$cache_nonestimable_estimate("kk_clogit_plus_glmm_no_discordant_pairs")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			fit = if (private$combine_reservoir_into_glmm()) {
				private$fit_clogit_plus_glmm(data_parts, estimate_only = estimate_only)
			} else {
				private$fit_clogit_plus_glmm_ivwc(data_parts, estimate_only = estimate_only)
			}
			if (is.null(fit)){
				private$cache_nonestimable_estimate("kk_clogit_plus_glmm_fit_failed")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(fit$beta_T)
			if (!is.finite(private$cached_values$beta_hat_T) ||
			    abs(private$cached_values$beta_hat_T) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_clogit_plus_glmm_extreme_estimate")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			if (!estimate_only){
				se = as.numeric(fit$se_beta_T)
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef) se else NA_real_
				if (!is.finite(private$cached_values$s_beta_hat_T)){
					private$cache_nonestimable_se("kk_clogit_plus_glmm_standard_error_unavailable")
					private$cached_values$is_z = TRUE
					return(invisible(NULL))
				}
			}
			private$clear_nonestimable_state()
			private$cached_values$is_z = TRUE
			invisible(NULL)
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

			p = if (private$include_covariates()) ncol(private$get_X()) else 0L
			X_m = if (p > 0L) as.matrix(private$get_X()[i_matched, , drop = FALSE]) else matrix(nrow = length(i_matched), ncol = 0L)
			y_m = as.numeric(private$y[i_matched])
			w_m = as.numeric(private$w[i_matched])
			strata_m = as.integer(m_vec[i_matched])

			pair_rows = split(seq_along(i_matched), strata_m)
			pair_rows = pair_rows[vapply(pair_rows, length, integer(1)) == 2L]
			if (length(pair_rows) == 0L) return(NULL)

			y_sum = vapply(pair_rows, function(ii) sum(y_m[ii]), numeric(1))
			discordant_rows = unlist(pair_rows[y_sum == 1], use.names = FALSE)
			concordant_pairs = pair_rows[y_sum == 0 | y_sum == 2]

			X_disc = NULL
			y_disc = NULL
			if (length(discordant_rows) > 0L){
				disc_strata = strata_m[discordant_rows]
				disc = collect_discordant_pairs_cpp(
					as.double(y_m[discordant_rows]),
					as.double(w_m[discordant_rows]),
					X_m[discordant_rows, , drop = FALSE],
					as.integer(disc_strata)
				)
				if (disc$nd > 0L){
					X_disc = if (p > 0L) cbind(disc$t_diffs, disc$X_diffs) else matrix(disc$t_diffs, ncol = 1L)
					y_disc = as.numeric(disc$y_01)
				}
			}

			X_conc = NULL
			y_conc = NULL
			group_conc = NULL
			if (length(concordant_pairs) > 0L){
				conc_rows = unlist(concordant_pairs, use.names = FALSE)
				X_conc = cbind(Intercept = 1, w = w_m[conc_rows], X_m[conc_rows, , drop = FALSE])
				y_conc = y_m[conc_rows]
				group_conc = rep(seq_along(concordant_pairs), each = 2L)
			}
			if (include_reservoir && KKstats$nRT > 0 && KKstats$nRC > 0){
				X_r = if (p > 0L) as.matrix(KKstats$X_reservoir) else matrix(nrow = length(KKstats$y_reservoir), ncol = 0L)
				X_r = cbind(Intercept = 1, w = KKstats$w_reservoir, X_r)
				next_group = if (is.null(group_conc)) 0L else max(group_conc)
				X_conc = if (is.null(X_conc)) X_r else rbind(X_conc, X_r)
				y_conc = c(y_conc, as.numeric(KKstats$y_reservoir))
				group_conc = c(group_conc, next_group + seq_along(KKstats$y_reservoir))
			}

			beta_names = c("w", if (p > 0L) colnames(private$get_X()) else character())
			if (is.null(beta_names) || any(!nzchar(beta_names))){
				beta_names = c("w", if (p > 0L) paste0("x", seq_len(p)) else character())
			}
			list(
				X_disc = X_disc,
				y_disc = y_disc,
				X_conc = X_conc,
				y_conc = y_conc,
				group_conc = group_conc,
				beta_names = beta_names,
				has_discordant = !is.null(X_disc) && nrow(X_disc) > 0L,
				has_concordant = !is.null(X_conc) && nrow(X_conc) > 0L
			)
		},

		diagnose_clogit_plus_glmm_components = function(){
			matched_data = private$build_clogit_plus_glmm_data(include_reservoir = FALSE)
			combined_data = private$build_clogit_plus_glmm_data(include_reservoir = TRUE)
			matched_fit = if (!is.null(matched_data)) private$fit_clogit_plus_glmm(matched_data, estimate_only = FALSE) else NULL
			combined_fit = if (!is.null(combined_data)) private$fit_clogit_plus_glmm(combined_data, estimate_only = FALSE) else NULL
			res_fit = private$fit_reservoir_logistic(estimate_only = FALSE)

			matched_se = if (!is.null(matched_fit)) matched_fit$se_beta_T else NA_real_
			reservoir_se = if (!is.null(res_fit)) res_fit$se_beta_T else NA_real_
			combined_se = if (!is.null(combined_fit)) combined_fit$se_beta_T else NA_real_
			matched_info = if (is.finite(matched_se) && matched_se > 0) 1 / matched_se^2 else NA_real_
			reservoir_info = if (is.finite(reservoir_se) && reservoir_se > 0) 1 / reservoir_se^2 else NA_real_
			combined_info = if (is.finite(combined_se) && combined_se > 0) 1 / combined_se^2 else NA_real_

			list(
				matched_component_estimate = if (!is.null(matched_fit)) matched_fit$beta_T else NA_real_,
				matched_component_se = matched_se,
				matched_observed_info_beta_T = matched_info,
				reservoir_component_estimate = if (!is.null(res_fit)) res_fit$beta_T else NA_real_,
				reservoir_component_se = reservoir_se,
				reservoir_observed_info_beta_T = reservoir_info,
				combined_component_estimate = if (!is.null(combined_fit)) combined_fit$beta_T else NA_real_,
				combined_component_se = combined_se,
				combined_observed_info_beta_T = combined_info,
				combined_minus_matched_observed_info_beta_T =
					if (is.finite(combined_info) && is.finite(matched_info)) combined_info - matched_info else NA_real_,
				combined_minus_matched_plus_reservoir_observed_info_beta_T =
					if (is.finite(combined_info) && is.finite(matched_info) && is.finite(reservoir_info)) combined_info - matched_info - reservoir_info else NA_real_
			)
		},

		fit_clogit_plus_glmm_ivwc = function(data_parts, estimate_only = FALSE){
			matched_fit = private$fit_clogit_plus_glmm(data_parts, estimate_only = FALSE)
			if (is.null(matched_fit)) return(NULL)
			res_fit = private$fit_reservoir_logistic(estimate_only = FALSE)
			if (is.null(res_fit)){
				if (estimate_only) return(list(beta_T = matched_fit$beta_T, se_beta_T = NA_real_))
				return(matched_fit)
			}
			ssq_m = matched_fit$se_beta_T^2
			ssq_r = res_fit$se_beta_T^2
			if (!is.finite(ssq_m) || !is.finite(ssq_r) || ssq_m <= 0 || ssq_r <= 0) return(matched_fit)
			w_star = ssq_r / (ssq_r + ssq_m)
			beta = w_star * matched_fit$beta_T + (1 - w_star) * res_fit$beta_T
			se = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			list(beta_T = beta, se_beta_T = se)
		},

		fit_reservoir_logistic = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats) || KKstats$nRT <= 0 || KKstats$nRC <= 0) return(NULL)
			X_r = as.matrix(KKstats$X_reservoir)
			X_full = cbind(Intercept = 1, w = KKstats$w_reservoir, X_r)
			mod = tryCatch(
				if (estimate_only) fast_logistic_regression(X_full, KKstats$y_reservoir) else fast_logistic_regression_with_var(X_full, KKstats$y_reservoir, j = 2L),
				error = function(e) NULL
			)
			if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2L])) return(NULL)
			if (estimate_only) return(list(beta_T = as.numeric(mod$b[2L]), se_beta_T = NA_real_))
			ssq = as.numeric(mod$ssq_b_j)
			if (!is.finite(ssq) || ssq <= 0) return(NULL)
			list(beta_T = as.numeric(mod$b[2L]), se_beta_T = sqrt(ssq))
		},

		fit_clogit_plus_glmm = function(data_parts, estimate_only = FALSE){
			q = length(data_parts$beta_names)
			has_concordant = data_parts$has_concordant
			has_discordant = data_parts$has_discordant

			start = private$get_clogit_plus_glmm_start(data_parts)
			if (is.null(start)) return(NULL)
			X_disc = if (has_discordant) data_parts$X_disc else matrix(numeric(), nrow = 0L, ncol = q)
			y_disc = if (has_discordant) data_parts$y_disc else numeric()
			X_conc = if (has_concordant) data_parts$X_conc else matrix(numeric(), nrow = 0L, ncol = q + 1L)
			y_conc = if (has_concordant) data_parts$y_conc else numeric()
			group_conc = if (has_concordant) as.integer(data_parts$group_conc) else integer()

			starts = private$get_clogit_plus_glmm_starts(start, q, has_concordant)
			best = NULL
			for (start_par in starts){
				fit = tryCatch(
					fast_clogit_plus_glmm_cpp(
						X_disc = X_disc,
						y_disc = y_disc,
						X_conc = X_conc,
						y_conc = y_conc,
						group_conc = group_conc,
						start = start_par,
						has_discordant = has_discordant,
						has_concordant = has_concordant,
						estimate_only = estimate_only,
						max_abs_log_sigma = private$max_abs_log_sigma
					),
					error = function(e) NULL
				)
				if (is.null(fit) || !isTRUE(fit$converged) || !is.finite(fit$neg_loglik)) next
				if (is.null(best) || fit$neg_loglik < best$neg_loglik) best = fit
			}
			if (is.null(best)) return(NULL)

			beta_T = if (has_concordant) best$b[2L] else best$b[1L]
			if (estimate_only) return(list(beta_T = beta_T, se_beta_T = NA_real_))
			ssq = as.numeric(best$ssq_b_j)
			if (!is.finite(ssq) || ssq <= 0) return(NULL)
			list(beta_T = beta_T, se_beta_T = sqrt(ssq))
		},

		get_clogit_plus_glmm_start = function(data_parts){
			if (data_parts$has_concordant){
				X = data_parts$X_conc
				y = data_parts$y_conc
				mod = tryCatch(fast_logistic_regression(X, y), error = function(e) NULL)
				beta = if (!is.null(mod) && length(mod$b) == ncol(X) && all(is.finite(mod$b))){
					as.numeric(mod$b)
				} else {
					rep(0, ncol(X))
				}
				return(c(beta, log_sigma = 0))
			}
			if (data_parts$has_discordant){
				X = data_parts$X_disc
				y = data_parts$y_disc
				mod = tryCatch(fast_logistic_regression(X, y), error = function(e) NULL)
				beta = if (!is.null(mod) && length(mod$b) == ncol(X) && all(is.finite(mod$b))) {
					as.numeric(mod$b)
				} else {
					rep(0, ncol(X))
				}
				return(beta)
			}
			NULL
		},

		get_clogit_plus_glmm_starts = function(start, q, has_concordant){
			if (!has_concordant) return(list(start))
			list(
				start,
				replace(start, q + 2L, log(0.5)),
				replace(start, q + 2L, log(2))
			)
		},

		neg_clogit_loglik = function(beta_no_intercept, data_parts){
			eta = as.vector(data_parts$X_disc %*% beta_no_intercept)
			ll = data_parts$y_disc * eta - private$log1pexp(eta)
			if (any(!is.finite(ll))) return(1e100)
			-sum(ll)
		},

		neg_concordant_glmm_loglik = function(par, data_parts){
			q = length(data_parts$beta_names)
			beta = par[seq_len(q + 1L)]
			sigma = exp(par[q + 2L])
			quad = private$gauss_hermite_rule(20L)
			log_norm_weights = log(quad$weights) - 0.5 * log(pi)
			b_vals = sqrt(2) * sigma * quad$nodes
			total = 0
			for (g in unique(data_parts$group_conc)){
				ii = which(data_parts$group_conc == g)
				eta0 = as.vector(data_parts$X_conc[ii, , drop = FALSE] %*% beta)
				log_terms = vapply(seq_along(b_vals), function(k){
					eta = eta0 + b_vals[k]
					sum(data_parts$y_conc[ii] * eta - private$log1pexp(eta)) + log_norm_weights[k]
				}, numeric(1))
				ll_g = private$log_sum_exp(log_terms)
				if (!is.finite(ll_g)) return(1e100)
				total = total + ll_g
			}
			-total
		},

		clogit_plus_glmm_observed_info = function(par, data_parts){
			q = length(data_parts$beta_names)
			has_concordant = data_parts$has_concordant
			p_full = if (has_concordant) q + 2L else q
			info = matrix(0, nrow = p_full, ncol = p_full)

			if (data_parts$has_discordant){
				neg_cond = function(beta_no_intercept) private$neg_clogit_loglik(beta_no_intercept, data_parts)
				beta_no_intercept = if (has_concordant) par[2L:(q + 1L)] else par[seq_len(q)]
				h_cond = tryCatch(numDeriv::hessian(neg_cond, beta_no_intercept), error = function(e) NULL)
				if (is.null(h_cond) || any(!is.finite(h_cond))) return(NULL)
				if (has_concordant) {
					idx = 2L:(q + 1L)
					info[idx, idx] = info[idx, idx] + h_cond
				} else {
					info = info + h_cond
				}
			}

			if (has_concordant){
				neg_glmm = function(par_full) private$neg_concordant_glmm_loglik(par_full, data_parts)
				h_glmm = tryCatch(numDeriv::hessian(neg_glmm, par), error = function(e) NULL)
				if (is.null(h_glmm) || any(!is.finite(h_glmm))) return(NULL)
				info = info + h_glmm
			}
			info = (info + t(info)) / 2
			if (any(!is.finite(diag(info)))) return(NULL)
			info
		},

		gauss_hermite_rule = function(n){
			cache_name = paste0("gh_rule_", n)
			if (!is.null(private$cached_values[[cache_name]])) return(private$cached_values[[cache_name]])
			i = seq_len(n - 1L)
			J = matrix(0, n, n)
			J[cbind(i, i + 1L)] = sqrt(i / 2)
			J[cbind(i + 1L, i)] = sqrt(i / 2)
			eig = eigen(J, symmetric = TRUE)
			nodes = eig$values
			weights = sqrt(pi) * eig$vectors[1L, ]^2
			o = order(nodes)
			rule = list(nodes = nodes[o], weights = weights[o])
			private$cached_values[[cache_name]] = rule
			rule
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
