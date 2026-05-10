#' Abstract class for Weibull Frailty / Standard Weibull Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Weibull AFT GLMM for matched pairs.
#' The matched-pair component uses a shared log-normal random intercept per pair,
#' fitted by the package's native Rcpp likelihood optimizer. The reservoir component uses standard Weibull AFT regression.
#' The two treatment-effect estimates are combined by inverse-variance weighting.
#'
#' @details
#' This compound estimator accounts for the dependence within matched pairs by
#' modeling it as a shared frailty.
#'
#' \strong{Univariate} (\code{ncol(as.matrix(private$X)) == 0}): uses the native
#' Rcpp Weibull frailty likelihood with \code{formula = survival::Surv(y, dead) ~ w}
#' and a pair-level random intercept.
#'
#' \strong{Multivariate} (\code{ncol(as.matrix(private$X)) > 0}): fits the same
#' native Rcpp likelihood with covariate adjustment, dropping rank-deficient
#' columns when needed.
#'
#' @keywords internal
InferenceAbstractKKWeibullFrailtyIVWC = R6::R6Class("InferenceAbstractKKWeibullFrailtyIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		},

		#' @description
		#' Returns the estimated treatment effect (log-time ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
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

		#' @description
		#' Computes the asymptotic p-value.
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

		#' @description
		#' Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			inf_obj = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			inf_obj
		}
	),

	private = list(

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design and parameters from the original data
			if (is.null(private$best_X_colnames_matched) && is.null(private$best_X_colnames_reservoir)){
				private$shared()
			}

			# If we still don't have enough (e.g., initial fit failed), fall back to standard
			if (is.null(private$best_X_colnames_matched) && is.null(private$best_X_colnames_reservoir)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			X_data = private$get_X()

			# Matched pairs component
			beta_m = NA_real_
			ssq_m = NA_real_
			if (m > 0 && !is.null(private$best_X_colnames_matched)){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_matched = which(m_vec > 0L)

				X_cov = X_data[i_matched, intersect(private$best_X_colnames_matched, colnames(X_data)), drop = FALSE]
				X = cbind(w = private$w[i_matched], X_cov)

				fit_m = .fit_weibull_frailty(
					y = private$y[i_matched],
					dead = private$dead[i_matched],
					X = X,
					pair_id = m_vec[i_matched],
					estimate_only = estimate_only,
					optimization_alg = private$optimization_alg
				)
				if (!is.null(fit_m) && is.finite(fit_m$beta)){
					beta_m = fit_m$beta
					if (!estimate_only && !is.null(fit_m$ssq) && is.finite(fit_m$ssq) && fit_m$ssq > 0){
						ssq_m = fit_m$ssq
					}
				}
			}

			# Reservoir component
			beta_r = NA_real_
			ssq_r = NA_real_
			if (nRT > 0 && nRC > 0 && !is.null(private$best_X_colnames_reservoir)){
				m_vec = private$m
				if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
				m_vec[is.na(m_vec)] = 0L
				i_reservoir = which(m_vec == 0L)

				X_cov = X_data[i_reservoir, intersect(private$best_X_colnames_reservoir, colnames(X_data)), drop = FALSE]
				X = cbind(w = private$w[i_reservoir], X_cov)

				fit_r = .fit_standard_weibull_aft_from_matrix(
					y = private$y[i_reservoir],
					dead = private$dead[i_reservoir],
					X = X,
					estimate_only = estimate_only
				)
				if (!is.null(fit_r) && is.finite(fit_r$beta)){
					beta_r = fit_r$beta
					if (!estimate_only && !is.null(fit_r$ssq) && is.finite(fit_r$ssq) && fit_r$ssq > 0){
						ssq_r = fit_r$ssq
					}
				}
			}

			# Pooling
			m_ok = is.finite(beta_m) && (estimate_only || is.finite(ssq_m))
			r_ok = is.finite(beta_r) && (estimate_only || is.finite(ssq_r))

			if (m_ok && r_ok){
				if (estimate_only) {
					ssq_m_orig = private$cached_values$ssq_beta_T_matched
					ssq_r_orig = private$cached_values$ssq_beta_T_reservoir
					if (!is.null(ssq_m_orig) && !is.null(ssq_r_orig) && is.finite(ssq_m_orig) && is.finite(ssq_r_orig)){
						w_star = ssq_r_orig / (ssq_r_orig + ssq_m_orig)
						return(w_star * beta_m + (1 - w_star) * beta_r)
					}
					return(0.5 * beta_m + 0.5 * beta_r)
				}
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		},

		best_X_colnames_matched = NULL,
		best_X_colnames_reservoir = NULL,

		frailty_for_matched_pairs = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0L)
			if (length(i_matched) == 0L) return(invisible(NULL))

			X_full = if (ncol(as.matrix(private$X)) == 0L){
				matrix(private$w[i_matched], ncol = 1L, dimnames = list(NULL, "w"))
			} else {
				cbind(w = private$w[i_matched], private$get_X()[i_matched, , drop = FALSE])
			}

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				fit_fun = function(X_fit, keep){
					res = .fit_weibull_frailty(
						y = private$y[i_matched],
						dead = private$dead[i_matched],
						X = X_fit,
						pair_id = m_vec[i_matched],
						estimate_only = estimate_only,
						optimization_alg = private$optimization_alg
					)
					res
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || !is.finite(mod$beta)) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq) && mod$ssq > 0
				}
			)
			
			if (!is.null(attempt$fit)){
				private$cached_values$beta_T_matched = attempt$fit$beta
				private$cached_values$ssq_beta_T_matched = attempt$fit$ssq
				private$best_X_colnames_matched = setdiff(colnames(attempt$X_fit), "w")
			}
		},

		weibull_for_reservoir = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_reservoir = which(m_vec == 0L)
			if (length(i_reservoir) == 0L) return(invisible(NULL))

			X_full = if (ncol(as.matrix(private$X)) == 0L){
				matrix(private$w[i_reservoir], ncol = 1L, dimnames = list(NULL, "w"))
			} else {
				cbind(w = private$w[i_reservoir], private$get_X()[i_reservoir, , drop = FALSE])
			}

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				fit_fun = function(X_fit, keep){
					res = .fit_standard_weibull_aft_from_matrix(
						y = private$y[i_reservoir],
						dead = private$dead[i_reservoir],
						X = X_fit,
						estimate_only = estimate_only
					)
					res
				},
				fit_ok = function(mod, X_fit, keep){
					if (is.null(mod) || !is.finite(mod$beta)) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(mod$ssq) && mod$ssq > 0
				}
			)
			
			if (!is.null(attempt$fit)){
				private$cached_values$beta_T_reservoir = attempt$fit$beta
				private$cached_values$ssq_beta_T_reservoir = attempt$fit$ssq
				private$best_X_colnames_reservoir = setdiff(colnames(attempt$X_fit), "w")
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			if (m > 0){
				private$frailty_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
			       (!estimate_only && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0 || estimate_only)

			if (nRT > 0 && nRC > 0){
				private$weibull_for_reservoir(estimate_only = estimate_only)
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) &&
			       (!estimate_only && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0 || estimate_only)

			if (m_ok && r_ok){
				if (estimate_only){
					if (!is.null(ssq_m) && !is.null(ssq_r) && is.finite(ssq_m) && is.finite(ssq_r)){
						w_star = ssq_r / (ssq_r + ssq_m)
						private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
					} else {
						private$cached_values$beta_hat_T = 0.5 * beta_m + 0.5 * beta_r
					}
					return(invisible(NULL))
				}
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
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
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)
)

#' Abstract class for Weibull Frailty Combined-Likelihood Inference
#'
#' Fits a single joint Weibull frailty model over all KK design data for survival
#' responses. Matched subjects share their pair ID as a frailty cluster; reservoir
#' subjects are assigned to singleton clusters. The treatment effect beta_T and
#' covariate slopes beta_xs are shared across all subjects.
#'
#' @keywords internal
InferenceAbstractKKWeibullFrailtyOneLik = R6::R6Class("InferenceAbstractKKWeibullFrailtyOneLik",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param use_rcpp Whether to use the custom Rcpp likelihood optimizer.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg The optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			private$use_rcpp = use_rcpp
			self$set_optimization_alg(optimization_alg, allow_irls = FALSE)
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		},

		#' @description
		#' Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only Whether to skip standard-error calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an asymptotic confidence interval for the treatment effect.
		#' @param alpha Significance level.
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

		#' @description
		#' Returns a 2-sided p-value for H0: beta_T = delta.
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
		}
	),

	private = list(

			use_rcpp = TRUE,

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
				isTRUE(private$use_rcpp)
			},

			get_likelihood_test_spec = function(){
				private$shared_combined_likelihood(estimate_only = FALSE)
				ctx = private$cached_values$likelihood_test_context
				if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
				X_fit = ctx$X
				y = as.numeric(ctx$y)
				dead = as.numeric(ctx$dead)
				group_id = as.integer(ctx$group_id)
				n_gh = as.integer(ctx$n_gh %||% 20L)
				max_abs_log_sigma = as.numeric(ctx$max_abs_log_sigma %||% 8.0)
				list(
					j = 1L,
					full_fit = private$cached_mod,
					fit_null = function(delta, start = NULL){
						fast_weibull_frailty_cpp(
							y = y,
							dead = dead,
							X = X_fit,
							group_id = group_id,
							start = start %||% private$get_fit_warm_start_for_length("params", length(ctx$start)) %||% as.numeric(ctx$start),
							estimate_only = FALSE,
							n_gh = n_gh,
							max_abs_log_sigma = max_abs_log_sigma,
							fixed_idx = 1L,
							fixed_values = delta,
							optimization_alg = private$optimization_alg
						)
					},
					extract_start = function(fit){
						as.numeric(fit$params %||% c(as.numeric(fit$b), as.numeric(fit$log_sigma_eps), as.numeric(fit$log_sigma_u)))
					},
					score = function(fit){
						as.numeric(fit$score %||% get_weibull_frailty_score_cpp(X_fit, y, dead, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
						},
						observed_information = function(fit){
						as.matrix(fit$observed_information %||% fit$information %||% -get_weibull_frailty_hessian_cpp(X_fit, y, dead, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
						},
						fisher_information = function(fit){
						as.matrix(fit$information %||% fit$observed_information %||% -get_weibull_frailty_hessian_cpp(X_fit, y, dead, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
						},
						information = function(fit){
						as.matrix(fit$information %||% fit$observed_information %||% -get_weibull_frailty_hessian_cpp(X_fit, y, dead, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
						},
				)
			},

			build_design_matrix = function(){
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			if (ncol(as.matrix(private$X)) > 0){
				X_full = cbind(X_full, as.matrix(private$get_X()))
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(1L %in% keep)) keep[r_full] = 1L
					keep = sort(keep)
					X_full = X_full[, keep, drop = FALSE]
				}
			}
			X_full
		},

			shared_rcpp = function(estimate_only = FALSE){
				if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
				if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
				private$cached_values$likelihood_test_context = NULL
				private$cached_mod = NULL

				if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			cluster_id = m_vec
			reservoir_idx = which(cluster_id == 0L)
			if (length(reservoir_idx) > 0L){
				cluster_id[reservoir_idx] = max(cluster_id) + seq_along(reservoir_idx)
			}

			if (sum(private$dead) == 0L){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}

			X_full = private$build_design_matrix()

			res = tryCatch(
				fast_weibull_frailty_cpp(
					y               = private$y,
					dead            = as.numeric(private$dead),
					X               = X_full,
					group_id        = as.integer(cluster_id),
					start           = private$get_fit_warm_start("params"),
					estimate_only   = estimate_only,
					n_gh            = 20L,
					max_abs_log_sigma = 8.0,
					maxit           = 300L,
					eps_g           = 1e-6,
					optimization_alg = private$optimization_alg
				),
				error = function(e) NULL
			)

			if (is.null(res) || !isTRUE(res$converged)){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}

				beta = as.numeric(res$b)[1L]
			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			private$cached_mod = res
			private$set_fit_warm_start(as.numeric(res$params %||% c(as.numeric(res$b), as.numeric(res$log_sigma_eps), as.numeric(res$log_sigma_u))), "params")
			private$cached_values$likelihood_test_context = list(
					X = X_full,
					y = as.numeric(private$y),
					dead = as.numeric(private$dead),
					group_id = as.integer(cluster_id),
					start = as.numeric(res$params %||% c(as.numeric(res$b), as.numeric(res$log_sigma_eps), as.numeric(res$log_sigma_u))),
					n_gh = 20L,
					max_abs_log_sigma = 8.0
				)

				if (!estimate_only){
				ssq = as.numeric(res$ssq_b_T)
				se  = if (is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
				private$cached_values$s_beta_hat_T = se
				}
				private$cached_values$df = NA_real_
				invisible(NULL)
			},

		shared_combined_likelihood = function(estimate_only = FALSE){
			if (private$use_rcpp) return(private$shared_rcpp(estimate_only = estimate_only))
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			cluster_id = m_vec
			reservoir_idx = which(cluster_id == 0L)
			if (length(reservoir_idx) > 0L){
				cluster_id[reservoir_idx] = max(cluster_id) + seq_along(reservoir_idx)
			}

			if (sum(private$dead) == 0L){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}

			X_full = private$build_design_matrix()
			dat = data.frame(y = private$y, dead = private$dead, w = X_full[, "w"], cluster = factor(cluster_id))
			formula_str = "survival::Surv(y, dead) ~ w"

			X_covs = X_full[, colnames(X_full) != "w", drop = FALSE]
			if (ncol(X_covs) > 0){
				colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
				dat = cbind(dat, X_covs)
				formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
			}
			formula_str = paste(formula_str, "+ cluster(cluster)")

			mod = tryCatch(
				suppressWarnings(survival::survreg(as.formula(formula_str), data = dat, dist = "weibull")),
				error = function(e) NULL
			)

			if (is.null(mod)){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}

			coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
			if (is.null(coefs) || !("w" %in% names(coefs))){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}

			beta = as.numeric(coefs["w"])
			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			
			if (!estimate_only) {
				vcv = tryCatch(stats::vcov(mod), error = function(e) NULL)
				if (is.null(vcv) || !("w" %in% rownames(vcv))){
					private$cached_values$s_beta_hat_T = NA_real_
				} else {
					se = sqrt(as.numeric(vcv["w", "w"]))
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
			}
			invisible(NULL)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)
)

#' Weibull Frailty IVWC Inference for KK Designs
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a native Rcpp Weibull AFT frailty model for matched
#' pairs and a standard Weibull AFT model for the reservoir. The estimates are
#' combined via IVWC.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKWeibullFrailtyIVWC$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKWeibullFrailtyIVWC = R6::R6Class("InferenceSurvivalKKWeibullFrailtyIVWC",
	inherit = InferenceAbstractKKWeibullFrailtyIVWC,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL} (default), 
		#'   covariates from the design object are included. Use \code{~ 1} for univariate.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg Optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, verbose = FALSE, optimization_alg = NULL){
			self$set_optimization_alg(optimization_alg)
			super$initialize(des_obj, model_formula = model_formula, verbose = verbose)
		}
	)
)

#' Weibull Frailty Combined-Likelihood Inference for KK Designs
#'
#' This class fits a single joint Weibull frailty model over all KK design data
#' (matched pairs + reservoir) for survival responses.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKWeibullFrailtyOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKWeibullFrailtyOneLik = R6::R6Class("InferenceSurvivalKKWeibullFrailtyOneLik",
	inherit = InferenceAbstractKKWeibullFrailtyOneLik,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param model_formula Optional formula for covariate adjustment. If \code{NULL} (default), 
		#'   covariates from the design object are included. Use \code{~ 1} for univariate.
		#' @param use_rcpp Whether to use the custom Rcpp likelihood optimizer.
		#' @param verbose Whether to print progress messages.
		#' @param optimization_alg The optimization algorithm to use. Default is dispatched via policy.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, optimization_alg = NULL){
			self$set_optimization_alg(optimization_alg)
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose)
		}
	)
)
