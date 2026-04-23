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
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
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
					fit_null = function(delta){
						fast_weibull_frailty_cpp(
							y = y,
							dead = dead,
							X = X_fit,
							group_id = group_id,
							start = as.numeric(ctx$start),
							estimate_only = FALSE,
							n_gh = n_gh,
							max_abs_log_sigma = max_abs_log_sigma,
							fixed_idx = 1L,
							fixed_values = delta,
							optimization_alg = private$optimization_alg
						)
					},
					score = function(fit){
						as.numeric(fit$score %||% get_weibull_frailty_score_cpp(y, dead, X_fit, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
					},
					observed_information = function(fit){
						as.matrix(fit$observed_information %||% fit$information %||% -get_weibull_frailty_hessian_cpp(y, dead, X_fit, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
					},
					information = function(fit){
						as.matrix(fit$information %||% fit$observed_information %||% -get_weibull_frailty_hessian_cpp(y, dead, X_fit, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
					},
					neg_loglik = function(fit){
						as.numeric(fit$neg_loglik %||% fit$neg_ll %||% get_weibull_frailty_neg_loglik_cpp(y, dead, X_fit, group_id, as.numeric(fit$params), n_gh, max_abs_log_sigma))
					}
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
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			X_full = private$build_design_matrix()

			res = tryCatch(
				fast_weibull_frailty_cpp(
					y               = private$y,
					dead            = as.numeric(private$dead),
					X               = X_full,
					group_id        = as.integer(cluster_id),
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
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

				beta = as.numeric(res$b)[1L]
				private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
				private$cached_mod = res
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
				private$cached_values$is_z = TRUE
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
				private$cached_values$is_z         = TRUE
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
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
			if (is.null(coefs) || !("w" %in% names(coefs))){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
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
			private$cached_values$is_z         = TRUE
			invisible(NULL)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		}
	)
)
