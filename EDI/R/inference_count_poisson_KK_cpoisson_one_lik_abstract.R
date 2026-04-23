#' Abstract class for Conditional Poisson Combined-Likelihood Compound Inference
#'
#' Fits a single joint likelihood over all KK design data for count responses.
#' The matched-pair component uses the conditional Poisson likelihood
#' (which simplifies to a binomial contribution per pair), and the reservoir
#' uses an ordinary Poisson likelihood.
#'
#' @keywords internal
InferenceAbstractKKPoissonCPoissonOneLik = R6::R6Class("InferenceAbstractKKPoissonCPoissonOneLik",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		#' @return The treatment estimate.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an asymptotic confidence interval.
		#' @param alpha Confidence level.
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
		#' Computes an asymptotic p-value.
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
			# Use the same joint-likelihood logic for the point estimate
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		compute_treatment_estimate_during_power_simulation = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = 1e-4){
			private$compute_treatment_estimate_during_power_simulation_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

			assert_finite_se = function(){
				if (!is.finite(private$cached_values$s_beta_hat_T)){
					return(invisible(NULL))
				}
			},

			get_standard_error = function(){
				private$shared_combined_likelihood(estimate_only = FALSE)
				private$cached_values$s_beta_hat_T
			},

			get_degrees_of_freedom = function(){
				private$cached_values$df %||% NA_real_
			},

			supports_likelihood_tests = function(){
				!is.null(private$cached_values$likelihood_test_context) || TRUE
			},

			get_likelihood_test_spec = function(){
				private$shared_combined_likelihood(estimate_only = FALSE)
				ctx = private$cached_values$likelihood_test_context
				if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
				yT_v = as.numeric(ctx$yT_v)
				n_k_v = as.numeric(ctx$n_k_v)
				X_diff_v = as.matrix(ctx$X_diff_v)
				y_r_v = as.numeric(ctx$y_r_v)
				w_r_v = as.numeric(ctx$w_r_v)
				X_r_v = as.matrix(ctx$X_r_v)
				list(
					j = 2L,
					full_fit = private$cached_mod,
					fit_null = function(delta){
						fast_cpoisson_combined_with_var_cpp(
							yT_v = yT_v,
							n_k_v = n_k_v,
							X_diff_v = X_diff_v,
							y_r = y_r_v,
							w_r = w_r_v,
							X_r = X_r_v,
							fixed_idx = 2L,
							fixed_values = delta
						)
					},
					score = function(fit){
						as.numeric(fit$score %||% get_cpoisson_combined_score_cpp(yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v, as.numeric(fit$params %||% fit$b)))
					},
					information = function(fit){
						as.matrix(fit$information %||% -get_cpoisson_combined_hessian_cpp(yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v, as.numeric(fit$params %||% fit$b)))
					},
					neg_loglik = function(fit){
						as.numeric(fit$neg_loglik %||% fit$neg_ll)
					}
				)
			},

		set_failed_combined_cache = function(){
			private$cached_values$beta_hat_T   = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z         = TRUE
			private$cache_nonestimable_estimate("kk_cpoisson_combined_fit_failed")
		},

		reduce_combined_covariates = function(X_diff, X_r, w_r){
			# Use a heuristic to ensure X_diff and X_r are jointly full-rank when combined
			# with treatment and reservoir intercept.
			n_p = nrow(X_diff)
			n_r = nrow(X_r)
			p   = ncol(X_diff)
			
			# Layout: [beta_0, beta_T, beta_xs]
			# Pairs:     [0, 1, X_diff]
			# Reservoir: [1, w_r, X_r]
			
			# Handle empty components to avoid cbind recycling warnings
			pairs_part = if (n_p > 0) {
				cbind(Intercept = rep(0, n_p), treatment = rep(1, n_p), X_diff)
			} else {
				matrix(0, nrow = 0, ncol = p + 2)
			}
			
			reservoir_part = if (n_r > 0) {
				cbind(Intercept = rep(1, n_r), treatment = w_r, X_r)
			} else {
				matrix(0, nrow = 0, ncol = p + 2)
			}
			
			X_stack = rbind(pairs_part, reservoir_part)
			
			if (nrow(X_stack) == 0) return(integer(0))
			
			qr_res = qr(X_stack)
			if (qr_res$rank < ncol(X_stack)){
				keep = qr_res$pivot[seq_len(qr_res$rank)]
				# Always keep intercept(1) and treatment(2)
				if (!(1L %in% keep)) keep = c(1L, keep)
				if (!(2L %in% keep)) keep = c(2L, keep)
				keep = sort(unique(keep))
				# Extract original covariate indices (offset by 2)
				return(keep[keep > 2L] - 2L)
			}
			seq_len(p)
		},

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
			try_combined_fit = function(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v){
				mod = tryCatch(
				fast_cpoisson_combined_with_var_cpp(
					yT_v = as.numeric(yT_v),
					n_k_v = as.numeric(n_k_v),
					X_diff_v = as.matrix(X_diff_v),
					y_r = as.numeric(y_r_v),
					w_r = as.numeric(w_r_v),
					X_r = as.matrix(X_r_v)
				),
				error = function(e) NULL
			)
				if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
				if (!estimate_only && (!is.finite(mod$ssq_b_j) || mod$ssq_b_j < 0)) return(FALSE)

				private$cached_mod = mod
				private$cached_values$likelihood_test_context = list(
					yT_v = as.numeric(yT_v),
					n_k_v = as.numeric(n_k_v),
					X_diff_v = as.matrix(X_diff_v),
					y_r_v = as.numeric(y_r_v),
					w_r_v = as.numeric(w_r_v),
					X_r_v = as.matrix(X_r_v)
				)
				private$cached_values$beta_hat_T = as.numeric(mod$b[2])
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
				private$cached_values$is_z = TRUE
				private$cached_values$df = NA_real_
				TRUE
			},

		try_pairs_only = function(estimate_only, yT_v, n_k_v, X_diff_v){
			if (length(yT_v) == 0L) return(FALSE)
			y_prop = yT_v / n_k_v
			Xmm = if (ncol(X_diff_v) > 0L) cbind(1, X_diff_v) else matrix(1, nrow = length(yT_v), ncol = 1L)
			mod = tryCatch({
				res = fast_logistic_regression_weighted_cpp(X = Xmm, y = y_prop, weights = n_k_v)
				list(b = res$b, ssq_b_1 = NA_real_, X_fit = Xmm) # weights version doesn't return var yet
			}, error = function(e) NULL)
			if (is.null(mod) || length(mod$b) < 1L || !is.finite(mod$b[1])) return(FALSE)
			private$cached_values$beta_hat_T = as.numeric(mod$b[1])
			private$cached_values$is_z = TRUE
			# Standard error remains NA for now as it requires specific Rcpp implementation
			TRUE
		},

		try_reservoir_only = function(estimate_only, y_r_v, w_r_v, X_r_v){
			if (length(y_r_v) == 0L) return(FALSE)
			X_full = if (ncol(X_r_v) > 0L) cbind(1, w_r_v, X_r_v) else cbind(1, w_r_v)
			j_treat = 2L
			if (ncol(X_full) > 2L){
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(2L %in% keep)) keep[r_full] = 2L
					keep = sort(unique(keep))
					X_full = X_full[, keep, drop = FALSE]
					j_treat = which(colnames(X_full) == "w_r_v")
					if (length(j_treat) == 0) j_treat = 2L
				}
			}
			mod = tryCatch(
				fast_neg_bin_with_var_cpp(X = X_full, y = as.integer(y_r_v)),
				error = function(e) NULL
			)
			if (is.null(mod) || length(mod$b) < j_treat || !is.finite(mod$b[j_treat])) return(FALSE)
			if (!estimate_only && (!is.finite(mod$ssq_b_j) || mod$ssq_b_j < 0)) return(FALSE)
			private$cached_values$beta_hat_T = as.numeric(mod$b[j_treat])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z = TRUE
			TRUE
		},

		# The combined case is handled by fast_cpoisson_combined_with_var_cpp
		# (Newton's method with analytic Fisher-information Hessian).
		shared_combined_likelihood = function(estimate_only = FALSE){
				print(paste("DEBUG: Cond Poisson OneLik shared called for", class(self)[1]))
				if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
				if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
				private$cached_values$likelihood_test_context = NULL
				private$cached_mod = NULL

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
				private$cached_values$combined_cov_keep = NULL
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			p             = ncol(as.matrix(private$X))
			has_reservoir = nRT > 0 && nRC > 0

			# ---- Pair data (conditional Poisson, zero-total pairs discarded) ----
			yT_v     = numeric(0)
			n_k_v    = numeric(0)
			X_diff_v = matrix(nrow = 0L, ncol = p)

			if (m > 0){
				yT    = KKstats$yTs_matched
				yC    = KKstats$yCs_matched
				n_k   = yT + yC
				valid = which(n_k > 0)
				if (length(valid) > 0){
					yT_v  = yT[valid]
					n_k_v = n_k[valid]
					if (p > 0L) {
						# Use the full-width pair-difference matrix here so any later
						# covariate reduction stays aligned with the reservoir matrix.
						X_diff_v = as.matrix(KKstats$X_matched_diffs_full[valid, drop = FALSE])
					}
				}
			}
			has_pairs = length(yT_v) > 0

			# ---- Reservoir data (marginal Poisson) ----
			y_r_v = numeric(0)
			w_r_v = numeric(0)
			X_r_v = matrix(nrow = 0L, ncol = p)
			if (has_reservoir){
				y_r_v = KKstats$y_reservoir
				w_r_v = KKstats$w_reservoir
				if (p > 0L) X_r_v = as.matrix(KKstats$X_reservoir)
			}

			if (!has_pairs && !has_reservoir){
				private$set_failed_combined_cache()
				return(invisible(NULL))
			}

			# ---- Covariate Reduction (Aligned) -----------------------------------
			if (p > 0L && is.null(private$cached_values$combined_cov_keep)){
				private$cached_values$combined_cov_keep = private$reduce_combined_covariates(X_diff_v, X_r_v, w_r_v)
			}
			if (p > 0L){
				keep = private$cached_values$combined_cov_keep
				if (length(keep) > 0L){
					X_diff_v = X_diff_v[, keep, drop = FALSE]
					X_r_v    = X_r_v[, keep, drop = FALSE]
				} else {
					X_diff_v = matrix(nrow = nrow(X_diff_v), ncol = 0L)
					X_r_v    = matrix(nrow = nrow(X_r_v), ncol = 0L)
				}
			}

			# ---- Joint Likelihood Fit --------------------------------------------
			success = FALSE
			if (has_pairs && has_reservoir){
				success = private$try_combined_fit(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v)
			}
			
			if (!success){
				# Fallback: one or the other
				fallback_success = FALSE
				if (has_pairs){
					fallback_success = private$try_pairs_only(estimate_only, yT_v, n_k_v, X_diff_v)
				}
				if (!fallback_success && has_reservoir){
					fallback_success = private$try_reservoir_only(estimate_only, y_r_v, w_r_v, X_r_v)
				}
				if (!fallback_success){
					private$set_failed_combined_cache()
				}
			}

			invisible(NULL)
		}
	)
)
