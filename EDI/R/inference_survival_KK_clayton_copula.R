#' Abstract class for Clayton Copula / Standard Weibull Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses using a Clayton copula with Weibull AFT margins for matched
#' pairs and a standard Weibull AFT model for the reservoir. The two treatment-effect
#' estimates (on the log-time ratio scale) are combined by inverse-variance weighting.
#'
#' @keywords internal
InferenceAbstractKKClaytonCopulaIVWC = R6::R6Class("InferenceAbstractKKClaytonCopulaIVWC",
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
		max_abs_reasonable_coef = 1e4,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Re-read design variables which might have been transformed during randomization
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			
			# Recompute basic match data for the new w/y/dead
			private$compute_basic_match_data()
			
			# Clear cached design candidates to allow recalculation if needed
			private$cached_values$clayton_design_candidates = NULL
			
			# Use the same joint-likelihood logic for the point estimate
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		filtered_covariate_candidates = function(X = as.matrix(private$X)){
			if (ncol(X) == 0L) return(list(matrix(nrow = nrow(X), ncol = 0L)))
			
			# Ensure no linearly dependent columns first
			res = drop_linearly_dependent_cols(X)
			X_reduced = res$M
			if (ncol(X_reduced) == 0L) return(list(matrix(nrow = nrow(X), ncol = 0L)))
			
			# Generate candidates by dropping highly correlated columns at various thresholds
			thresholds = c(0.99, 0.9, 0.7)
			candidates = list(X_reduced)
			for (thresh in thresholds) {
				X_try = drop_highly_correlated_cols(X_reduced, threshold = thresh)$M
				# Simple way to avoid duplicates in the list
				is_new = TRUE
				for (c in candidates) {
					if (ncol(c) == ncol(X_try) && all(colnames(c) == colnames(X_try))) {
						is_new = FALSE
						break
					}
				}
				if (is_new) candidates[[length(candidates) + 1L]] = X_try
			}
			candidates
		},

		design_matrix_candidates = function(){
			if (!is.null(private$cached_values$clayton_design_candidates)){
				return(private$cached_values$clayton_design_candidates)
			}

			candidates = list()
			# Case 1: no covariates
			if (ncol(as.matrix(private$X)) == 0L){
				Xcand = matrix(private$w, ncol = 1)
				colnames(Xcand) = "w"
				candidates = list(Xcand)
			} else {
				cov_candidates = private$filtered_covariate_candidates()
				for (X in cov_candidates){
					M = matrix(private$w, ncol = 1)
					colnames(M) = "w"
					if (ncol(X) > 0){
						M = cbind(M, X)
					}
					# Ensure we drop any additional linear dependencies
					qr_res = qr(M)
					if (qr_res$rank < ncol(M)){
						keep = qr_res$pivot[seq_len(qr_res$rank)]
						if (!(1L %in% keep)) keep = c(1L, keep) # Keep treatment
						keep = sort(unique(keep))
						M = M[, keep, drop = FALSE]
					}
					candidates[[length(candidates) + 1L]] = M
				}
			}

			private$cached_values$clayton_design_candidates = candidates
			candidates
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			if (m > 0){
				private$clayton_copula_for_matched_pairs(estimate_only = estimate_only)
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
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
		},

		clayton_copula_for_matched_pairs = function(estimate_only = FALSE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0L)
			if (length(i_matched) == 0L) return(invisible(NULL))

			X_data = private$get_X()
			X_matched = X_data[i_matched, , drop = FALSE]

			# Fit using optimized helper
			fit = .fit_clayton_weibull_aft(
				y = private$y[i_matched],
				dead = private$dead[i_matched],
				X = cbind(w = private$w[i_matched], X_matched),
				pair_id = m_vec[i_matched],
				estimate_only = estimate_only,
				optimization_alg = private$optimization_alg
			)
			
			if (!is.null(fit) && is.finite(fit$beta)){
				private$cached_values$beta_T_matched = fit$beta
				private$cached_values$ssq_beta_T_matched = if (estimate_only) NA_real_ else fit$ssq
			}
		},

		weibull_for_reservoir = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			y_r    = KKstats$y_reservoir
			w_r    = KKstats$w_reservoir
			m_vec_safe = private$m
			if (is.null(m_vec_safe)) m_vec_safe = rep(0L, private$n)
			m_vec_safe[is.na(m_vec_safe)] = 0L
			dead_r = private$dead[m_vec_safe == 0]
			X_r    = as.matrix(KKstats$X_reservoir)

			# Candidate reduced design matrices for the reservoir
			candidates = list()
			if (ncol(X_r) == 0L){
				Xcand = matrix(w_r, ncol = 1)
				colnames(Xcand) = "w"
				candidates = list(Xcand)
			} else {
				cov_candidates = private$filtered_covariate_candidates(X_r)
				for (X in cov_candidates){
					M = matrix(w_r, ncol = 1)
					colnames(M) = "w"
					if (ncol(X) > 0){
						M = cbind(M, X)
					}
					qr_res = qr(M)
					if (qr_res$rank < ncol(M)){
						keep = qr_res$pivot[seq_len(qr_res$rank)]
						if (!(1L %in% keep)) keep = c(1L, keep)
						keep = sort(unique(keep))
						M = M[, keep, drop = FALSE]
					}
					candidates[[length(candidates) + 1L]] = M
				}
			}

			for (Xcand in candidates){
				fit = .fit_standard_weibull_aft_from_matrix(
					y = y_r,
					dead = dead_r,
					X = Xcand,
					estimate_only = estimate_only
				)
				if (!is.null(fit) && is.finite(fit$beta) && (estimate_only || (is.finite(fit$ssq) && fit$ssq > 0))){
					private$cached_values$beta_T_reservoir = fit$beta
					private$cached_values$ssq_beta_T_reservoir = fit$ssq
					return(invisible(NULL))
				}
			}
		}
	)
)

#' Abstract class for Clayton Copula Combined-Likelihood Inference
#'
#' This class implements a combined-likelihood estimator for KK matching-on-the-fly
#' designs with survival responses using a Clayton copula with Weibull AFT margins.
#' The likelihood includes contributions from both matched pairs (bivariate) and
#' reservoir subjects (univariate).
#'
#' @keywords internal
InferenceAbstractKKClaytonCopulaOneLik = R6::R6Class("InferenceAbstractKKClaytonCopulaOneLik",
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

		best_X_colnames = NULL,
		best_par = NULL,
		optimization_alg = "lbfgs",

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Re-read design variables which might have been transformed during randomization
			private$w = private$des_obj_priv_int$w
			private$y = private$des_obj_priv_int$y
			private$dead = private$des_obj_priv_int$dead
			
			# Recompute basic match data for the new w/y/dead
			private$compute_basic_match_data()
			
			# Use the same joint-likelihood logic for the point estimate
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		get_standard_error = function(){
			private$shared()
			as.numeric(private$cached_values$s_beta_hat_T)
		},

		get_degrees_of_freedom = function() Inf,

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		supports_likelihood_tests = function(){
			TRUE
		},

		get_likelihood_test_spec = function(){
			private$shared(estimate_only = FALSE)
			ctx = private$cached_values$likelihood_test_context
			if (is.null(ctx) || is.null(private$cached_mod)) return(NULL)
			X = ctx$X
			y = as.numeric(ctx$y)
			dead = as.numeric(ctx$dead)
			pair_idx = ctx$pair_idx
			singleton_rows = ctx$singleton_rows
			j_treat = as.integer(ctx$j_treat %||% 1L)
			start_len = length(ctx$start %||% numeric(0))
			list(
				X = X,
				y = y,
				j = j_treat,
				full_fit = private$cached_mod$best_fit %||% private$cached_mod,
				fit_null = function(delta, start = NULL){
					start_par = start %||% private$get_fit_warm_start_for_length("params", start_len)
					if (length(start_par) == 0L) start_par = ctx$start
					fast_clayton_weibull_aft_optim_cpp(
						X = X,
						y = y,
						dead = dead,
						pair_idx = pair_idx,
						singleton_rows = singleton_rows,
						start_params = start_par,
						estimate_only = FALSE,
						optimization_alg = private$optimization_alg,
						fixed_idx = j_treat,
						fixed_values = delta
					)
				},
				extract_start = function(fit){
					as.numeric(fit$params %||% fit$par %||% fit$best_par)
				},
				score = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.numeric(fit$score %||% get_clayton_weibull_aft_score_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				observed_information = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.matrix(fit$observed_information %||% fit$information %||% -get_clayton_weibull_aft_hessian_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				fisher_information = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.matrix(fit$information %||% fit$observed_information %||% -get_clayton_weibull_aft_hessian_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				information = function(fit){
					params = as.numeric(fit$params %||% fit$par %||% fit$best_par)
					as.matrix(fit$information %||% fit$observed_information %||% -get_clayton_weibull_aft_hessian_cpp(X, y, dead, pair_idx, singleton_rows, params))
				},
				neg_loglik = function(fit){
					as.numeric(fit$neg_loglik %||% fit$neg_ll %||% fit$value)
				}
			)
		},

		filtered_covariate_candidates = function(){
			X = as.matrix(private$X)
			if (ncol(X) == 0L) return(list(matrix(nrow = private$n, ncol = 0L)))
			
			# Ensure no linearly dependent columns first
			res = drop_linearly_dependent_cols(X)
			X_reduced = res$M
			if (ncol(X_reduced) == 0L) return(list(matrix(nrow = private$n, ncol = 0L)))
			
			# Generate candidates by dropping highly correlated columns at various thresholds
			thresholds = c(0.99, 0.9, 0.7)
			candidates = list(X_reduced)
			for (thresh in thresholds) {
				X_try = drop_highly_correlated_cols(X_reduced, threshold = thresh)$M
				# Simple way to avoid duplicates in the list
				is_new = TRUE
				for (c in candidates) {
					if (ncol(c) == ncol(X_try) && all(colnames(c) == colnames(X_try))) {
						is_new = FALSE
						break
					}
				}
				if (is_new) candidates[[length(candidates) + 1L]] = X_try
			}
			candidates
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Optimization candidates
			w_only = matrix(private$w, ncol = 1); colnames(w_only) = "w"
			if (ncol(as.matrix(private$X)) == 0L){
				candidates = list(w_only)
			} else {
				cov_candidates = private$filtered_covariate_candidates()
				candidates = list()
				for (X_cov in cov_candidates){
					if (ncol(X_cov) == 0L){
						M = w_only
					} else {
						M = matrix(private$w, ncol = 1)
						colnames(M) = "w"
						M = cbind(M, X_cov)
						qr_M = qr(M)
						if (qr_M$rank < ncol(M)){
							keep = qr_M$pivot[seq_len(qr_M$rank)]
							if (!(1L %in% keep)) keep = c(1L, keep)
							keep = sort(unique(keep))
							M = M[, keep, drop = FALSE]
						}
						colnames(M)[1] = "w"
					}
					candidates[[length(candidates) + 1L]] = M
				}
				candidates[[length(candidates) + 1L]] = w_only
			}

			for (X in candidates){
				fit = .fit_clayton_weibull_aft(
					y = private$y,
					dead = private$dead,
					X = X,
					pair_id = m_vec,
					include_singletons = TRUE,
					estimate_only = estimate_only
				)
				if (!is.null(fit) && is.finite(fit$beta) && (isTRUE(estimate_only) || (is.finite(fit$ssq) && fit$ssq > 0))){
					private$cached_values$beta_hat_T = fit$beta
					private$cached_values$s_beta_hat_T = if (is.finite(fit$ssq) && fit$ssq > 0) sqrt(fit$ssq) else NA_real_
					private$cached_values$theta = fit$theta
					private$best_par = fit$best_par
					private$best_X_colnames = colnames(X)
					private$cached_mod = fit
					pair_idx = .complete_pair_index_matrix(m_vec) - 1L
					singleton_rows = if (length(pair_idx) > 0L) setdiff(seq_len(length(private$y)), sort(unique(as.vector(pair_idx + 1L)))) - 1L else seq_len(length(private$y)) - 1L
					private$cached_values$likelihood_test_context = list(
						X = X,
						y = as.numeric(private$y),
						dead = as.numeric(private$dead),
						pair_idx = pair_idx,
						singleton_rows = singleton_rows,
						start = fit$best_par,
						j_treat = 1L
					)
					return(invisible(NULL))
				}
			}

			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
		}
	)
)

#' Clayton Copula IVWC Compound Inference for KK Designs
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with survival responses
#' using a Clayton copula with Weibull AFT margins for matched pairs and standard
#' Weibull AFT regression for reservoir subjects.
#'
#' @examples
#' \dontrun{
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKClaytonCopulaOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' }
#' @export
InferenceSurvivalKKClaytonCopulaIVWC = R6::R6Class("InferenceSurvivalKKClaytonCopulaIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClaytonCopulaIVWC,
	public = list(
	)
)
InferenceIncidKKClaytonCopulaIVWC = InferenceSurvivalKKClaytonCopulaIVWC

#' Clayton Copula Combined-Likelihood Inference for KK Designs
#'
#' Fits a joint copula-based likelihood for KK matching-on-the-fly designs with
#' survival responses using a Clayton copula with Weibull AFT margins.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'survival')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(runif(10))
#' inf = InferenceSurvivalKKClaytonCopulaOneLik$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceSurvivalKKClaytonCopulaOneLik = R6::R6Class("InferenceSurvivalKKClaytonCopulaOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClaytonCopulaOneLik,
	public = list(
	)
)
InferenceIncidKKClaytonCopulaOneLik = InferenceSurvivalKKClaytonCopulaOneLik
