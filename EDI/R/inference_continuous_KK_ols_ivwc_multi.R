#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
#' @export
InferenceContinMultOLSKKIVWC = R6::R6Class("InferenceContinMultOLSKKIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompound,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{FALSE}.
		#' @param ols_sanity_kappa_max Numeric scalar. Maximum acceptable condition number for the
		#'   OLS design matrix (computed via \code{kappa(..., exact = FALSE)}). When the condition
		#'   number exceeds this threshold the design matrix is considered near-singular and the
		#'   estimator falls back to the simple mean-difference formula without fitting OLS at all.
		#'   Default: \code{1e6}.
		#' @param ols_sanity_var_max_factor Numeric scalar \eqn{\geq 1}. The OLS variance estimate
		#'   must be strictly less than \code{ols_sanity_var_max_factor} times the simple
		#'   (variance-of-mean-difference) variance. Values larger than this indicate a near-singular
		#'   \eqn{X^\top X} whose inverse is inflated. Default: \code{10}.
		#' @param ols_sanity_var_min_factor Numeric scalar in \eqn{(0, 1)}. The OLS variance estimate
		#'   must be strictly greater than \code{ols_sanity_var_min_factor} times the simple variance.
		#'   Values smaller than this indicate \eqn{R^2 \approx 1}, meaning OLS has absorbed the
		#'   treatment effect into a covariate and the estimate is degenerate. Default: \code{1e-6}.
		#' @param ols_sanity_se_max_factor Numeric scalar \eqn{\geq 0}. The OLS point estimate must
		#'   lie within \code{ols_sanity_se_max_factor} simple standard errors of the simple
		#'   mean-difference estimate. Larger deviations suggest the treatment column is confounded
		#'   with a near-collinear covariate. Default: \code{5}.
		initialize = function(des_obj, verbose = FALSE,
				ols_sanity_kappa_max    = 1e6,
				ols_sanity_var_max_factor = 10,
				ols_sanity_var_min_factor = 1e-6,
				ols_sanity_se_max_factor  = 5){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "continuous")
			}
			if (should_run_asserts()) {
				assertNumber(ols_sanity_kappa_max,      lower = 1)
				assertNumber(ols_sanity_var_max_factor, lower = 1)
				assertNumber(ols_sanity_var_min_factor, lower = 0, upper = 1)
				assertNumber(ols_sanity_se_max_factor,  lower = 0)
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$ols_sanity_kappa_max      = ols_sanity_kappa_max
			private$ols_sanity_var_max_factor = ols_sanity_var_max_factor
			private$ols_sanity_var_min_factor = ols_sanity_var_min_factor
			private$ols_sanity_se_max_factor  = ols_sanity_se_max_factor
		},

		#' @description
		#' Computes the appropriate estimate
		#'
		#' @return	The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @examples
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#'
		#' seq_des_inf = InferenceContinMultOLS$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#'
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			if (is.null(private$cached_values$beta_T_reservoir) & is.null(private$cached_values$beta_T_matched)){
				private$compute_estimate_from_matched_and_reservoir(
					private$ols_for_matched_pairs,
					private$ols_for_reservoir
				)
			}
			if (is.null(private$cached_values[["beta_hat_T"]])){
				beta_hat_T = if (private$only_matches()){
								private$cached_values$beta_T_matched
							} else if (private$only_reservoir()){
								private$cached_values$beta_T_reservoir
							} else {
								# If one side has an unusable variance (NA/non-positive from a
								# rank-deficient OLS or too few matched pairs), fall back to the
								# other side's estimate so the result stays finite.
								ssq_m = private$cached_values$ssq_beta_T_matched
								ssq_r = private$cached_values$ssq_beta_T_reservoir
								if (!is.finite(ssq_m) || ssq_m <= 0) {
									private$cached_values$beta_T_reservoir
								} else if (!is.finite(ssq_r) || ssq_r <= 0) {
									private$cached_values$beta_T_matched
								} else {
									w_star = ssq_r / (ssq_r + ssq_m)
									w_star * private$cached_values$beta_T_matched + (1 - w_star) * private$cached_values$beta_T_reservoir
								}
							}
				# Avoid floating-point residue when the estimate is mathematically zero.
				if (is.finite(beta_hat_T) && abs(beta_hat_T) < sqrt(.Machine$double.eps)){
					beta_hat_T = 0
				}
				private$cached_values[["beta_hat_T"]] = beta_hat_T
			}
			private$cached_values[["beta_hat_T"]]
		},



	#' @description
	#' Computes a 1-alpha level frequentist confidence interval
	#' differently for all response types, estimate types, and
	#' test types.
	#'
	#' Here we use the theory that MLE's computed for GLM's are asymptotically normal.
	#' Hence these confidence intervals are asymptotically valid
	#' and thus approximate for any sample size.
	#'
	#' @param alpha The confidence level in the computed confidence
	#'   interval is 1 - \code{alpha}. The default is 0.05.
	#'
	#' @return	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
	#'
	#' @examples
	#' \dontrun{
	#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
	#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
	#'
	#' seq_des_inf = InferenceContinMultOLS$new(seq_des)
	#' seq_des_inf$compute_asymp_confidence_interval()
	#' }
	#'
	compute_asymp_confidence_interval = function(alpha = 0.05){
		if (should_run_asserts()) {
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
		}
		private$shared_for_inference()
		private$compute_z_or_t_ci_from_s_and_df(alpha)
	},

	#' @description
	#' Computes a 2-sided p-value
	#'
	#' @param delta The null difference to test against. For any
	#'   treatment effect at all this is set to zero (the default).
	#'
	#' @return	The approximate frequentist p-value
	#'
	#' @examples
	#' \dontrun{
	#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
	#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
	#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
	#'
	#' seq_des_inf = InferenceContinMultOLS$new(seq_des)
	#' seq_des_inf$compute_asymp_two_sided_pval_for_treatment_effect()
	#' }
	#'
	compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
		if (should_run_asserts()) {
			assertNumeric(delta)
		}
		private$shared_for_inference()
		private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
	}
	),

	private = list(
		ols_sanity_kappa_max      = 1e6,
		ols_sanity_var_max_factor = 10,
		ols_sanity_var_min_factor = 1e-6,
		ols_sanity_se_max_factor  = 5,

		shared_for_inference = function(){
			if (is.null(private$cached_values[["beta_hat_T"]])){
				self$compute_treatment_estimate()
			}
			ssq_beta_hat_T = if (private$only_matches()){
								ssq_m = private$cached_values$ssq_beta_T_matched
								if (!is.finite(ssq_m) || ssq_m <= 0) {
									# Fallback: use reservoir SE even though reservoir is small
									if (is.null(private$cached_values$ssq_beta_T_reservoir)) {
										private$ols_for_reservoir()
									}
									private$cached_values$ssq_beta_T_reservoir
								} else {
									ssq_m
								}
							} else if (private$only_reservoir()){
								ssq_r = private$cached_values$ssq_beta_T_reservoir
								if (!is.finite(ssq_r) || ssq_r <= 0) {
									# Fallback: use matched-pairs SE
									if (is.null(private$cached_values$ssq_beta_T_matched) &&
										private$cached_values$KKstats$m > 0) {
										private$ols_for_matched_pairs()
									}
									private$cached_values$ssq_beta_T_matched
								} else {
									ssq_r
								}
							} else {
								# Combined: if one side is unusable, fall back to the other
								ssq_m = private$cached_values$ssq_beta_T_matched
								ssq_r = private$cached_values$ssq_beta_T_reservoir
								if (!is.finite(ssq_m) || ssq_m <= 0) {
									ssq_r
								} else if (!is.finite(ssq_r) || ssq_r <= 0) {
									ssq_m
								} else {
									ssq_m * ssq_r / (ssq_m + ssq_r) #analogous eq's
								}
							}
			private$cached_values$s_beta_hat_T = sqrt(ssq_beta_hat_T)
			private$cached_values$is_z = TRUE #TO-DO: linear combination of degrees of freedom of t's
		},

		only_matches = function(){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(FALSE)
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (!is.finite(nRT) || !is.finite(nRC)) return(FALSE)
			nRT <= 2 || nRC <= 2 || (nRT + nRC <= ncol(private$get_X()) + 2)
		},

		only_reservoir = function(){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(FALSE)
			m = KKstats$m
			is.finite(m) && m <= 1
		},

		ols_for_matched_pairs = function(){
			yd = private$cached_values$KKstats$y_matched_diffs
			Xd = private$cached_values$KKstats$X_matched_diffs
			m  = length(yd)

			# Use QR with column pivoting to reduce Xd to a numerically full-rank subset.
			# This handles: (a) exact linear dependencies from factor dummy columns that
			# sum to zero in the differences, and (b) near-zero columns that would make
			# X'X near-singular and cause ConjugateGradient to fail in fast_ols_with_var_cpp.
			if (ncol(Xd) > 0) {
				qr_Xd = qr(Xd)
				r = qr_Xd$rank
				if (r < ncol(Xd)) {
					Xd = Xd[, qr_Xd$pivot[seq_len(r)], drop = FALSE]
				}
			}
			p_kept = ncol(Xd)

			# If the system is underdetermined after rank reduction (or no covariates),
			# fall back to the simple mean difference (intercept-only OLS).
			if (m <= p_kept + 1) {
				private$cached_values$beta_T_matched     = if (m >= 1) mean(yd) else NA_real_
				private$cached_values$ssq_beta_T_matched = if (m >= 2) var(yd) / m else NA_real_
				return(invisible(NULL))
			}

			design_mat = cbind(1, Xd)
			# Guard against near-singular bootstrap samples: duplicate rows reduce the
			# effective rank, making X'X ill-conditioned even after QR rank reduction.
			# LDLT on a near-singular X'X returns extreme but finite coefficients.
			# Fall back to simple mean difference when the design matrix is ill-conditioned.
			if (kappa(design_mat, exact = FALSE) > private$ols_sanity_kappa_max) {
				private$cached_values$beta_T_matched     = if (m >= 1) mean(yd) else NA_real_
				private$cached_values$ssq_beta_T_matched = if (m >= 2) var(yd) / m else NA_real_
				return(invisible(NULL))
			}

			beta_simple  = if (m >= 1) mean(yd) else NA_real_
			fallback_ssq = if (m >= 2) var(yd) / m else NA_real_
			mod = fast_ols_with_var_cpp(design_mat, yd, j = 1) #the only time you need the intercept's ssq
			# Guard against failure modes in ill-conditioned/degenerate bootstrap samples:
			# (1) OLS variance inflated (>ols_sanity_var_max_factor x simple variance): near-singular X'X.
			# (2) OLS estimate far from simple mean diff (>ols_sanity_se_max_factor SEs): treatment
			#     confounded with a near-perfect covariate predictor, giving extreme estimate.
			# (3) OLS variance near-zero (<ols_sanity_var_min_factor of simple variance): R^2 approx 1,
			#     which means OLS has absorbed treatment into a covariate, inflating the coefficient.
			#     Near-zero ssq_b_j also makes w_star -> 0 in the compound formula,
			#     giving all weight to the (degenerate) reservoir OLS estimate.
			se_simple = if (!is.na(fallback_ssq) && fallback_ssq > 0) sqrt(fallback_ssq) else NA_real_
			ols_ok = is.finite(mod$ssq_b_j) &&
			         (is.na(fallback_ssq) || (mod$ssq_b_j < private$ols_sanity_var_max_factor * fallback_ssq &&
			                                   mod$ssq_b_j > private$ols_sanity_var_min_factor * fallback_ssq)) &&
			         (is.na(se_simple)    || abs(mod$b[1] - beta_simple) <= private$ols_sanity_se_max_factor * se_simple)
			if (ols_ok) {
				private$cached_values$beta_T_matched     = mod$b[1]
				private$cached_values$ssq_beta_T_matched = mod$ssq_b_j
			} else {
				private$cached_values$beta_T_matched     = beta_simple
				private$cached_values$ssq_beta_T_matched = fallback_ssq
			}
		},

		ols_for_reservoir = function(){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			X_r = as.matrix(private$cached_values$KKstats$X_reservoir)

			# Build full design matrix: intercept (col 1), treatment (col 2), covariates (cols 3+).
			X_full  = cbind(1, w_r, X_r)
			j_treat = 2L

			# QR-reduce to full rank while always preserving the treatment column.
			# This handles datasets where model.matrix(~ 0 + .) includes all factor levels,
			# which when combined with the intercept creates a linear dependency.
			qr_full = qr(X_full)
			r_full  = qr_full$rank
			if (r_full < ncol(X_full)) {
				keep = qr_full$pivot[seq_len(r_full)]
				# Ensure treatment column is always included
				if (!(2L %in% keep)) {
					keep[r_full] = 2L
				}
				keep    = sort(keep)
				X_full  = X_full[, keep, drop = FALSE]
				j_treat = which(keep == 2L)
			}

			# Guard against near-singular bootstrap samples: duplicate rows reduce the
			# effective rank, making X'X ill-conditioned even after QR rank reduction.
			# LDLT on a near-singular X'X returns extreme but finite coefficients.
			# Fall back to simple mean difference when the design matrix is ill-conditioned.
			if (kappa(X_full, exact = FALSE) > private$ols_sanity_kappa_max) {
				y_rT = y_r[w_r == 1]; y_rC = y_r[w_r == 0]
				nRT_local = length(y_rT); nRC_local = length(y_rC)
				private$cached_values$beta_T_reservoir = mean(y_rT) - mean(y_rC)
				private$cached_values$ssq_beta_T_reservoir =
					if (nRT_local > 1 && nRC_local > 1) var(y_rT)/nRT_local + var(y_rC)/nRC_local else NA_real_
				return(invisible(NULL))
			}

			y_rT = y_r[w_r == 1]; y_rC = y_r[w_r == 0]
			nRT_local = length(y_rT); nRC_local = length(y_rC)
			beta_simple  = mean(y_rT) - mean(y_rC)
			fallback_ssq = if (nRT_local > 1 && nRC_local > 1) var(y_rT)/nRT_local + var(y_rC)/nRC_local else NA_real_

			mod = fast_ols_with_var_cpp(X_full, y_r, j = j_treat)
			# Guard against failure modes in ill-conditioned/degenerate bootstrap samples:
			# (1) OLS variance inflated (>ols_sanity_var_max_factor x simple variance): near-singular X'X.
			# (2) OLS estimate far from simple mean diff (>ols_sanity_se_max_factor SEs): treatment
			#     confounded with a near-perfect covariate predictor, giving extreme estimate.
			# (3) OLS variance near-zero (<ols_sanity_var_min_factor of simple variance): R^2 approx 1,
			#     which means OLS has absorbed treatment into a covariate, inflating the coefficient.
			#     Near-zero ssq_b_j also makes w_star -> 0 in the compound formula,
			#     giving all weight to the (degenerate) reservoir OLS estimate.
			se_simple = if (!is.na(fallback_ssq) && fallback_ssq > 0) sqrt(fallback_ssq) else NA_real_
			ols_ok = is.finite(mod$ssq_b_j) &&
			         (is.na(fallback_ssq) || (mod$ssq_b_j < private$ols_sanity_var_max_factor * fallback_ssq &&
			                                   mod$ssq_b_j > private$ols_sanity_var_min_factor * fallback_ssq)) &&
			         (is.na(se_simple)    || abs(mod$b[j_treat] - beta_simple) <= private$ols_sanity_se_max_factor * se_simple)
			if (ols_ok) {
				private$cached_values$beta_T_reservoir     = mod$b[j_treat]
				private$cached_values$ssq_beta_T_reservoir = mod$ssq_b_j
			} else {
				private$cached_values$beta_T_reservoir     = beta_simple
				private$cached_values$ssq_beta_T_reservoir = fallback_ssq
			}
		}
	)
)
