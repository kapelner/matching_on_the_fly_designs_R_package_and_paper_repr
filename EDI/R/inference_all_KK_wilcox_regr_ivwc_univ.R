#' Univariate Wilcox Rank-based Regression Compound Inference for KK Designs
#'
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with only the treatment indicator (no additional covariates).
#' For matched pairs, it uses R-estimation on the within-pair response differences.
#' For reservoir subjects, it uses a standard rank-based linear model with the treatment
#' assignment as the sole predictor. This is the univariate (covariate-free) variant.
#'
#' @export
InferenceAllKKWilcoxRegrUnivIVWC = R6::R6Class("InferenceAllKKWilcoxRegrUnivIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			res_type = des_obj$get_response_type()
			if (should_run_asserts()) {
				if (res_type == "incidence"){
					stop("Rank-based regression is not recommended for incidence data; clogit and compound mean diff is recommended.")
				}
			}
			if (should_run_asserts()) {
				assertResponseType(res_type, c("continuous", "count", "proportion", "survival", "ordinal"))
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				if (private$any_censoring){
					stop(class(self)[1], " does not support censored survival data.")
				}
				if (!check_package_installed("Rfit")) {
					stop("Package 'Rfit' is required for ", class(self)[1], ". Please install it.")
				}
			}
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha The confidence level is 1 - \code{alpha}. The default is 0.05.
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
		#' Computes the asymptotic two-sided p-value.
		#' @param delta The null difference to test against. Default is zero.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				if (should_run_asserts()) {
					stop("Testing non-zero delta is not yet implemented for the combined rank-regression estimator.")
				}
				NA_real_
			}
		}
	),

	private = list(
		best_Xmm_colnames_matched = NULL,
		best_Xmm_colnames_reservoir = NULL,
		best_Xmm_j_treat_reservoir = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames_matched) && is.null(private$best_Xmm_colnames_reservoir)){
				private$shared(estimate_only = TRUE)
			}

			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs component ---
			beta_m = NA_real_
			if (m > 0){
				dy = KKstats$y_matched_diffs
				if (private$include_covariates() && !is.null(private$best_Xmm_colnames_matched)){
					dX = KKstats$X_matched_diffs[, intersect(private$best_Xmm_colnames_matched, colnames(KKstats$X_matched_diffs)), drop = FALSE]
					Xmm = cbind(`(Intercept)` = 1, dX)
					beta_m = private$fit_rank_regression_exact(dy, Xmm, 1L, estimate_only = TRUE)$beta
				} else {
					beta_m = wilcox_hl_signed_rank_point_estimate_cpp(as.numeric(dy))
				}
			}

			# --- Reservoir component ---
			beta_r = NA_real_
			if (nRT > 0 && nRC > 0){
				y_r = KKstats$y_reservoir
				w_r = KKstats$w_reservoir
				if (private$include_covariates() && !is.null(private$best_Xmm_colnames_reservoir)){
					X_r = as.matrix(KKstats$X_reservoir)
					X_cov = X_r[, intersect(private$best_Xmm_colnames_reservoir, colnames(X_r)), drop = FALSE]
					Xmm = cbind(`(Intercept)` = 1, w = w_r, X_cov)
					j_treat = private$best_Xmm_j_treat_reservoir %||% 2L
					beta_r = private$fit_rank_regression_exact(y_r, Xmm, j_treat, estimate_only = TRUE)$beta
				} else {
					beta_r = wilcox_hl_point_estimate_cpp(as.numeric(y_r), as.integer(w_r))
				}
			}

			combo = private$combine_component_estimates(
				beta_m = beta_m,
				ssq_m = private$cached_values$ssq_beta_T_matched, # use orig ssq for weights
				beta_r = beta_r,
				ssq_r = private$cached_values$ssq_beta_T_reservoir
			)
			combo$beta
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			# Preserve original SSQ values across rand worker cache resets so that
			# combine_component_estimates() can use them for IVWC variance-weighting.
			private$compute_fast_randomization_distr_via_reused_worker(
				y, permutations, delta, transform_responses,
				preserve_cache_keys = c("ssq_beta_T_matched", "ssq_beta_T_reservoir"),
				zero_one_logit_clamp = zero_one_logit_clamp
			)
		},

		include_covariates = function(){
			FALSE
		},

		rfit_internal_fns = NULL,

		get_rfit_internal_fns = function(){
			if (is.null(private$rfit_internal_fns)){
				private$rfit_internal_fns = list(
					disp = getFromNamespace("disp", "Rfit"),
					jaeckel = getFromNamespace("jaeckel", "Rfit"),
					gettauF0 = getFromNamespace("gettauF0", "Rfit"),
					taustar = getFromNamespace("taustar", "Rfit")
				)
			}
			private$rfit_internal_fns
		},

		fit_rank_regression_exact = function(y, X, coef_index, estimate_only = FALSE){
			tryCatch({
				x = as.matrix(X)
				storage.mode(x) = "double"
				y = as.numeric(y)
				if (nrow(x) != length(y) || ncol(x) < coef_index){
					return(list(beta = NA_real_, ssq = NA_real_))
				}
				if (is.null(colnames(x))){
					colnames(x) = c("(Intercept)", paste0("x", seq_len(ncol(x) - 1L)))
				}
				if (abs(max(x) - min(x)) < .Machine$double.eps^0.5){
					return(list(beta = NA_real_, ssq = NA_real_))
				}

				# Reduce any collinear columns while keeping the treatment coefficient.
				# This mirrors the QR-based guard used by the continuous robust-regression
				# KK estimator and prevents singular X from reaching chol()/lsfit().
				qr_X = qr(x)
				if (qr_X$rank < ncol(x)){
					keep = qr_X$pivot[seq_len(qr_X$rank)]
					if (!(coef_index %in% keep)) keep[qr_X$rank] = coef_index
					keep = sort(unique(keep))
					x = x[, keep, drop = FALSE]
					coef_index = which(keep == coef_index)
				}

				fns = private$get_rfit_internal_fns()
				scores = Rfit::wscores
				x1 = as.matrix(x[, colnames(x) != "(Intercept)", drop = FALSE])
				x1 = as.matrix(cbind(`(Intercept)` = rep(1, nrow(x1)), x1))
				qrx = qr(x1)
				Q = as.matrix(qr.Q(qrx))
				q1 = Q[, 1, drop = FALSE]
				xq = if (qrx$rank > 1L) as.matrix(Q[, 2:qrx$rank, drop = FALSE]) else matrix(0, nrow = nrow(x), ncol = 0)
				betahat0 = if (ncol(xq) > 0L) drop(crossprod(xq, y)) else numeric(0)
				if (length(betahat0) > 0L &&
					fns$disp(betahat0, xq, y, scores) > fns$disp(rep(0, length(betahat0)), xq, y, scores)){
					betahat0 = rep(0, length(betahat0))
				}

				if (length(betahat0) > 0L){
					ord = order(y - xq %*% betahat0)
					fit = fns$jaeckel(as.matrix(xq[ord, , drop = FALSE]), y[ord], betahat0, scores = scores)
					if (fit$convergence != 0){
						ord = order(y - xq %*% fit$par)
						fit2 = fns$jaeckel(as.matrix(xq[ord, , drop = FALSE]), y[ord], jitter(fit$par), scores = scores)
						if (fit$convergence != 0){
							warning("rfit: Convergence status not zero in jaeckel")
							if (fit2$value < fit$value){
								fit = fit2
							}
						} else {
							fit = fit2
						}
					}
					betahat = fit$par
					yhat = drop(xq %*% betahat)
				} else {
					yhat = rep(0, length(y))
				}

				ehat = y - yhat
				alphahat = stats::median(ehat)
				ehat = ehat - alphahat
				yhat = yhat + alphahat
				bhat = stats::lsfit(x, yhat, intercept = FALSE)$coef
				alphahat0 = stats::median(y)
				bhat0 = c(alphahat0, rep(0, length(bhat) - 1L))
				if (fns$disp(bhat, x, y, scores) > fns$disp(bhat0, x, y, scores)){
					bhat = bhat0
					ehat = y - alphahat0
				}

				beta = unname(as.numeric(bhat[coef_index]))
				if (!is.finite(beta)) return(list(beta = NA_real_, ssq = NA_real_))
				if (estimate_only) return(list(beta = beta, ssq = NA_real_))

				tauhat = fns$gettauF0(ehat, ncol(xq), scores)
				taushat = fns$taustar(ehat, qrx$rank)
				if (coef_index > length(bhat)){
					return(list(beta = NA_real_, ssq = NA_real_))
				}

				xxpxi = x %*% chol2inv(chol(crossprod(x)))
				A1 = crossprod(q1, xxpxi)
				vcov_mat = taushat^2 * crossprod(A1)
				if (qrx$rank > 1L){
					Q2 = Q[, 2:qrx$rank, drop = FALSE]
					A2 = crossprod(Q2, xxpxi)
					vcov_mat = vcov_mat + tauhat^2 * crossprod(A2)
				}

				ssq = unname(as.numeric(vcov_mat[coef_index, coef_index]))
				if (!is.finite(ssq) || ssq <= 0){
					return(list(beta = NA_real_, ssq = NA_real_))
				}
				list(beta = beta, ssq = ssq)
			}, error = function(e){
				list(beta = NA_real_, ssq = NA_real_)
			})
		},

		combine_component_estimates = function(beta_m, ssq_m, beta_r, ssq_r){
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
			       !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0
			r_ok = !is.null(beta_r) && is.finite(beta_r) &&
			       !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				return(list(
					beta = w_star * beta_m + (1 - w_star) * beta_r,
					se = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
				))
			}
			if (m_ok){
				return(list(beta = beta_m, se = sqrt(ssq_m)))
			}
			if (r_ok){
				return(list(beta = beta_r, se = sqrt(ssq_r)))
			}
			list(beta = NA_real_, se = NA_real_)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			# Recompute KKstats if cache was cleared (e.g., after y transformation for rand CI)
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: Rfit on differences ---
			if (m > 0){
				private$rfit_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched

			# --- Reservoir: Rfit on level data ---
			if (nRT > 0 && nRC > 0){
				private$rfit_for_reservoir(estimate_only = estimate_only)
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir

			combo = private$combine_component_estimates(beta_m, ssq_m, beta_r, ssq_r)
			private$cached_values$beta_hat_T = combo$beta
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = combo$se
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (should_run_asserts()) {
				if (!is.finite(private$cached_values$s_beta_hat_T)){
					stop("Rank-regression compound estimator: could not compute a finite standard error.")
				}
			}
		},

		rfit_for_matched_pairs = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			fit = private$fit_matched_component(
				dy = KKstats$y_matched_diffs,
				dX = KKstats$X_matched_diffs,
				estimate_only = estimate_only
			)
			private$cached_values$beta_T_matched = fit$beta
			private$cached_values$ssq_beta_T_matched = fit$ssq
			if (private$include_covariates() && !is.null(KKstats$X_matched_diffs)){
				private$best_Xmm_colnames_matched = colnames(KKstats$X_matched_diffs)
			}
		},

		rfit_for_reservoir = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			fit = private$fit_reservoir_component(
				y_r = KKstats$y_reservoir,
				w_r = KKstats$w_reservoir,
				X_r = KKstats$X_reservoir,
				estimate_only = estimate_only
			)
			private$cached_values$beta_T_reservoir = fit$beta
			private$cached_values$ssq_beta_T_reservoir = fit$ssq
			if (private$include_covariates() && !is.null(KKstats$X_reservoir)){
				private$best_Xmm_colnames_reservoir = colnames(KKstats$X_reservoir)
				# j_treat is 2 in fit_reservoir_component when covariates are included
				private$best_Xmm_j_treat_reservoir = 2L
			}
		},

		fit_matched_component = function(dy, dX = NULL, estimate_only = FALSE){
			if (length(dy) == 0L){
				return(list(beta = NA_real_, ssq = NA_real_))
			}

			if (private$include_covariates()){
				dX = as.matrix(dX)
				x = cbind(`(Intercept)` = rep(1, length(dy)), dX)
				return(private$fit_rank_regression_exact(dy, x, 1L))
			}

			# Univariate case: Rfit fails on dy ~ 1. Use Wilcoxon HL estimate.
			if (all(dy == 0)){
				return(list(beta = 0, ssq = NA_real_))
			}
			
			if (estimate_only) {
				return(list(beta = wilcox_hl_signed_rank_point_estimate_cpp(as.numeric(dy)), ssq = NA_real_))
			}

			mod = tryCatch(
				stats::wilcox.test(dy, conf.int = TRUE),
				error = function(e) NULL
			)

			if (is.null(mod)){
				return(list(beta = NA_real_, ssq = NA_real_))
			}

			beta = as.numeric(mod$estimate)
			se = (as.numeric(mod$conf.int[2]) - as.numeric(mod$conf.int[1])) / (2 * 1.96)
			list(
				beta = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_,
				ssq = if (length(se) == 1L && is.finite(se) && se > 0) se^2 else NA_real_
			)
		},

		fit_reservoir_component = function(y_r, w_r, X_r = NULL, estimate_only = FALSE){
			if (length(y_r) == 0L || !any(w_r == 1) || !any(w_r == 0)){
				return(list(beta = NA_real_, ssq = NA_real_))
			}

			if (private$include_covariates()){
				X_r = as.matrix(X_r)
				x = cbind(`(Intercept)` = rep(1, length(y_r)), w = w_r, X_r)
				return(private$fit_rank_regression_exact(y_r, x, 2L))
			}

			if (estimate_only) {
				return(list(beta = wilcox_hl_point_estimate_cpp(as.numeric(y_r), as.integer(w_r)), ssq = NA_real_))
			}

			yT = y_r[w_r == 1]
			yC = y_r[w_r == 0]
			mod = tryCatch(
				stats::wilcox.test(yT, yC, conf.int = TRUE),
				error = function(e) NULL
			)

			if (is.null(mod)){
				return(list(beta = NA_real_, ssq = NA_real_))
			}

			beta = as.numeric(mod$estimate)
			se = (as.numeric(mod$conf.int[2]) - as.numeric(mod$conf.int[1])) / (2 * 1.96)
			list(
				beta = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_,
				ssq = if (length(se) == 1L && is.finite(se) && se > 0) se^2 else NA_real_
			)
		},

		compute_fast_bootstrap_distr = function(B, i_reservoir, n_reservoir, m, y, w, m_vec){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}

			matched_dy = KKstats$y_matched_diffs
			matched_dX = if (private$include_covariates()) as.matrix(KKstats$X_matched_diffs) else NULL
			y_res = KKstats$y_reservoir
			w_res = KKstats$w_reservoir
			X_res = if (private$include_covariates()) as.matrix(KKstats$X_reservoir) else NULL

			# Capture original SSQ values for variance-weighting.
			# fit_matched_component / fit_reservoir_component always return ssq = NA_real_
			# under estimate_only = TRUE, so we fix the weights at the original-data values.
			orig_ssq_m = private$cached_values$ssq_beta_T_matched
			orig_ssq_r = private$cached_values$ssq_beta_T_reservoir

			compute_one_bootstrap_estimate = function(b){
				reservoir_fit = list(beta = NA_real_, ssq = NA_real_)
				if (n_reservoir > 0){
					res_idx = sample.int(n_reservoir, n_reservoir, replace = TRUE)
					reservoir_fit = private$fit_reservoir_component(
						y_r = y_res[res_idx],
						w_r = w_res[res_idx],
						X_r = if (private$include_covariates()) X_res[res_idx, , drop = FALSE] else NULL,
						estimate_only = TRUE
					)
				}

				matched_fit = list(beta = NA_real_, ssq = NA_real_)
				if (m > 0){
					pair_idx = sample.int(m, m, replace = TRUE)
					matched_fit = private$fit_matched_component(
						dy = matched_dy[pair_idx],
						dX = if (private$include_covariates()) matched_dX[pair_idx, , drop = FALSE] else NULL,
						estimate_only = TRUE
					)
				}

				private$combine_component_estimates(
					beta_m = matched_fit$beta,
					ssq_m = orig_ssq_m,
					beta_r = reservoir_fit$beta,
					ssq_r = orig_ssq_r
				)$beta
			}

			if (self$num_cores == 1L){
				return(vapply(seq_len(B), compute_one_bootstrap_estimate, numeric(1)))
			}

			unlist(private$par_lapply(seq_len(B), function(b){
				compute_one_bootstrap_estimate(b)
			}, n_cores = self$num_cores, show_progress = private$verbose))
		}
	)
)
