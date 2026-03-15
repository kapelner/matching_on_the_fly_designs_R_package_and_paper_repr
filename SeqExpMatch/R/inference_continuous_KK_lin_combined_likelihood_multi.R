#' Lin (2013) Combined-Likelihood Inference for KK Designs with Continuous Responses
#'
#' @description
#' Fits a stacked Lin (2013) style covariate-adjusted estimator for KK
#' matching-on-the-fly designs with continuous responses using both matched pairs
#' and reservoir observations in one combined fit.
#'
#' @details
#' The class name uses `CombinedLikelihood` to match the package terminology used
#' for other KK stacked estimators. Statistically, this implementation uses a
#' single stacked OLS fit with HC2 heteroskedasticity-robust standard errors,
#' not a literal homoscedastic Gaussian likelihood.
#'
#' The stacked parameterization uses:
#' \itemize{
#'   \item reservoir rows with columns \eqn{(\alpha, \tau, \eta, \gamma)}
#'         corresponding to \eqn{1}, treatment, centered covariates, and
#'         \eqn{(w - 1/2) x_c};
#'   \item matched-pair rows with columns \eqn{(0, 1, \Delta x, \bar x_c)},
#'         scaled by \eqn{1 / \sqrt{2}} to parallel the existing KK combined OLS
#'         weighting.
#' }
#'
#' @export
SeqDesignInferenceContinMultiKKLinCombinedLikelihood = R6::R6Class("SeqDesignInferenceContinMultiKKLinCombinedLikelihood",
	inherit = SeqDesignInferenceKKPassThroughCompound,
	public = list(

		#' @description
		#' Initialize the inference object.
		#'
		#' @param seq_des_obj A completed KK \code{SeqDesign} object with a continuous
		#'   response.
		#' @param num_cores The number of CPU cores to use for bootstrap and
		#'   randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "continuous")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the stacked combined estimate of the treatment effect.
		#'
		#' @return The estimated treatment effect.
		compute_treatment_estimate = function(){
			private$fit_combined()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval.
		#'
		#' @param alpha The confidence level in the computed confidence interval is
		#'   1 - \code{alpha}. The default is 0.05.
		#'
		#' @return A confidence interval for the treatment effect.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$fit_combined()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a two-sided p-value.
		#'
		#' @param delta The null treatment effect. Defaults to 0.
		#'
		#' @return The approximate p-value.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$fit_combined()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			private$cached_values$KKstats = .compute_kk_lin_basic_match_data(
				X = private$X,
				n = private$n,
				y = private$y,
				w = private$w,
				match_indic = private$match_indic
			)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				stop("KK Lin combined estimator: could not compute a finite HC2 standard error.")
			}
		},

		reduce_design_matrix_preserving_required = function(X_full, required){
			qr_X = qr(X_full)
			target_rank = qr_X$rank
			candidate_order = c(required, setdiff(qr_X$pivot, required))
			keep = integer(0)

			for (j in candidate_order){
				trial_keep = c(keep, j)
				trial_rank = qr(X_full[, trial_keep, drop = FALSE])$rank
				if (trial_rank > length(keep)){
					keep = trial_keep
				}
				if (length(keep) >= target_rank){
					break
				}
			}

			sort(unique(keep))
		},

		fit_combined = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}

			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			nR = nRT + nRC
			p = ncol(private$get_X())
			xbar = colMeans(private$get_X())
			scale_pair = 1 / sqrt(2)

			pair_rows =
				if (m > 0){
					Xd = as.matrix(KKstats$X_matched_diffs_full)
					Xm_centered = sweep(as.matrix(KKstats$X_matched_means_full), 2, xbar, "-")
					cbind(
						0,
						rep(scale_pair, m),
						scale_pair * Xd,
						scale_pair * Xm_centered
					)
				} else {
					matrix(0, nrow = 0, ncol = 2 + 2 * p)
				}

			reservoir_rows =
				if (nR > 0){
					Xr = as.matrix(KKstats$X_reservoir)
					Xc = sweep(Xr, 2, xbar, "-")
					cbind(
						1,
						KKstats$w_reservoir,
						Xc,
						(KKstats$w_reservoir - 0.5) * Xc
					)
				} else {
					matrix(0, nrow = 0, ncol = 2 + 2 * p)
				}

			X_full = rbind(pair_rows, reservoir_rows)
			y_full = c(
				if (m > 0) scale_pair * KKstats$y_matched_diffs else numeric(0),
				if (nR > 0) KKstats$y_reservoir else numeric(0)
			)
			j_target = 2L

			if (nrow(X_full) == 0){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			keep = private$reduce_design_matrix_preserving_required(X_full, required = c(1L, 2L))
			if (!(j_target %in% keep)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			X_fit = X_full[, keep, drop = FALSE]
			j_fit = match(j_target, keep)
			if (nrow(X_fit) <= ncol(X_fit)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			mod = stats::lm.fit(X_fit, y_full)
			coef_hat = as.numeric(mod$coefficients)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			post_fit = tryCatch(
				ols_hc2_post_fit_cpp(X_fit, y_full, coef_hat, j_fit),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			beta_hat = post_fit$beta_hat
			vcov_hc2 = post_fit$vcov
			ssq_hat = post_fit$ssq_hat
			if (!is.finite(beta_hat) || !is.finite(ssq_hat) || ssq_hat <= 0){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			coef_names = c(
				"(Intercept)",
				"treatment",
				paste0("eta:", colnames(private$get_X())),
				paste0("gamma:", colnames(private$get_X()))
			)
			if (length(coef_names) != ncol(X_full)){
				coef_names = paste0("V", seq_len(ncol(X_full)))
				coef_names[1:2] = c("(Intercept)", "treatment")
			}
			coef_names = coef_names[keep]
			names(coef_hat) = coef_names
			colnames(vcov_hc2) = rownames(vcov_hc2) = coef_names

			private$cached_values$beta_hat_T = beta_hat
			private$cached_values$s_beta_hat_T = sqrt(ssq_hat)
			private$cached_values$is_z = TRUE
			private$cached_values$full_coefficients = coef_hat
			private$cached_values$full_vcov = vcov_hc2
		}
	)
)
