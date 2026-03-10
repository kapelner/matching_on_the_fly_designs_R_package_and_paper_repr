#' OLS Compound Combined-Likelihood Inference for KK Designs
#'
#' @description
#' Fits a single joint normal likelihood over all KK design data — matched-pair
#' differences and reservoir subjects — estimating the treatment effect \eqn{\beta_T}
#' in one combined OLS.
#'
#' \strong{Model:} \eqn{y = \beta_0 + \beta_T w + \mathbf{x}'\boldsymbol{\beta_{xs}} + \varepsilon},
#' \eqn{\varepsilon \sim N(0, \sigma^2)}. Differencing within pairs eliminates the
#' pair-specific intercept and yields
#' \eqn{d_k = \beta_T + \mathbf{x}_{d,k}'\boldsymbol{\beta_{xs}} + \varepsilon_{d,k}},
#' \eqn{\varepsilon_{d,k} \sim N(0, 2\sigma^2)}.
#' The same \eqn{\boldsymbol{\beta_{xs}}} appears in both components.
#'
#' \strong{Combined design matrix} (column layout: \eqn{\beta_0,\,\beta_T,\,\boldsymbol{\beta_{xs}}}):
#' \itemize{
#'   \item Matched-pair rows (scaled by \eqn{1/\sqrt{2}}):
#'         \eqn{[0,\; 1/\sqrt{2},\; \mathbf{x}_{d,k}'/\sqrt{2}]}.
#'   \item Reservoir rows:
#'         \eqn{[1,\; w_i,\; \mathbf{x}_{r,i}']}.
#' }
#' Scaling pair rows by \eqn{1/\sqrt{2}} equalises all residual variances to
#' \eqn{\sigma^2}, giving a standard homoscedastic OLS. The reservoir intercept
#' \eqn{\beta_0} is a nuisance parameter (zero for pair rows).
#'
#' The Hessian of the log-likelihood with respect to
#' \eqn{(\beta_0, \beta_T, \boldsymbol{\beta_{xs}})} is
#' \eqn{H = -\mathbf{X}'\mathbf{X}/\hat\sigma^2}, stored in
#' \code{private$cached_values$hessian} for downstream use (CI, p-value).
#'
#' @export
SeqDesignInferenceContinMultOLSKKCombinedLikelihood = R6::R6Class("SeqDesignInferenceContinMultOLSKKCombinedLikelihood",
	inherit = SeqDesignInferenceKKPassThroughCompound,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "continuous")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the combined-likelihood OLS estimate of the treatment effect.
		#'
		#' Fits a single joint OLS over scaled matched-pair differences and reservoir
		#' observations with shared covariate effects. The Hessian
		#' \eqn{H = -\mathbf{X}'\mathbf{X}/\hat\sigma^2} is stored in
		#' \code{private$cached_values$hessian} as a side-effect.
		#'
		#' @return	Numeric scalar \eqn{\hat\beta_T}.
		compute_treatment_estimate = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				private$fit_combined_likelihood()
			}
			private$cached_values$beta_hat_T
		},

		compute_mle_confidence_interval = function(alpha = 0.05)
			stop(paste(class(self)[1], ": combined-likelihood confidence interval not yet implemented.")),
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0)
			stop(paste(class(self)[1], ": combined-likelihood p-value not yet implemented."))
	),

	private = list(

		# Fit the combined normal likelihood over matched-pair differences (scaled) and
		# reservoir observations with SHARED covariate effects beta_xs.
		#
		# Column layout of X_comb: [beta_0 | beta_T | beta_xs (p cols)]
		#   col 1       : beta_0   (0 for pair rows; 1 for reservoir rows)
		#   col 2       : beta_T   (1/sqrt2 for pair rows; w_i for reservoir rows)
		#   cols 3..p+2 : beta_xs  (Xd/sqrt2 for pair rows; X_r for reservoir rows)
		#
		# The combined case uses build_kk_combined_ols_design_cpp to construct
		# X_comb and y_comb in a single C++ pass (no intermediate R allocations).
		# Pair rows are scaled by 1/sqrt(2) so that all residuals are N(0, sigma^2).
		# When only one component is present the design reduces accordingly (no beta_0
		# column for pairs-only; same layout minus the zero-padding for reservoir-only).
		fit_combined_likelihood = function(){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}

			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			nR  = nRT + nRC

			# ---- Build design matrix and response --------------------------------
			if (m > 0 && nRT > 0 && nRC > 0){
				# Both components: shared beta_xs, layout [beta_0 | beta_T | beta_xs]
				# Built in one C++ pass — no intermediate R matrix allocations.
				design   = build_kk_combined_ols_design_cpp(
					KKstats$y_matched_diffs,
					as.matrix(KKstats$X_matched_diffs),
					KKstats$y_reservoir,
					KKstats$w_reservoir,
					as.matrix(KKstats$X_reservoir)
				)
				X_comb   = design$X_comb
				y_comb   = design$y_comb
				j_beta_T = 2L

			} else if (m > 0){
				# Matched pairs only: d_k = beta_T + Xd_k' beta_xs + eps_k
				# Layout: [beta_T | beta_xs]  (no beta_0 needed)
				X_comb   = cbind(1, as.matrix(KKstats$X_matched_diffs))
				y_comb   = KKstats$y_matched_diffs
				j_beta_T = 1L

			} else if (nRT > 0 && nRC > 0){
				# Reservoir only: layout [beta_0 | beta_T | beta_xs]
				X_comb   = cbind(rep(1, nR), KKstats$w_reservoir, as.matrix(KKstats$X_reservoir))
				y_comb   = KKstats$y_reservoir
				j_beta_T = 2L

			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- QR-reduce to full rank, always preserving beta_T column --------
			qr_X = qr(X_comb)
			if (qr_X$rank < ncol(X_comb)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_beta_T %in% keep)){
					keep[qr_X$rank] = j_beta_T
				}
				keep     = sort(keep)
				X_comb   = X_comb[, keep, drop = FALSE]
				j_beta_T = which(keep == j_beta_T)
			}

			n_total  = nrow(X_comb)
			n_params = ncol(X_comb)

			if (n_total <= n_params){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- Fit OLS --------------------------------------------------------
			mod = tryCatch(
				fast_ols_with_var_cpp(X_comb, y_comb, j = j_beta_T),
				error = function(e) NULL
			)
			if (is.null(mod) || !is.finite(mod$b[j_beta_T])){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- Hessian of log-likelihood w.r.t. (beta_0, beta_T, beta_xs) -----
			# H = -X'X / sigma^2  (negative definite at MLE)
			# mod$XtX and mod$sigma2_hat are already computed inside fast_ols_with_var_cpp;
			# mod$ssq_b_j already holds sigma2_hat * (XtX^{-1})_{jj}.
			private$cached_values$hessian      = -as.matrix(mod$XtX) / mod$sigma2_hat
			private$cached_values$j_beta_T     = j_beta_T
			private$cached_values$beta_hat_T   = as.numeric(mod$b[j_beta_T])
			private$cached_values$s_beta_hat_T = sqrt(mod$ssq_b_j)
			private$cached_values$is_z         = TRUE
			invisible(NULL)
		}
	)
)
