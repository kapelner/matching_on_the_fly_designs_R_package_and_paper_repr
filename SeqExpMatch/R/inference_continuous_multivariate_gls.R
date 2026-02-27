#' Generalized Least Squares Inference based on Maximum Likelihood for KK designs
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in continuous response types for sequential experimental designs using Generalized Least Squares (GLS).
#' This model assumes a shared correlation within matched pairs but independence for subjects in the reservoir.
#' When no matched pairs exist (pure reservoir design), the model degenerates to OLS.
#'
#' @details
#' The GLS model is fit by REML, which provides unbiased estimates of the within-pair correlation
#' \eqn{\rho} and the residual variance \eqn{\sigma^2}. The residual degrees of freedom for
#' t-based inference on the treatment effect are set to \eqn{n - p - 1}, where \eqn{p} is the
#' number of fixed-effect parameters (intercept + treatment + covariates) and the extra \eqn{-1}
#' accounts for the estimated \eqn{\rho}. Treating \eqn{\rho} as known would yield \eqn{n - p},
#' which slightly overstates precision; the conservative \eqn{n - p - 1} correction is used instead.
#' When the model falls back to OLS (no matched pairs), the standard OLS residual df \eqn{n - p}
#' is used since no correlation parameter is estimated.
#'
#' @export
SeqDesignInferenceContinMultGLS = R6::R6Class("SeqDesignInferenceContinMultGLS",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference 
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs), 
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "continuous")			
			super$initialize(seq_des_obj, num_cores, verbose)	
			assertNoCensoring(private$any_censoring)
			
			if (!requireNamespace("nlme", quietly = TRUE)) {
				stop("Package 'nlme' is required for SeqDesignInferenceContinMultGLS. Please install it.")
			}
		},
		
		#' @description
		#' Computes the appropriate GLS estimate
		#' 
		#' @return 	The setting-appropriate numeric estimate of the treatment effect
		compute_treatment_estimate = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				private$shared()
			}
			private$cached_values$beta_hat_T
		},
		
		#' @description		
		#' Computes a 1-alpha level frequentist confidence interval differently for all response types, estimate types and test types.
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' 
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}			
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},		
		
		#' @description
		#' Computes a 2-sided p-value
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$df)){
				private$shared()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	
	private = list(
		shared = function(){
			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			full_X_matrix = private$create_design_matrix()
			# Full matrix has Intercept, w, and covariates.
			# Drop col 1 (manual intercept); gls/lm add their own.
			X_model = full_X_matrix[, -1, drop = FALSE]
			# Pin the treatment column name so coefficient lookup is unambiguous.
			colnames(X_model)[1] = "w"

			dat = data.frame(y = private$y, X_model)

			m = max(match_indic)
			if (m == 0L) {
				# No matched pairs — degenerate GLS equals OLS, so run OLS directly.
				mod_ols = lm(y ~ ., data = dat)
				cs = coef(summary(mod_ols))
				private$cached_values$beta_hat_T   = cs["w", "Estimate"]
				private$cached_values$s_beta_hat_T = cs["w", "Std. Error"]
				private$cached_values$df            = mod_ols$df.residual
				private$cached_values$is_z          = FALSE
				return(invisible(NULL))
			}

			# Matched pairs share a group_id; reservoir subjects each get a unique one.
			group_id = match_indic
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L) {
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			}
			dat$group_id = factor(group_id)

			# Fit GLS with compound symmetry within matched pairs (REML for unbiased variance components).
			mod = tryCatch({
				nlme::gls(y ~ . - group_id, data = dat,
					correlation = nlme::corCompSymm(form = ~ 1 | group_id),
					method = "REML")
			}, error = function(e) NULL)

			if (is.null(mod)) {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$df            = NA_real_
				private$cached_values$is_z          = FALSE
				return(invisible(NULL))
			}

			coef_table = summary(mod)$tTable
			treat_idx  = which(rownames(coef_table) == "w")
			if (length(treat_idx) == 0L) treat_idx = 2L

			private$cached_values$beta_hat_T   = coef_table[treat_idx, "Value"]
			private$cached_values$s_beta_hat_T = coef_table[treat_idx, "Std.Error"]
			# Subtract 1 extra df for the estimated within-pair correlation rho.
			# GLS uses REML to estimate rho and sigma^2; treating rho as known (df = N - p)
			# slightly overstates precision, so we use N - p - 1 as a conservative correction.
			private$cached_values$df            = mod$dims$N - mod$dims$p - 1L
			private$cached_values$is_z          = FALSE
		}
	)
)
