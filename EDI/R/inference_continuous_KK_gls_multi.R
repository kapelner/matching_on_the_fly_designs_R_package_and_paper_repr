#' Generalized Least Squares Inference based on Maximum Likelihood for KK designs
#'
#' The methods that support confidence intervals and testing for the mean difference
#' in continuous response types for sequential experimental designs using Generalized Least
#' Squares (GLS).
#' This model assumes a shared correlation within matched pairs but independence for subjects in
#' the reservoir.
#' When no matched pairs exist (pure reservoir design), the model degenerates to OLS.
#'
#' @details
#' The GLS model is fit by REML, which provides unbiased estimates of the within-pair
#' correlation
#' \eqn{\rho} and the residual variance \eqn{\sigma^2}. The residual degrees of freedom for
#' t-based inference on the treatment effect are set to \eqn{n - p - 1}, where \eqn{p} is the
#' number of fixed-effect parameters (intercept + treatment + covariates) and the extra \eqn{-1}
#' accounts for the estimated \eqn{\rho}. Treating \eqn{\rho} as known would yield \eqn{n - p},
#' which slightly overstates precision; the conservative \eqn{n - p - 1} correction is used
#' instead.
#' When the model falls back to OLS (no matched pairs), the standard OLS residual df \eqn{n - p}
#' is used since no correlation parameter is estimated.
#'
#' This class requires the \pkg{nlme} package, which is listed in Suggests and is not
#' installed automatically with \pkg{EDI}.
#' Install \pkg{nlme} before using this class.
#'
#' @export
InferenceContinMultGLS = R6::R6Class("InferenceContinMultGLS",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		#' @examples
		#' set.seed(1)
		#' x_dat <- data.frame(
		#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
		#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
		#' )
		#' seq_des <- DesignSeqOneByOneKK14$new(n = nrow(x_dat), response_type = "continuous", verbose =
		#' FALSE)
		#' for (i in seq_len(nrow(x_dat))) {
		#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
		#' }
		#' seq_des$add_all_subject_responses(c(1.2, 0.9, 1.5, 1.8, 2.1, 1.7, 2.6, 2.2))
		#' infer <- InferenceContinMultGLS$new(seq_des, verbose = FALSE)
		#' infer
		#'
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "continuous")
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)

			if (!requireNamespace("nlme", quietly = TRUE)) {
				stop("Package 'nlme' is required for InferenceContinMultGLS. Please install it.")
			}
		},

		#' @description
		#' Computes the appropriate GLS estimate
		#'
		#' @return	The setting-appropriate numeric estimate of the treatment effect
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			if (is.null(private$cached_values$beta_hat_T)){
				private$shared(estimate_only = estimate_only)
			}
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval
		#' differently for all response types, estimate types, and
		#' test types.
		#'
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		#'
		#' @return	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a 2-sided p-value
		#'
		#' @param delta The null difference to test against. For any
		#'   treatment effect at all this is set to zero (the default).
		#'
		#' @return	The approximate frequentist p-value
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$df)){
				private$shared()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		fit_warm_rho = NULL,

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},
		fit_group_id_signature = NULL,
		fit_group_id_cached = NULL,

		get_ci_fit_controls = function(){
			ctrl = private$randomization_mc_control
			list(
				warm_start = !is.null(ctrl) && isTRUE(ctrl$fit_warm_start_enable),
				reuse_factorizations = !is.null(ctrl) && isTRUE(ctrl$fit_reuse_factorizations)
			)
		},

		set_failed_fit_cache = function(){
			private$cache_nonestimable_estimate("kk_gls_fit_unavailable")
			private$cached_values$is_z = FALSE
		},

		get_group_id = function(m_vec, reuse_factorizations = FALSE){
			signature = paste(m_vec, collapse = ",")
			if (reuse_factorizations && identical(signature, private$fit_group_id_signature) && !is.null(private$fit_group_id_cached)) {
				return(private$fit_group_id_cached)
			}
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L) {
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			}
			if (reuse_factorizations) {
				private$fit_group_id_signature = signature
				private$fit_group_id_cached = group_id
			}
			group_id
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			fit_controls = private$get_ci_fit_controls()
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			full_X_matrix = private$create_design_matrix()
			X_model = full_X_matrix[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			dat = data.frame(y = private$y, X_model)

			m = max(m_vec)
			if (m == 0L) {
				mod_ols = stats::lm(y ~ ., data = dat)
				coef_vec = tryCatch(stats::coef(mod_ols), error = function(e) NULL)
				beta = if (!is.null(coef_vec) && "w" %in% names(coef_vec)) as.numeric(coef_vec[["w"]]) else NA_real_
				private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
				if (estimate_only) return(invisible(NULL))
				cs = stats::coef(summary(mod_ols))
				private$cached_values$s_beta_hat_T = cs["w", "Std. Error"]
				private$cached_values$df = mod_ols$df.residual
				private$cached_values$is_z = FALSE
				return(invisible(NULL))
			}

			dat$group_id = factor(private$get_group_id(m_vec, reuse_factorizations = fit_controls$reuse_factorizations))
			cor_struct = if (fit_controls$warm_start && estimate_only && is.finite(private$fit_warm_rho)) {
				nlme::corCompSymm(value = max(min(private$fit_warm_rho, 0.95), -0.95), form = ~ 1 | group_id)
			} else {
				nlme::corCompSymm(form = ~ 1 | group_id)
			}

			mod = tryCatch({
				nlme::gls(y ~ . - group_id, data = dat,
					correlation = cor_struct,
					method = "REML")
			}, error = function(e) NULL)

			if (is.null(mod)) {
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			coef_vec = tryCatch(stats::coef(mod), error = function(e) NULL)
			beta = if (!is.null(coef_vec) && "w" %in% names(coef_vec)) as.numeric(coef_vec[["w"]]) else NA_real_
			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			rho = tryCatch(as.numeric(stats::coef(mod$modelStruct$corStruct, unconstrained = FALSE)[1]), error = function(e) NA_real_)
			if (is.finite(rho)) private$fit_warm_rho = rho
			if (estimate_only) return(invisible(NULL))

			coef_table = summary(mod)$tTable
			treat_idx = which(rownames(coef_table) == "w")
			if (length(treat_idx) == 0L) treat_idx = 2L
			private$cached_values$s_beta_hat_T = coef_table[treat_idx, "Std.Error"]
			private$cached_values$df = mod$dims$N - mod$dims$p - 1L
			private$cached_values$is_z = FALSE
		}
	)
)
