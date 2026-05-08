#' Abstract class for ordinal CLMM-based Inference in KK designs
#'
#' @keywords internal
InferenceAbstractKKOrdinalCLMM = R6::R6Class("InferenceAbstractKKOrdinalCLMM",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula   Optional formula for covariate adjustment.
		#' @param use_rcpp Logical. If \code{TRUE} (default), use the internal Rcpp
		#'   implementation (no external packages required). Set \code{FALSE} to fall
		#'   back to \pkg{ordinal::clmm}.
		#' @param verbose A flag indicating whether messages should be displayed.
		#' @param harden Whether to apply robustness measures.
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE, harden = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "ordinal")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			if (!use_rcpp && should_run_asserts()) {
				if (!check_package_installed("ordinal")){
					stop("Package 'ordinal' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, verbose = verbose, harden = harden, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			private$use_rcpp = use_rcpp
		},

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha The significance level (default 0.05).
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
		#' Compute asymp two sided pval for treatment effect
		#' @param delta The null treatment effect (default 0).
		compute_asymp_two_sided_pval = function(delta = 0){
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
					stop("TO-DO")
				}
				NA_real_
			}
		}
	),

	private = list(
		use_rcpp = TRUE,
		max_abs_reasonable_coef = 1e4,
		best_X_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (private$use_rcpp) {
				return(private$compute_ri_estimate_rcpp())
			}
			# Ensure we have the best design from the original data
			if (is.null(private$best_X_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_X_colnames)){
				return(self$compute_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			X_cols = private$best_X_colnames
			X_data = private$get_X()

			X_fit = if (length(X_cols) == 0L){
				matrix(private$w, ncol = 1, dimnames = list(NULL, "treatment"))
			} else {
				X_cov = X_data[, intersect(X_cols, colnames(X_data)), drop = FALSE]
				cbind(treatment = private$w, X_cov)
			}
			# Add intercept column for clmm internal expectations
			X_fit_full = cbind("(Intercept)" = 1, X_fit)

			mod = private$fit_clmm(X_fit_full)
			if (is.null(mod)) mod = private$fit_clm_fallback(X_fit_full)
			if (is.null(mod)) return(NA_real_)

			as.numeric(stats::coef(mod)["w"])
		},

		compute_ri_estimate_rcpp = function(){
			X_fit = private$clmm_X_for_rcpp()
			group_id = private$clmm_group_id()
			y_levels = sort(unique(private$y))
			K = length(y_levels)
			y = as.integer(match(private$y, y_levels))
			fit = tryCatch(
				fast_ordinal_clmm_cpp(
					X          = X_fit,
					y          = y,
					group_id   = as.integer(group_id),
					K          = K,
					j_T        = 0L,
					link       = private$clmm_link(),
					estimate_only = TRUE,
					eps_g      = 1e-3
				),
				error = function(e) NULL
			)
			if (is.null(fit) || !isTRUE(fit$converged)) return(NA_real_)
			as.numeric(fit$b[1L])
		},

		clmm_link = function() stop(class(self)[1], " must implement clmm_link()"),

		clmm_predictors_df = function(){
			full_X = private$create_design_matrix()
			private$clmm_predictors_df_from_design(full_X)
		},

		clmm_predictors_df_from_design = function(full_X){
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},

		# Build X matrix for Rcpp (no intercept, treatment at col 1)
		clmm_X_for_rcpp = function(){
			as.matrix(private$clmm_predictors_df())
		},

		# Build group_id vector (reservoir singletons get unique IDs)
		clmm_group_id = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)
			group_id
		},

		# Warm start: fixed-effects ordinal MLE with appropriate link
		clmm_warm_start = function(X_fit, y, n_alpha){
			tryCatch({
				warm_fn = switch(private$clmm_link(),
					logit   = fast_ordinal_regression_cpp,
					probit  = fast_ordinal_probit_regression_cpp,
					cauchit = fast_ordinal_cauchit_regression_cpp,
					cloglog = fast_ordinal_cloglog_regression_cpp,
					stop("Unknown link: ", private$clmm_link())
				)
				nore = warm_fn(X_fit, as.numeric(y) - 1L)
				alpha_direct = as.numeric(nore$alpha)
				beta_nore    = as.numeric(nore$b)
				alpha_par = numeric(n_alpha)
				if (n_alpha >= 1L) alpha_par[1L] = alpha_direct[1L]
				if (n_alpha >= 2L) {
					for (k in 2L:n_alpha) {
						diff_k = alpha_direct[k] - alpha_direct[k - 1L]
						alpha_par[k] = if (diff_k > 0) log(diff_k) else 0.0
					}
				}
				c(alpha_par, beta_nore, -3.0)
			}, error = function(e) NULL)
		},

		shared = function(estimate_only = FALSE){
			if (private$use_rcpp) {
				private$shared_rcpp(estimate_only)
			} else {
				private$shared_clmm(estimate_only)
			}
		},

		shared_rcpp = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

			group_id = private$clmm_group_id()
			X_fit    = private$clmm_X_for_rcpp()

			y_levels = sort(unique(private$y))
			K        = length(y_levels)
			y    = as.integer(match(private$y, y_levels))
			n_alpha  = K - 1L
			j_T      = 0L  # treatment is always first column of X_fit

			start = private$clmm_warm_start(X_fit, y, n_alpha)

			fit = tryCatch(
				fast_ordinal_clmm_cpp(
					X             = X_fit,
					y         = y,
					group_id      = as.integer(group_id),
					K             = K,
					j_T           = j_T,
					link          = private$clmm_link(),
					estimate_only = estimate_only,
					start         = start,
					eps_g         = 1e-3
				),
				error = function(e) NULL
			)

			if (is.null(fit) || !isTRUE(fit$converged)) {
				private$cache_nonestimable_estimate("kk_clmm_rcpp_failed")
				return(invisible(NULL))
			}

			beta_hat_T = as.numeric(fit$b[j_T + 1L])

			if (!is.finite(beta_hat_T) || abs(beta_hat_T) > private$max_abs_reasonable_coef) {
				private$cache_nonestimable_estimate("kk_clmm_rcpp_nonestimable")
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$df         = Inf

			if (estimate_only) return(invisible(NULL))

			ssq = fit$ssq_b_T
			private$cached_values$s_beta_hat_T = if (!is.null(ssq) && is.finite(ssq) && ssq > 0) sqrt(ssq) else NA_real_
		},

		shared_clmm = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			full_X = private$create_design_matrix()
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = full_X,
				required_cols = c(1L, 2L),
				fit_fun = function(X_fit){
					mod = private$fit_clmm(X_fit)
					summ = if (!is.null(mod)) tryCatch(summary(mod), error = function(e) NULL) else NULL
					se = if (!is.null(summ)) as.numeric(summ$coefficients["w", "Std. Error"]) else NA_real_
					if (is.null(mod) || (!estimate_only && (!is.finite(se) || se <= 0))){
						mod = private$fit_clm_fallback(X_fit)
						summ = if (!is.null(mod)) tryCatch(summary(mod), error = function(e) NULL) else NULL
						se = if (!is.null(summ)) as.numeric(summ$coefficients["w", "Std. Error"]) else NA_real_
					}
					list(mod = mod, summ = summ, se = se)
				},
				fit_ok = function(fit, X_fit, keep){
					if (is.null(fit) || is.null(fit$mod)) return(FALSE)
					beta = tryCatch(as.numeric(stats::coef(fit$mod)["w"]), error = function(e) NA_real_)
					if (!is.finite(beta)) return(FALSE)
					if (estimate_only) return(TRUE)
					is.finite(fit$se) && fit$se > 0
				}
			)
			mod = attempt$fit$mod
			summ = attempt$fit$summ
			se = attempt$fit$se
			if (!is.null(mod)){
				private$best_X_colnames = setdiff(colnames(attempt$X_fit), c("(Intercept)", "treatment"))
			}
			if (is.null(mod) || is.null(summ)){
				private$cached_values$beta_hat_T   = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(stats::coef(mod)["w"])
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				return(invisible(NULL))
		},

		fit_clmm = function(full_X = private$create_design_matrix()){
			group_id = private$clmm_group_id()

			dat = data.frame(
				y = factor(private$y, ordered = TRUE),
				private$clmm_predictors_df_from_design(full_X),
				group_id = factor(group_id)
			)
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			clmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					ordinal::clmm(
						clmm_formula,
						data = dat,
						link = private$clmm_link()
					)
				)))
				mod
			}, error = function(e) NULL)
		},

		fit_clm_fallback = function(full_X = private$create_design_matrix()){
			dat = data.frame(
				y = factor(private$y, ordered = TRUE),
				private$clmm_predictors_df_from_design(full_X)
			)
			fixed_terms = setdiff(colnames(dat), "y")
			clm_formula = stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					ordinal::clm(
						clm_formula,
						data = dat,
						link = private$clmm_link()
					)
				)))
				mod
			}, error = function(e) NULL)
		}
	)
)

#' Ordinal KK CLMM (Proportional Odds / logit link)
#'
#' Cumulative-logit mixed model for ordinal KK designs.
#' Random intercept per matched pair, treatment + optional covariates as fixed effects.
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKCLMM$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKCLMM = R6::R6Class("InferenceOrdinalKKCLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKOrdinalCLMM,
	public = list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Use internal Rcpp implementation (default \code{TRUE}).
		#' @param verbose Print messages?
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose)
		}
	),
	private = list(
		clmm_link = function() "logit"
	)
)

#' Ordinal KK CLMM (Probit link)
#'
#' Cumulative-probit mixed model for ordinal KK designs.
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKCLMMCauchit$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKCLMMProbit = R6::R6Class("InferenceOrdinalKKCLMMProbit",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKOrdinalCLMM,
	public = list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Use internal Rcpp implementation (default \code{TRUE}).
		#' @param verbose Print messages?
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose)
		}
	)
)

#' Ordinal KK CLMM (Cauchit link)
#'
#' Cumulative-cauchit mixed model for ordinal KK designs.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKCLMMCauchit$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKCLMMCauchit = R6::R6Class("InferenceOrdinalKKCLMMCauchit",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKOrdinalCLMM,
	public = list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Use internal Rcpp implementation (default \code{TRUE}).
		#' @param verbose Print messages?
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose)
		}
	),
	private = list(
		clmm_link = function() "cauchit"
	)
)

#' Ordinal KK CLMM (Complementary log-log link)
#'
#' Cumulative complementary-log-log mixed model for ordinal KK designs.
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'ordinal')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(sample(1:4, 10, replace = TRUE))
#' inf = InferenceOrdinalKKCLMMCloglog$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceOrdinalKKCLMMCloglog = R6::R6Class("InferenceOrdinalKKCLMMCloglog",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKOrdinalCLMM,
	public = list(
		#' @description Initialize
		#' @param des_obj A completed \code{Design} object.
		#' @param model_formula Optional formula for covariate adjustment.
		#' @param use_rcpp Use internal Rcpp implementation (default \code{TRUE}).
		#' @param verbose Print messages?
		initialize = function(des_obj, model_formula = NULL, use_rcpp = TRUE, verbose = FALSE){
			super$initialize(des_obj, model_formula = model_formula, use_rcpp = use_rcpp, verbose = verbose)
		}
	),
	private = list(
		clmm_link = function() "cloglog"
	)
)
