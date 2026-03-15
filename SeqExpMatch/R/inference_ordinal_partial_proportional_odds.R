#' Partial Proportional-Odds Inference for Ordinal Responses
#'
#' @description
#' Fits a partial proportional-odds model for ordinal responses. Treatment remains a
#' parallel predictor while user-specified covariates can vary across thresholds.
#' This class works on any completed sequential design and reuses the fast Rcpp
#' logistic solver after expanding the cumulative logits.
#'
#' @param seq_des_obj A completed \code{SeqDesign} object whose response_type is "ordinal".
#' @param nonparallel Character vector of predictor names (from \code{SeqDesign$get_X()})
#'   that should have threshold-specific slopes. Defaults to \code{character(0)}.
#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
#' @param verbose Whether to print progress messaging.
#'
#' @export
SeqDesignInferenceOrdinalPartialProportionalOdds = R6::R6Class(
	"SeqDesignInferenceOrdinalPartialProportionalOdds",
	inherit = SeqDesignInference,
	public = list(
		#' @description
		#' Initialize the partial proportional-odds inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an ordinal response.
		#' @param nonparallel Predictor names that should vary across thresholds.
		#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
		#' @param verbose Whether to emit progress messaging.
		initialize = function(seq_des_obj,
				nonparallel = character(0),
				num_cores = 1,
				verbose = FALSE) {
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			super$initialize(seq_des_obj, num_cores, verbose)
			private$nonparallel = unique(nonparallel)
			available = colnames(private$get_X())
			private$nonparallel = intersect(private$nonparallel, available)
			if ("treatment" %in% private$nonparallel) {
				private$nonparallel = setdiff(private$nonparallel, "treatment")
			}
		},

		#' @description
		#' Retrieves the estimated log-odds shift for the treatment arm.
		compute_treatment_estimate = function() {
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a Wald-style CI using the beta concentration estimate.
		#' @param alpha Significance level for the interval.
		compute_mle_confidence_interval = function(alpha = 0.05) {
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a Wald p-value for the (parallel) treatment effect.
		#' @param delta Null treatment effect to test.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0) {
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Randomization-based CIs are not supported for this estimator.
		#' @param ... Unused.
		compute_confidence_interval_rand = function(...) {
			stop("Randomization confidence intervals are not supported for the partial proportional-odds estimator.")
		}
	),

	private = list(
		nonparallel = character(0),

		shared = function() {
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			y = as.numeric(private$y)
			levels = sort(unique(y))
			K = length(levels)
			if (K < 2) {
				private$set_failed_cache()
				return(invisible(NULL))
			}

			y_int = match(y, levels)
			parallel_covars = private$get_parallels()
			nonparallel_covars = private$get_nonparallels()

			stacked = private$build_stacked_design(y_int, parallel_covars, nonparallel_covars)
			if (is.null(stacked)) {
				private$set_failed_cache()
				return(invisible(NULL))
			}

			mod = tryCatch(
				fast_logistic_regression_with_var_cpp(stacked$X_design, stacked$y_stack, stacked$j_treat),
				error = function(e) NULL
			)
			if (is.null(mod)) {
				private$set_failed_cache()
				return(invisible(NULL))
			}

			coef = as.numeric(mod$b)
			ssq = as.numeric(mod$ssq_b_2)
			if (length(coef) < stacked$treat_idx || !is.finite(coef[stacked$treat_idx]) || !is.finite(ssq) || ssq <= 0) {
				private$set_failed_cache()
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = coef[stacked$treat_idx]
			private$cached_values$s_beta_hat_T = sqrt(ssq)
			private$cached_values$is_z = TRUE
		},

		build_stacked_design = function(y_int, X_par_additional, X_non) {
			K = length(unique(y_int))
			n = length(y_int)
			nrep = n * (K - 1)
			if (nrep <= 0) return(NULL)

			intercepts = matrix(0, nrow = nrep, ncol = K - 1)
			y_stack = integer(nrep)
			for (j in seq_len(K - 1)) {
				rows = seq(j, by = K - 1, length.out = n)
				intercepts[rows, j] = 1
				y_stack[rows] = as.integer(y_int <= j)
			}

			X_par = cbind(treatment = private$w)
			if (ncol(X_par_additional) > 0) {
				X_par = cbind(X_par, as.matrix(X_par_additional))
			}

			X_par_stack = X_par[rep(seq_len(n), each = K - 1), , drop = FALSE]
			intercept_cols = ncol(intercepts)

			if (ncol(X_non) == 0) {
				X_non_stack = matrix(numeric(0), nrow = nrep, ncol = 0)
			} else {
				m = ncol(X_non)
				X_non_stack = matrix(0, nrow = nrep, ncol = (K - 1) * m)
				for (j in seq_len(K - 1)) {
					rows = seq(j, by = K - 1, length.out = n)
					cols = ((j - 1) * m + seq_len(m))
					X_non_stack[rows, cols] = X_non
				}
			}

			X_design = cbind(intercepts, X_par_stack, X_non_stack)
			list(
				X_design = X_design,
				y_stack = as.numeric(y_stack),
				j_treat = intercept_cols + 1,
				treat_idx = intercept_cols + 1
			)
		},

		get_parallels = function() {
			X = private$get_X()
			if (ncol(X) == 0) return(matrix(0, nrow = private$n, ncol = 0))
			parallel_names = setdiff(colnames(X), private$nonparallel)
			if (length(parallel_names) == 0) return(matrix(0, nrow = private$n, ncol = 0))
			as.matrix(X[, parallel_names, drop = FALSE])
		},

		get_nonparallels = function() {
			X = private$get_X()
			cand = intersect(private$nonparallel, colnames(X))
			if (length(cand) == 0) return(matrix(0, nrow = private$n, ncol = 0))
			as.matrix(X[, cand, drop = FALSE])
		},

		assert_finite_se = function() {
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0) {
				stop("Partial proportional-odds estimator: could not compute a finite standard error.")
			}
		},

		set_failed_cache = function() {
			private$cached_values$beta_hat_T = NA_real_
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$is_z = TRUE
		}
	)
)
