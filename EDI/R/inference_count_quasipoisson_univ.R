#' Univariate Quasi-Poisson Regression Inference for Count Responses
#'
#' Fits a quasi-Poisson log-link regression for count responses using only the
#' treatment indicator. The treatment effect is reported on the log-rate scale and
#' inference uses the model-based quasi-Poisson variance with an estimated
#' dispersion parameter.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "count",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- InferenceCountUnivQuasiPoissonRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceCountUnivQuasiPoissonRegr = R6::R6Class("InferenceCountUnivQuasiPoissonRegr",
	lock_objects = FALSE,
	inherit = InferenceCountUnivRobustPoissonRegr,
	public = list(

	),

	private = list(
		fit_count_model_with_var = function(Xmm, estimate_only = FALSE){
			reduced = private$reduce_design_matrix_preserving_treatment_fixed_covariates(Xmm)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			mod = tryCatch(
				if (estimate_only) {
					fast_poisson_regression_cpp(X_fit, private$y)
				} else {
					fast_quasipoisson_regression_with_var_cpp(X_fit, private$y, j = reduced$j_treat)
				},
				error = function(e) NULL
			)
			if (is.null(mod) || !isTRUE(mod$converged)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			coef_hat = as.numeric(mod$b)
			ssq_b_2 = if (estimate_only) NA_real_ else as.numeric(mod$ssq_b_j)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}
			if (!estimate_only && (!is.finite(ssq_b_2) || ssq_b_2 < 0)){
				ssq_b_2 = NA_real_
			}

			b_full = rep(NA_real_, ncol(Xmm))
			b_full[reduced$keep] = coef_hat
			names(b_full) = colnames(Xmm)
			list(b = b_full, ssq_b_2 = ssq_b_2)
		}
	)
)
