#' Univariate Quasi-Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a quasi-Poisson log-link regression for count responses using only the
#' treatment indicator. The treatment effect is reported on the log-rate scale and
#' inference uses the model-based quasi-Poisson variance with an estimated
#' dispersion parameter.
#'
#' @export
SeqDesignInferenceCountUnivQuasiPoissonRegr = R6::R6Class("SeqDesignInferenceCountUnivQuasiPoissonRegr",
	inherit = SeqDesignInferenceCountUnivRobustPoissonRegr,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param	seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param	num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs),
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param	verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		fit_count_model_with_var = function(Xmm){
			reduced = private$reduce_design_matrix_preserving_treatment(Xmm)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			mod = tryCatch(
				fast_quasipoisson_regression_with_var_cpp(X_fit, private$y, j = reduced$j_treat),
				error = function(e) NULL
			)
			if (is.null(mod) || !isTRUE(mod$converged)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			coef_hat = as.numeric(mod$b)
			ssq_b_2 = as.numeric(mod$ssq_b_j)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat)) || !is.finite(ssq_b_2) || ssq_b_2 < 0){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			b_full = rep(NA_real_, ncol(Xmm))
			b_full[reduced$keep] = coef_hat
			names(b_full) = colnames(Xmm)
			list(b = b_full, ssq_b_2 = ssq_b_2)
		}
	)
)
