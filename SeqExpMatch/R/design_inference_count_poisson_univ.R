#' Univariate Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a Poisson log-link regression for count responses using only the treatment
#' indicator. The treatment effect is reported on the log-rate scale and inference
#' uses the model-based Poisson variance.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "count",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- DesignInferenceCountUnivPoissonRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceCountUnivPoissonRegr = R6::R6Class("DesignInferenceCountUnivPoissonRegr",
	inherit = DesignInferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object
		#' after the sequential design is completed.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects
		#'   are assigned and response y is recorded within.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param verbose A flag indicating whether messages should be
		#'   displayed to the user. Default is \code{TRUE}.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "count")
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the appropriate estimate.
		#'
		#' @return	The log-rate treatment-effect estimate.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		fit_poisson_with_var = function(Xmm){
			reduced = private$reduce_design_matrix_preserving_treatment(Xmm)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			mod = tryCatch(
				fast_poisson_regression_with_var_cpp(X_fit, private$y, j = reduced$j_treat),
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
		},

		generate_mod = function(){
			Xmm = cbind(1, private$w)
			colnames(Xmm) = c("(Intercept)", "treatment")
			private$fit_poisson_with_var(Xmm)
		}
	)
)
