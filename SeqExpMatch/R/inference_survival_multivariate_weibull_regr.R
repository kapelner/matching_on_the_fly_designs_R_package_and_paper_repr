#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#'
#'
#' @export
SeqDesignInferenceSurvivalMultiWeibullRegr = R6::R6Class("SeqDesignInferenceSurvivalMultiWeibullRegr",
	inherit = SeqDesignInferenceSurvivalUniWeibullRegr,
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
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		generate_mod = function(){
			# Multivariate: covariates + treatment (mirrors SeqDesignInferenceSurvivalMultiCoxPHRegr)
			X_cov_orig = private$get_X()
			# Try fast C++ path with progressively lower correlation thresholds.
			# The robust survreg fallback (slow: up to 50 random restarts) is invoked at most
			# ONCE after all fast paths are exhausted, so bootstrap iterations stay fast.
			thresholds = c(Inf, 0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			full_X_matrix_last = NULL
			for (thresh in thresholds) {
				if (ncol(X_cov_orig) > 0) {
					if (is.finite(thresh)) {
						X_cov = drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M
					} else {
						X_cov = X_cov_orig
					}
					full_X_matrix = cbind(X_cov, private$w)
					colnames(full_X_matrix) = c(colnames(X_cov), "treatment")
				} else {
					full_X_matrix = matrix(private$w, ncol = 1)
					colnames(full_X_matrix) = "treatment"
				}
				full_X_matrix_last = full_X_matrix
				mod = tryCatch(private$weibull_generate_mod_from_X(full_X_matrix), error = function(e) NULL)
				if (!is.null(mod)) return(mod)
			}
			# All fast C++ paths failed: one robust survreg fallback on the most-reduced matrix
			surv_mod = robust_survreg_with_surv_object(survival::Surv(private$y, private$dead), full_X_matrix_last)
			if (!is.null(surv_mod)) {
				full_coefficients = c(surv_mod$coefficients, "log(scale)" = log(surv_mod$scale))
				full_vcov = surv_mod$var
				if (!is.null(full_vcov) && is.matrix(full_vcov) && all(is.finite(diag(full_vcov)))) {
					colnames(full_vcov) = rownames(full_vcov) = names(full_coefficients)
					return(list(coefficients = full_coefficients, vcov = full_vcov))
				}
			}
			stop("Weibull regression failed to converge even after progressive correlation dropping and robust retries.")
		}
	)
)
