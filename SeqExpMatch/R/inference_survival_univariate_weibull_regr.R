#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#'
#'
#' @export
SeqDesignInferenceSurvivalUniWeibullRegr = R6::R6Class("SeqDesignInferenceSurvivalUniWeibullRegr",
	inherit = SeqDesignInferenceMLEorKMSummaryTable,
	public = list(

		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 								(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 								for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			assertResponseType(seq_des_obj$get_response_type(), "survival")
		}
	),

	private = list(
		generate_mod = function(){
			# Univariate: treatment only, no covariates (mirrors SeqDesignInferenceSurvivalUniCoxPHRegr)
			full_X_matrix = matrix(private$w, ncol = 1)
			colnames(full_X_matrix) = "treatment"
			private$weibull_generate_mod_from_X(full_X_matrix)
		},

		weibull_generate_mod_from_X = function(full_X_matrix){
			weibull_regr_mod = fast_weibull_regression(
				private$y,
				private$dead,
				as.matrix(full_X_matrix)
			)
			# fast_weibull_regression already names coefficients correctly (only retained columns after
			# collinearity dropping), so no name re-assignment is needed here.

			if (is.null(weibull_regr_mod$coefficients) || is.null(weibull_regr_mod$vcov) || !is.matrix(weibull_regr_mod$vcov)){
				stop("fast_weibull_regression failed to return valid coefficients or vcov.")
			}

			full_coefficients = c(weibull_regr_mod$coefficients, "log(scale)" = weibull_regr_mod$log_sigma)
			full_vcov = weibull_regr_mod$vcov
			colnames(full_vcov) = rownames(full_vcov) = names(full_coefficients)

			list(coefficients = full_coefficients, vcov = full_vcov)
		}
	)
)
