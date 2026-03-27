#' Simple Mean Difference Inference based on Maximum Likelihood
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring)
#' sequential experimental design estimation and test object
#' after the sequential design is completed.
#'
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
#'   response_type = "survival",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(
#'   c(1.2, 2.4, 1.8, 3.1, 2.7, 4.0, 3.3, 4.5),
#'   c(1, 1, 0, 1, 0, 1, 1, 0)
#' )
#' infer <- InferenceSurvivalUniWeibullRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceSurvivalUniWeibullRegr = R6::R6Class("InferenceSurvivalUniWeibullRegr",
	inherit = InferenceMLEorKMSummaryTable,
	public = list(

	),

	private = list(
		generate_mod = function(){
			# Univariate: treatment only, no covariates (mirrors InferenceSurvivalUniCoxPHRegr)
			full_X_matrix = matrix(private$w, ncol = 1)
			colnames(full_X_matrix) = "treatment"
			mod = tryCatch(private$weibull_generate_mod_from_X(full_X_matrix), error = function(e) NULL)
			if (!is.null(mod)) return(mod)
			# Fast path failed: fall back to robust_survreg with multiple random initializations
			surv_mod = robust_survreg_with_surv_object(survival::Surv(private$y, private$dead), full_X_matrix)
			if (!is.null(surv_mod)) {
				full_coefficients = c(surv_mod$coefficients, "log(scale)" = log(surv_mod$scale))
				full_vcov = surv_mod$var
				if (!is.null(full_vcov) && is.matrix(full_vcov) && all(is.finite(diag(full_vcov)))) {
					colnames(full_vcov) = rownames(full_vcov) = names(full_coefficients)
					return(list(coefficients = full_coefficients, vcov = full_vcov))
				}
			}
			stop("Weibull regression failed to converge even after robust retries.")
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
