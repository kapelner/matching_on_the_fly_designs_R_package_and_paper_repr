#' Marginal Standardization / G-Computation for Ordinal Responses
#'
#' @description
#' Internal base class for ordinal-outcome g-computation (G-Comp) estimators. A
#' proportional odds model is fit, then potential-outcome expected values under
#' all-treated and all-control assignments are standardized over the empirical
#' covariate distribution. Inference currently supports bootstrap.
#'
#' @keywords internal
#' @noRd
SeqDesignInferenceOrdinalGCompAbstract = R6::R6Class("SeqDesignInferenceOrdinalGCompAbstract",
	inherit = SeqDesignInference,
	public = list(

		# @description
		# Initialize the g-computation (G-Comp) inference object.
		# @param seq_des_obj A completed \code{SeqDesign} object with an ordinal response.
		# @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		# @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		# @description
		# Computes the g-computation (G-Comp) treatment-effect estimate (mean difference).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$md
		}
	),

	private = list(
		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

		shared = function(){
			if (!is.null(private$cached_values$md)) return(invisible(NULL))

			X_full = private$build_design_matrix()
			# X_full expected to have treatment in column 2, no intercept
			X_fit = X_full[, -1, drop = FALSE] 
			j_treat = 1 # Treatment is now first column of X_fit

			# Fit proportional odds model
			fit = tryCatch(
				fast_ordinal_regression_cpp(X = X_fit, y = as.numeric(private$y)),
				error = function(e) NULL
			)

			if (is.null(fit) || length(fit) == 0){
				private$cached_values$md = NA_real_
				return(invisible(NULL))
			}

			# Post-fit standardization
			res = gcomp_ordinal_proportional_odds_post_fit_cpp(
				X_fit = X_fit,
				coef_hat = as.numeric(fit$b),
				alpha_hat = as.numeric(fit$alpha),
				j_treat = j_treat
			)

			private$cached_values$mean1 = res$mean1
			private$cached_values$mean0 = res$mean0
			private$cached_values$md = res$md
		}
	)
)
