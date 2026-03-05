#' Abstract class for GLMM-based Inference
#' 
#' @keywords internal
SeqDesignInferenceAbstractKKGLMM = R6::R6Class("SeqDesignInferenceAbstractKKGLMM",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), private$glmm_response_type())
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			if (!requireNamespace("lme4", quietly = TRUE)){
				stop("Package 'lme4' is required for ", class(self)[1], ". Please install it.")
			}
		},

		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the MLE-based confidence interval.
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the MLE-based p-value.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("TO-DO")
			}
		}
	),

	private = list(

		# Abstract: subclasses must return the expected response type string.
		glmm_response_type = function() stop(class(self)[1], " must implement glmm_response_type()"),

		# Abstract: subclasses must return the glm family object for glmer.
		glmm_family = function() stop(class(self)[1], " must implement glmm_family()"),

		# Default (multivariate): intercept dropped, treatment column named "w".
		# Univariate subclasses override this to return data.frame(w = private$w).
		glmm_predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			mod = private$fit_glmm()
			if (is.null(mod)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}
			coef_table = summary(mod)$coefficients
			private$cached_values$beta_hat_T   = as.numeric(coef_table["w", "Estimate"])
			se = as.numeric(coef_table["w", "Std. Error"])
			# Store NA when the SE is non-finite; SE-dependent methods detect this via assert_finite_se()
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z         = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				stop("GLMM: non-finite standard error (possible separation or singular random-effect fit).")
		},

		fit_glmm = function(){
			match_indic = private$match_indic
			if (is.null(match_indic)) match_indic = rep(0L, private$n)
			match_indic[is.na(match_indic)] = 0L

			# Build group ID: matched pairs share their match_indic value;
			# reservoir subjects (match_indic == 0) each get a unique singleton ID.
			group_id = match_indic
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			dat = data.frame(y = private$y, private$glmm_predictors_df(), group_id = factor(group_id))

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					lme4::glmer(
						y ~ . - group_id + (1 | group_id),
						family = private$glmm_family(),
						data   = dat
					)
				)))
				mod
			}, error = function(e) NULL)
		}
	)
)
