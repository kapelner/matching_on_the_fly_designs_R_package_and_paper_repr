#' Abstract class for GLMM-based Inference
#'
#' @keywords internal
InferenceAbstractKKGLMM = R6::R6Class("InferenceAbstractKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		#' @param make_fork_cluster Whether to use a fork cluster for parallelization.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), private$glmm_response_type())
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
			if (!requireNamespace("glmmTMB", quietly = TRUE)){
				stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
			}
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("TO-DO")
			}
		},

		#' @description
		#' Overridden to provide a warning about slowness.
		#' @param r		Number of randomization iterations.
		#' @param ... 					Additional arguments passed to super.
		compute_confidence_interval_rand = function(r = 501, ...){
			warning("Randomization-based confidence intervals for GLMM models are extremely slow because they require many hundreds or thousands of model fits. Consider using asymptotic or bootstrap intervals instead.")
			super$compute_confidence_interval_rand(r = r, ...)
		},

		#' @description
		#' Overridden to provide a warning about slowness.
		#' @param r		Number of randomization iterations.
		#' @param ... 					Additional arguments passed to super.
		compute_two_sided_pval_for_treatment_effect_rand = function(r = 501, ...){
			warning("Randomization-based p-values for GLMM models are slow because they require many model fits.")
			super$compute_two_sided_pval_for_treatment_effect_rand(r = r, ...)
		}
	),

	private = list(

		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the fixed-effect coefficient for "w" directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(){
			mod = private$fit_glmm(se = FALSE)
			if (is.null(mod)) return(NA_real_)
			
			# glmmTMB fixed effects for the conditional model
			beta = glmmTMB::fixef(mod)$cond
			if ("w" %in% names(beta)){
				return(as.numeric(beta["w"]))
			}
			NA_real_
		},

		# Abstract: subclasses must return the expected response type string.
		glmm_response_type = function() stop(class(self)[1], " must implement glmm_response_type()"),

		# Abstract: subclasses must return the glm family object for glmmTMB.
		glmm_family = function() stop(class(self)[1], " must implement glmm_family()"),

		# Default (multivariate): intercept dropped, treatment column named "w".
		# Univariate subclasses override this to return data.frame(w = private$w).
		glmm_predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			mod = private$fit_glmm(se = TRUE)
			if (is.null(mod)){
				private$cached_values$beta_hat_T   = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}
			coef_table = summary(mod)$coefficients$cond
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

		fit_glmm = function(se = TRUE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Build group ID: matched pairs share their m_vec value;
			# reservoir subjects (m_vec == 0) each get a unique singleton ID.
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			dat = data.frame(y = private$y, private$glmm_predictors_df(), group_id = factor(group_id))
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			glmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))

			# Respect the core budget provided by the user. If we are in an outer 
			# parallel loop (e.g. mclapply), private$num_cores will be 1L.
			glmm_control = glmmTMB::glmmTMBControl(parallel = private$num_cores)

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					glmmTMB::glmmTMB(
						glmm_formula,
						family  = private$glmm_family(),
						data    = dat,
						control = glmm_control,
						se      = se
					)
				)))
				mod
			}, error = function(e) NULL)
		}
	)
)
