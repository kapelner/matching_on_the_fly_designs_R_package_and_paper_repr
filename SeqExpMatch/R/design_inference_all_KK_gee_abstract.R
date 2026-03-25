# Abstract class for GEE-based Inference
#
# @keywords internal
DesignInferenceAbstractKKGEE = R6::R6Class("DesignInferenceAbstractKKGEE",
	inherit = DesignInferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), private$gee_response_type())
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			if (!requireNamespace("geepack", quietly = TRUE)){
				stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
			}
		},

		# @description
		# Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes the asymptotic confidence interval.
		# @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes the asymptotic p-value.
		# @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
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

		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the coefficient for "w" directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(){
			mod = private$fit_gee(std_err = FALSE)
			if (is.null(mod)) return(NA_real_)
			
			beta = coef(mod)
			if ("w" %in% names(beta)){
				return(as.numeric(beta["w"]))
			}
			NA_real_
		},

		# Abstract: subclasses must return the expected response type string.
		gee_response_type = function() stop(class(self)[1], " must implement gee_response_type()"),

		# Abstract: subclasses must return the glm family object for geeglm.
		gee_family = function() stop(class(self)[1], " must implement gee_family()"),

		# Default (multivariate): intercept dropped, treatment column named "w".
		# Univariate subclasses override this to return data.frame(w = private$w).
		gee_predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			mod = private$fit_gee(std_err = TRUE)
			if (is.null(mod)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}
			
			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table))){
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			} else {
				private$cached_values$beta_hat_T   = as.numeric(coef_table["w", "Estimate"])
				se = as.numeric(coef_table["w", "Std.err"])
				private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				stop("GEE: non-finite standard error (possible separation or insufficient data).")
		},

		fit_gee = function(std_err = TRUE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Build group ID: matched pairs share their m_vec value;
			# reservoir subjects (m_vec == 0) each get a unique singleton ID.
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			dat = data.frame(y = private$y, private$gee_predictors_df(), group_id = group_id)
			# geeglm requires data sorted by id
			dat = dat[order(dat$group_id), ]
			id_sorted = dat$group_id
			dat$group_id = NULL

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					geepack::geeglm(
						y ~ .,
						family = private$gee_family(),
						data   = dat,
						id     = id_sorted,
						corstr = "exchangeable",
						std.err = std_err
					)
				)))
				mod
			}, error = function(e) NULL)
		}
	)
)
