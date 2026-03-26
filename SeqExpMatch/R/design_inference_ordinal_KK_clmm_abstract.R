# Abstract class for ordinal CLMM-based Inference in KK designs
#
# @keywords internal
DesignInferenceAbstractKKOrdinalCLMM = R6::R6Class("DesignInferenceAbstractKKOrdinalCLMM",
	inherit = DesignInferenceKKPassThrough,
	public = list(

		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "ordinal")
			if (!is(des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			if (!requireNamespace("ordinal", quietly = TRUE)){
				stop("Package 'ordinal' is required for ", class(self)[1], ". Please install it.")
			}
		},

		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

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
		clmm_link = function() stop(class(self)[1], " must implement clmm_link()"),

		clmm_predictors_df = function(){
			full_X = private$create_design_matrix()
			X_model = full_X[, -1, drop = FALSE]
			colnames(X_model)[1] = "w"
			as.data.frame(X_model)
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			mod = private$fit_clmm()
			summ = if (!is.null(mod)) tryCatch(summary(mod), error = function(e) NULL) else NULL
			se = if (!is.null(summ)) as.numeric(summ$coefficients["w", "Std. Error"]) else NA_real_
			if (is.null(mod) || !is.finite(se) || se <= 0){
				mod = private$fit_clm_fallback()
				summ = if (!is.null(mod)) tryCatch(summary(mod), error = function(e) NULL) else NULL
				se = if (!is.null(summ)) as.numeric(summ$coefficients["w", "Std. Error"]) else NA_real_
			}
			if (is.null(mod) || is.null(summ)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = as.numeric(stats::coef(mod)["w"])
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				stop("CLMM: non-finite standard error (possible separation or singular random-effect fit).")
		},

		fit_clmm = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			dat = data.frame(
				y = factor(private$y, ordered = TRUE),
				private$clmm_predictors_df(),
				group_id = factor(group_id)
			)
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			clmm_formula = stats::as.formula(paste("y ~", paste(c(fixed_terms, "(1 | group_id)"), collapse = " + ")))

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					ordinal::clmm(
						clmm_formula,
						data = dat,
						link = private$clmm_link()
					)
				)))
				mod
			}, error = function(e) NULL)
		},

		fit_clm_fallback = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			dat = data.frame(
				y = factor(private$y, ordered = TRUE),
				private$clmm_predictors_df()
			)
			fixed_terms = setdiff(colnames(dat), "y")
			clm_formula = stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))

			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					ordinal::clm(
						clm_formula,
						data = dat,
						link = private$clmm_link()
					)
				)))
				mod
			}, error = function(e) NULL)
		}
	)
)
