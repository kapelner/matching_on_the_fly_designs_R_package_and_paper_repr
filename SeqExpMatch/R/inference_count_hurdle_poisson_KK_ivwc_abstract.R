#' KK Hurdle Poisson IVWC Inference for Count Responses
#'
#' @description
#' Internal base class for KK hurdle-Poisson inverse-variance weighted combined
#' inference. The matched-pair component is fit with a hurdle-Poisson mixed model
#' using pair random intercepts, and the reservoir component is fit with an
#' ordinary Poisson log-link regression. The reported treatment effect is on the
#' log-rate scale.
#'
#' @keywords internal
#' @noRd
InferenceAbstractKKHurdlePoissonIVWC = R6::R6Class("InferenceAbstractKKHurdlePoissonIVWC",
	inherit = InferenceAsymp,
	public = list(

		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "count")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			if (!requireNamespace("glmmTMB", quietly = TRUE)){
				stop("Package 'glmmTMB' is required for ", class(self)[1], ". Please install it.")
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
				stop("Testing non-zero delta is not yet implemented for this class.")
			}
		},

		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the fixed-effect coefficient for "w" directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(){
			Xmm = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(0L, nrow(Xmm))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L

			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)

			beta_m = NA_real_
			ssq_m = NA_real_
			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(Xmm, matched_idx, m_vec, se = FALSE)
				beta_m = res_m$beta_hat
				ssq_m = res_m$se^2
			}
			m_ok = !is.na(beta_m) && is.finite(beta_m) && !is.na(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			beta_r = NA_real_
			ssq_r = NA_real_
			if (length(reservoir_idx) > 1L && length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(Xmm, reservoir_idx) # Already uses fast C++ which includes var, but we can't easily skip it
				beta_r = res_r$beta_hat
				ssq_r = res_r$ssq_hat
			}
			r_ok = !is.na(beta_r) && is.finite(beta_r) && !is.na(ssq_r) && is.finite(ssq_r) && ssq_r > 0

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				return(w_star * beta_m + (1 - w_star) * beta_r)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		}
	),

	private = list(
		# ... (other methods)

		build_model_matrix = function(){
			if (private$include_covariates()){
				Xmm = private$create_design_matrix()
				colnames(Xmm) = c("(Intercept)", "w", if (ncol(Xmm) > 2L) paste0("x", seq_len(ncol(Xmm) - 2L)) else NULL)
			} else {
				Xmm = cbind(1, private$w)
				colnames(Xmm) = c("(Intercept)", "w")
			}
			Xmm
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			Xmm = private$build_model_matrix()
			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(0L, nrow(Xmm))
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L

			matched_idx = which(m_vec > 0L)
			reservoir_idx = which(m_vec <= 0L)

			if (length(matched_idx) > 0L){
				res_m = private$fit_hurdle_for_matched_pairs(Xmm, matched_idx, m_vec, se = TRUE)
				private$cached_values$beta_T_matched = res_m$beta_hat
				private$cached_values$ssq_beta_T_matched = res_m$se^2
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
				!is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			if (length(reservoir_idx) > 1L &&
				length(unique(private$w[reservoir_idx])) > 1L){
				res_r = private$fit_poisson_for_reservoir(Xmm, reservoir_idx)
				private$cached_values$beta_T_reservoir = res_r$beta_hat
				private$cached_values$ssq_beta_T_reservoir = res_r$ssq_hat
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) &&
				!is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		build_glmm_formula = function(dat){
			fixed_terms = setdiff(colnames(dat), c("y", "pair_group"))
			rhs = paste(c(fixed_terms, "(1 | pair_group)"), collapse = " + ")
			stats::as.formula(paste("y ~", rhs))
		},

		fit_hurdle_for_matched_pairs = function(Xmm, matched_idx, m_vec, se = TRUE){
			X_matched = Xmm[matched_idx, , drop = FALSE]
			reduced = private$reduce_design_matrix_preserving_treatment(X_matched)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(beta_hat = NA_real_, se = NA_real_))
			}

			pred_df = as.data.frame(X_fit[, -1, drop = FALSE])
			colnames(pred_df)[1] = "w"
			dat = data.frame(
				y = private$y[matched_idx],
				pred_df,
				pair_group = factor(m_vec[matched_idx])
			)

			glmm_control = glmmTMB::glmmTMBControl(parallel = private$num_cores)

			formula_cond = private$build_glmm_formula(dat)
			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = stats::as.formula(sub("^y ~ ", "~ ", deparse(formula_cond))),
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat,
						control = glmm_control,
						se = se
					)
				)),
				error = function(e) NULL
			)
			if (is.null(mod) && ncol(dat) > 3L){
				dat = dat[, c("y", "w", "pair_group"), drop = FALSE]
				mod = tryCatch(
					suppressWarnings(suppressMessages(
						glmmTMB::glmmTMB(
							y ~ w + (1 | pair_group),
							ziformula = ~ w + (1 | pair_group),
							family = glmmTMB::truncated_poisson(link = "log"),
							data = dat,
							control = glmm_control,
							se = se
						)
					)),
					error = function(e) NULL
				)
			}
			if (is.null(mod)) return(list(beta_hat = NA_real_, se = NA_real_))

			if (!se){
				beta = glmmTMB::fixef(mod)$cond
				if ("w" %in% names(beta)){
					return(list(beta_hat = as.numeric(beta["w"]), se = 1.0)) # Return dummy SE > 0
				}
				return(list(beta_hat = NA_real_, se = NA_real_))
			}

			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table))) return(list(beta_hat = NA_real_, se = NA_real_))

			beta_hat = as.numeric(coef_table["w", "Estimate"])
			se_val = as.numeric(coef_table["w", "Std. Error"])
			if (!is.finite(beta_hat) || !is.finite(se_val) || se_val <= 0) return(list(beta_hat = NA_real_, se = NA_real_))

			list(beta_hat = beta_hat, se = se_val)
		},

		fit_poisson_for_reservoir = function(Xmm, reservoir_idx){
			X_res = Xmm[reservoir_idx, , drop = FALSE]
			reduced = private$reduce_design_matrix_preserving_treatment(X_res)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(beta_hat = NA_real_, ssq_hat = NA_real_))
			}

			mod = tryCatch(
				fast_poisson_regression_with_var_cpp(X_fit, private$y[reservoir_idx], j = reduced$j_treat),
				error = function(e) NULL
			)
			if (is.null(mod) || !isTRUE(mod$converged)) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))

			beta_hat = as.numeric(mod$b[reduced$j_treat])
			ssq_hat = as.numeric(mod$ssq_b_j)
			if (!is.finite(beta_hat) || !is.finite(ssq_hat) || ssq_hat <= 0) return(list(beta_hat = NA_real_, ssq_hat = NA_real_))

			list(beta_hat = beta_hat, ssq_hat = ssq_hat)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("KK hurdle-Poisson IVWC estimator: could not compute a finite standard error.")
			}
		}
	)
)
