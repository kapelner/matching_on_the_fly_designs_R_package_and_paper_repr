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
SeqDesignInferenceAbstractKKHurdlePoissonIVWC = R6::R6Class("SeqDesignInferenceAbstractKKHurdlePoissonIVWC",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "count")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
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
		}
	),

	private = list(
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()."),

		reduce_design_matrix_preserving_treatment = function(X_full){
			qr_X = qr(X_full)
			target_rank = qr_X$rank
			required = c(1L, 2L)
			candidate_order = c(required, setdiff(qr_X$pivot, required))
			keep = integer(0)

			for (j in candidate_order){
				trial_keep = c(keep, j)
				trial_rank = qr(X_full[, trial_keep, drop = FALSE])$rank
				if (trial_rank > length(keep)){
					keep = trial_keep
				}
				if (length(keep) >= target_rank){
					break
				}
			}

			keep = sort(unique(keep))
			if (!(2L %in% keep)){
				return(list(X = NULL, keep = keep, j_treat = NA_integer_))
			}

			list(
				X = X_full[, keep, drop = FALSE],
				keep = keep,
				j_treat = match(2L, keep)
			)
		},

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
			match_indic = private$match_indic
			if (is.null(match_indic)){
				match_indic = rep(0L, nrow(Xmm))
			}
			match_indic = as.integer(match_indic)
			match_indic[is.na(match_indic)] = 0L

			matched_idx = which(match_indic > 0L)
			reservoir_idx = which(match_indic <= 0L)

			if (length(matched_idx) > 0L){
				private$fit_hurdle_for_matched_pairs(Xmm, matched_idx, match_indic)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) &&
				!is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			if (length(reservoir_idx) > 1L &&
				length(unique(private$w[reservoir_idx])) > 1L){
				private$fit_poisson_for_reservoir(Xmm, reservoir_idx)
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

		fit_hurdle_for_matched_pairs = function(Xmm, matched_idx, match_indic){
			X_matched = Xmm[matched_idx, , drop = FALSE]
			reduced = private$reduce_design_matrix_preserving_treatment(X_matched)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(invisible(NULL))
			}

			pred_df = as.data.frame(X_fit[, -1, drop = FALSE])
			colnames(pred_df)[1] = "w"
			dat = data.frame(
				y = private$y[matched_idx],
				pred_df,
				pair_group = factor(match_indic[matched_idx])
			)

			glmm_control = if (private$num_cores > 1) {
				glmmTMB::glmmTMBControl(parallel = 1L)
			} else {
				glmmTMB::glmmTMBControl()
			}

			formula_cond = private$build_glmm_formula(dat)
			mod = tryCatch(
				suppressWarnings(suppressMessages(
					glmmTMB::glmmTMB(
						formula_cond,
						ziformula = stats::as.formula(sub("^y ~ ", "~ ", deparse(formula_cond))),
						family = glmmTMB::truncated_poisson(link = "log"),
						data = dat,
						control = glmm_control
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
							control = glmm_control
						)
					)),
					error = function(e) NULL
				)
			}
			if (is.null(mod)) return(invisible(NULL))

			coef_table = tryCatch(summary(mod)$coefficients$cond, error = function(e) NULL)
			if (is.null(coef_table) || !("w" %in% rownames(coef_table))) return(invisible(NULL))

			beta_hat = as.numeric(coef_table["w", "Estimate"])
			se = as.numeric(coef_table["w", "Std. Error"])
			if (!is.finite(beta_hat) || !is.finite(se) || se <= 0) return(invisible(NULL))

			private$cached_values$beta_T_matched = beta_hat
			private$cached_values$ssq_beta_T_matched = se^2
		},

		fit_poisson_for_reservoir = function(Xmm, reservoir_idx){
			X_res = Xmm[reservoir_idx, , drop = FALSE]
			reduced = private$reduce_design_matrix_preserving_treatment(X_res)
			X_fit = reduced$X
			if (is.null(X_fit) || !is.finite(reduced$j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(invisible(NULL))
			}

			mod = tryCatch(
				fast_poisson_regression_with_var_cpp(X_fit, private$y[reservoir_idx], j = reduced$j_treat),
				error = function(e) NULL
			)
			if (is.null(mod) || !isTRUE(mod$converged)) return(invisible(NULL))

			beta_hat = as.numeric(mod$b[reduced$j_treat])
			ssq_hat = as.numeric(mod$ssq_b_j)
			if (!is.finite(beta_hat) || !is.finite(ssq_hat) || ssq_hat <= 0) return(invisible(NULL))

			private$cached_values$beta_T_reservoir = beta_hat
			private$cached_values$ssq_beta_T_reservoir = ssq_hat
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("KK hurdle-Poisson IVWC estimator: could not compute a finite standard error.")
			}
		}
	)
)
