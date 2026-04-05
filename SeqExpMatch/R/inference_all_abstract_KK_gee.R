#' Abstract class for GEE-based Inference
#'
#' @keywords internal
InferenceAbstractKKGEE = R6::R6Class("InferenceAbstractKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), private$gee_response_type())
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, verbose)
			private$m = des_obj$.__enclos_env__$private$m
			if (identical(private$gee_response_type(), "proportion")) {
				private$y = .sanitize_proportion_response(private$y, interior = FALSE)
			}
			assertNoCensoring(private$any_censoring)
			if (!requireNamespace("geepack", quietly = TRUE)){
				stop("Package 'geepack' is required for ", class(self)[1], ". Please install it.")
			}
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
				private$shared(estimate_only = TRUE)
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
		}
	),

	private = list(

		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the coefficient for "w" directly from the fit.
		# Overridden to avoid the heavy summary() call during randomization iterations.
		# Extracts the coefficient for treatment directly from the fit.
		compute_treatment_estimate_during_randomization_inference = function(){
			mod = private$fit_gee_with_fallback(std_err = FALSE, estimate_only = TRUE)
			private$extract_gee_treatment_estimate(mod)
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

		gee_treatment_index = function(beta){
			if (is.null(beta) || !length(beta)) return(NA_integer_)
			beta_names = names(beta)
			if (!is.null(beta_names) && ("w" %in% beta_names)) return(match("w", beta_names))
			if (length(beta) >= 2L) return(2L)
			NA_integer_
		},

		extract_gee_treatment_estimate = function(mod){
			if (is.null(mod)) return(NA_real_)
			beta = tryCatch(stats::coef(mod), error = function(e) NULL)
			j_treat = private$gee_treatment_index(beta)
			if (!is.finite(j_treat) || is.na(j_treat) || j_treat < 1L || j_treat > length(beta)) return(NA_real_)
			as.numeric(beta[[j_treat]])
		},

		extract_gee_treatment_se = function(mod, j_treat = NA_integer_, coef_table = NULL){
			if (is.null(mod)) return(NA_real_)
			beta = tryCatch(stats::coef(mod), error = function(e) NULL)
			if (is.na(j_treat) || !is.finite(j_treat)) j_treat = private$gee_treatment_index(beta)
			if (!is.finite(j_treat) || is.na(j_treat) || j_treat < 1L) return(NA_real_)

			if (!is.null(coef_table)) {
				se_col = intersect(c("Std.err", "Std.error", "Robust S.E."), colnames(coef_table))
				if (length(se_col) > 0L) {
					row_idx = if (!is.null(rownames(coef_table)) && ("w" %in% rownames(coef_table))) match("w", rownames(coef_table)) else j_treat
					if (is.finite(row_idx) && !is.na(row_idx) && row_idx >= 1L && row_idx <= nrow(coef_table)) {
						se = suppressWarnings(as.numeric(coef_table[row_idx, se_col[1]]))
						if (is.finite(se) && se > 0) return(se)
					}
				}
			}

			vc = tryCatch(stats::vcov(mod), error = function(e) NULL)
			if (!is.null(vc) && is.matrix(vc) && j_treat <= nrow(vc) && j_treat <= ncol(vc)) {
				v = suppressWarnings(as.numeric(vc[j_treat, j_treat]))
				if (is.finite(v) && v > 0) return(sqrt(v))
			}
			NA_real_
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (estimate_only){
				mod = private$fit_gee_with_fallback(std_err = FALSE, estimate_only = TRUE)
				private$cached_values$beta_hat_T = private$extract_gee_treatment_estimate(mod)
				return(invisible(NULL))
			}

			mod = private$fit_gee_with_fallback(std_err = TRUE, estimate_only = FALSE)
			if (is.null(mod)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			beta = tryCatch(stats::coef(mod), error = function(e) NULL)
			j_treat = private$gee_treatment_index(beta)
			private$cached_values$beta_hat_T = if (is.finite(j_treat) && !is.na(j_treat) && j_treat >= 1L && j_treat <= length(beta)) as.numeric(beta[[j_treat]]) else NA_real_

			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			private$cached_values$s_beta_hat_T = private$extract_gee_treatment_se(mod, j_treat = j_treat, coef_table = coef_table)
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T))
				return(invisible(NULL))
		},
		gee_has_reservoir = function(){
			m_vec = private$m
			if (is.null(m_vec)) return(FALSE)
			m_vec[is.na(m_vec)] = 0L
			any(m_vec == 0L)
		},

		build_gee_fit_data = function(include_reservoir = TRUE){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			keep = if (isTRUE(include_reservoir)) rep(TRUE, length(m_vec)) else m_vec > 0L
			if (!any(keep)) return(NULL)

			m_keep = m_vec[keep]
			group_id = m_keep
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L) {
				max_group = if (all(group_id == 0L)) 0L else max(group_id)
				group_id[reservoir_idx] = max_group + seq_along(reservoir_idx)
			}

			pred_df = private$gee_predictors_df()
			if (!is.data.frame(pred_df)) pred_df = as.data.frame(pred_df)
			dat = data.frame(y = private$y[keep], pred_df[keep, , drop = FALSE], group_id = group_id)
			dat = dat[order(dat$group_id), , drop = FALSE]
			id_sorted = dat$group_id
			dat$group_id = NULL
			list(dat = dat, id_sorted = id_sorted)
		},

		fit_gee_on_data = function(dat, id_sorted, std_err = TRUE){
			std_err_arg = if (is.character(std_err)) std_err[1] else "san.se"
			tryCatch({
				utils::capture.output(mod <- suppressMessages(suppressWarnings(
					geepack::geeglm(
						y ~ .,
						family = private$gee_family(),
						data   = dat,
						id     = id_sorted,
						corstr = "exchangeable",
						std.err = std_err_arg
					)
				)))
				mod
			}, error = function(e) NULL)
		},

		fit_gee = function(std_err = TRUE, include_reservoir = TRUE){
			fit_data = private$build_gee_fit_data(include_reservoir = include_reservoir)
			if (is.null(fit_data)) return(NULL)
			private$fit_gee_on_data(fit_data$dat, fit_data$id_sorted, std_err = std_err)
		},

		fit_gee_with_fallback = function(std_err = TRUE, estimate_only = FALSE){
			mod = private$fit_gee(std_err = std_err, include_reservoir = TRUE)
			needs_fallback = is.null(mod)
			if (!needs_fallback) {
				beta = tryCatch(stats::coef(mod), error = function(e) NULL)
				beta_hat = private$extract_gee_treatment_estimate(mod)
				needs_fallback = !is.finite(beta_hat)
				if (!estimate_only && !needs_fallback) {
					j_treat = private$gee_treatment_index(beta)
					coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
					se_hat = private$extract_gee_treatment_se(mod, j_treat = j_treat, coef_table = coef_table)
					needs_fallback = !is.finite(se_hat) || se_hat <= 0
				}
			}
			if (needs_fallback && private$gee_has_reservoir()) {
				mod_fb = private$fit_gee(std_err = std_err, include_reservoir = FALSE)
				if (!is.null(mod_fb)) return(mod_fb)
			}
			mod
		}
	)
)
