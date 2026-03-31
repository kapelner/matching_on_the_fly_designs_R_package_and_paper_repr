#' Abstract class for Weibull Frailty Combined-Likelihood Inference
#'
#' Fits a single joint Weibull frailty model over all KK design data for survival
#' responses. Matched subjects share their pair ID as a frailty cluster; reservoir
#' subjects are assigned to singleton clusters. The treatment effect beta_T and
#' covariate slopes beta_xs are shared across all subjects.
#'
#' @keywords internal
InferenceAbstractKKWeibullFrailtyCombinedLikelihood = R6::R6Class("InferenceAbstractKKWeibullFrailtyCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK design object.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "survival")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, verbose)
		},

		#' @description
		#' Returns the combined-likelihood estimate of the treatment effect.
		#' @param estimate_only Whether to skip standard-error calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an asymptotic confidence interval for the treatment effect.
		#' @param alpha Significance level.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Returns a 2-sided p-value for H0: beta_T = delta.
		#' @param delta Null treatment effect value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		build_design_matrix = function(){
			X_full = matrix(private$w, ncol = 1)
			colnames(X_full) = "w"
			if (private$include_covariates()){
				X_full = cbind(X_full, as.matrix(private$X))
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(1L %in% keep)) keep[r_full] = 1L
					keep = sort(keep)
					X_full = X_full[, keep, drop = FALSE]
				}
			}
			X_full
		},

		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			cluster_id = m_vec
			reservoir_idx = which(cluster_id == 0L)
			if (length(reservoir_idx) > 0L){
				cluster_id[reservoir_idx] = max(cluster_id) + seq_along(reservoir_idx)
			}

			if (sum(private$dead) == 0L){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			X_full = private$build_design_matrix()
			dat = data.frame(y = private$y, dead = private$dead, w = X_full[, "w"], cluster = factor(cluster_id))
			formula_str = "survival::Surv(y, dead) ~ w"

			X_covs = X_full[, colnames(X_full) != "w", drop = FALSE]
			if (ncol(X_covs) > 0){
				colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
				dat = cbind(dat, X_covs)
				formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
			}
			formula_str = paste(formula_str, "+ cluster(cluster)")

			mod = tryCatch(
				suppressWarnings(survival::survreg(as.formula(formula_str), data = dat, dist = "weibull")),
				error = function(e) NULL
			)

			if (is.null(mod)){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
			if (is.null(coefs) || !("w" %in% names(coefs))){
				private$cached_values$beta_hat_T = NA_real_
				if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			beta = as.numeric(coefs["w"])
			private$cached_values$beta_hat_T = if (is.finite(beta)) beta else NA_real_
			
			if (!estimate_only) {
				vcv = tryCatch(stats::vcov(mod), error = function(e) NULL)
				if (is.null(vcv) || !("w" %in% rownames(vcv))){
					private$cached_values$s_beta_hat_T = NA_real_
				} else {
					se = sqrt(as.numeric(vcv["w", "w"]))
					private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
				}
			}
			private$cached_values$is_z         = TRUE
			invisible(NULL)
		}
	)
)
