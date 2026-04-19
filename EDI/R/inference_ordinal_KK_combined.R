#' GEE Inference for KK Designs with Ordinal Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{multgee})
#' for ordinal responses under a KK matching-on-the-fly design using the
#' treatment indicator and, optionally, all recorded covariates as predictors.
#'
#' @export
InferenceOrdinalKKGEE = R6::R6Class("InferenceOrdinalKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(
		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed \code{Design} object with an ordinal response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				if (!check_package_installed("multgee")){
					stop("Package 'multgee' is required for ", class(self)[1], ". Please install it.")
				}
			}
			super$initialize(des_obj, include_covariates, verbose)
		}
	),
	private = list(
		gee_response_type = function() "ordinal",

		# Override: ordinal response requires ordLORgee, not geeglm.
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			pred_df = private$gee_predictors_df()
			dat = data.frame(y = factor(private$y, ordered = TRUE), pred_df, group_id = group_id)
			dat = dat[order(dat$group_id), ]
			id_sorted = dat$group_id
			
			fixed_terms = setdiff(colnames(dat), c("y", "group_id"))
			formula_gee = stats::as.formula(paste("y ~", paste(fixed_terms, collapse = " + ")))

			mod = tryCatch({
				utils::capture.output(m <- suppressMessages(suppressWarnings(
					multgee::ordLORgee(
						formula_gee,
						data   = dat,
						id     = id_sorted,
						LORstr = "uniform",
						link   = "logit"
					)
				)))
				m
			}, error = function(e) NULL)

			if (is.null(mod)){
				private$cache_nonestimable_estimate("ordinal_kk_gee_fit_unavailable")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}

			beta = stats::coef(mod)
			j_treat = private$gee_treatment_index(beta)
			private$cached_values$beta_hat_T = as.numeric(beta[j_treat])

			if (estimate_only) return(invisible(NULL))

			vcov_robust = tryCatch(stats::vcov(mod), error = function(e) NULL)
			if (is.null(vcov_robust)) {
				private$cached_values$s_beta_hat_T = NA_real_
			} else {
				private$cached_values$s_beta_hat_T = sqrt(as.numeric(vcov_robust[j_treat, j_treat]))
			}
			private$cached_values$is_z = TRUE
			private$cached_values$df = Inf
			private$cached_values$summary_table = summary(mod)$coefficients
		}
	)
)

#' GLMM Inference for KK Designs with Ordinal Response
#'
#' Fits a cumulative-link mixed model (using \pkg{glmmTMB}) for ordinal responses
#' under a KK matching-on-the-fly design using the treatment indicator and,
#' optionally, all recorded covariates as fixed-effect predictors.
#'
#' @export
InferenceOrdinalKKGLMM = R6::R6Class("InferenceOrdinalKKGLMM",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGLMM,
	public = list(
	),
	private = list(
		glmm_response_type  = function() "ordinal",
		glmm_family         = function() glmmTMB::cumulative(link = "logit")
	)
)
