#' Multivariate GEE Inference for KK Designs with Ordinal Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (via \code{multgee::ordLORgee})
#' for ordinal responses under a KK matching-on-the-fly design, adjusting for
#' baseline covariates.
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$new(n = nrow(x_dat), response_type = "ordinal",
#' verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_one_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalMultiKKGEE$
#'   new(seq_des, verbose = FALSE)
#' infer
#' }
#'
InferenceOrdinalMultiKKGEE = R6::R6Class("InferenceOrdinalMultiKKGEE",
	lock_objects = FALSE,
	inherit = InferenceOrdinalUnivKKGEE,
	public = list(
	),
	private = list(
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			dat = data.frame(y = factor(private$y, ordered = TRUE), w = private$w,
				group_id = group_id)
			X = private$get_X()
			dat = cbind(dat, X)
			
			# Ensure data is sorted by group_id for multgee
			dat = dat[order(dat$group_id), ]
			id_sorted = dat$group_id

			# Build formula including all covariates
			x_names = colnames(X)
			if (is.null(x_names)) x_names = paste0("V", seq_len(ncol(X)))
			formula_str = paste("y ~ w +", paste(x_names, collapse = " + "))
			
			mod = tryCatch({
				utils::capture.output(m <- suppressMessages(suppressWarnings(
					multgee::ordLORgee(
						as.formula(formula_str),
						data   = dat,
						id     = id_sorted,
						LORstr = "uniform",
						link   = "logit"
					)
				)))
				m
			}, error = function(e) NULL)

			if (is.null(mod)){
				private$cached_values$beta_hat_T   = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}
			
			summ = summary(mod)
			private$cached_values$beta_hat_T   = as.numeric(stats::coef(mod)["w"])
			se = as.numeric(summ$coefficients["w", "san.se"])
			
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z         = TRUE
		}
	)
)
