#' Univariate GEE Inference for KK Designs with Ordinal Response
#'
#' Fits a Generalized Estimating Equations (GEE) model (using \pkg{multgee})
#' for ordinal responses under a KK matching-on-the-fly design using only
#' the treatment indicator as a predictor (treatment). Matched pairs are
#' treated as clusters (with exchangeable correlation structure); reservoir subjects
#' each form their own singleton cluster. Unlike conditional-logit-only matched-pair
#' analyses, all subjects (matched and reservoir) are included. Inference is based on
#' sandwich-robust standard errors, so the test statistic is Z-distributed.
#'
#' @details
#' This class requires the \pkg{multgee} package, which is listed in Suggests
#' and is not installed automatically with \pkg{SeqExpMatch}.
#' Install \pkg{multgee} before using this class.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUnivKKGEE$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceOrdinalUnivKKGEE = R6::R6Class("InferenceOrdinalUnivKKGEE",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKGEE,
	public = list(

		#' @description
		#' Initialize a univariate GEE inference object for a completed KK design
		#' with an ordinal response.
		#' @param	des_obj		A DesignSeqOneByOne object (must be a KK design) whose entire n subjects
		#' 							are assigned and whose ordinal response y is recorded.
		#' @param num_cores The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference and
		#'   bootstrap resampling.
		#'   The default is 1 for serial computation. For simple
		#'   estimators (e.g. mean difference and KK compound),
		#'   parallelization is achieved with zero-overhead C++ OpenMP.
		#'   For complex models (e.g. GLMs),
		#'   parallelization falls back to R's
		#'   \code{parallel::mclapply}, which incurs
		#'   session-forking overhead.
		#' @param	verbose			Whether to print progress messages. Default is \code{FALSE}.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
			
			if (!requireNamespace("multgee", quietly = TRUE)){
				stop("Package 'multgee' is required for ", class(self)[1], ". Please install it.")
			}
		}
	),
	private = list(
		gee_response_type = function() "ordinal",
		
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Build group ID: matched pairs share their m_vec value;
			# reservoir subjects (m_vec == 0) each get a unique singleton ID.
			group_id = m_vec
			reservoir_idx = which(group_id == 0L)
			if (length(reservoir_idx) > 0L)
				group_id[reservoir_idx] = max(group_id) + seq_along(reservoir_idx)

			dat = data.frame(y = factor(private$y, ordered = TRUE), w = private$w, group_id = group_id)
			# ordgee requires data sorted by id
			dat = dat[order(dat$group_id), ]
			id_sorted = dat$group_id

			mod = tryCatch({
				utils::capture.output(m <- suppressMessages(suppressWarnings(
					multgee::ordLORgee(
						y ~ w,
						data   = dat,
						id     = id_sorted,
						LORstr = "uniform",
						link   = "logit"
					)
				)))
				if (private$verbose) {
					print(summary(m))
				}
				m
			}, error = function(e) {
				if (private$verbose) print(e)
				NULL
			})

			if (is.null(mod)){
				private$cached_values$beta_hat_T   = NA_real_
			if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}
			
			# multgee puts intercepts first, then coefficients
			# The coefficient for 'w' is what we want.
			summ = summary(mod)
			private$cached_values$beta_hat_T   = as.numeric(stats::coef(mod)["w"])
			se = as.numeric(summ$coefficients["w", "san.se"]) # Use sandwich SE
			
			private$cached_values$s_beta_hat_T = if (is.finite(se) && se > 0) se else NA_real_
			private$cached_values$is_z         = TRUE
		}
	)
)
