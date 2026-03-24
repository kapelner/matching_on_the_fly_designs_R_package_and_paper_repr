#' Multivariate Conditional Continuation-Ratio Inference for KK Designs
#'
#' @description
#' Fits a conditional continuation-ratio logit model for ordinal responses under a KK
#' matching-on-the-fly design, adjusting for baseline covariates.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignKK14$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalMultiKKCondContRatioRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceOrdinalMultiKKCondContRatioRegr = R6::R6Class(
	"SeqDesignInferenceOrdinalMultiKKCondContRatioRegr",
	inherit = SeqDesignInferenceOrdinalUnivKKCondContRatioRegr,
	public = list(
		#' @description
		#' Initialize a multivariate conditional continuation-ratio inference object.
		#' @param	seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param	num_cores			Number of CPU cores.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		}
	),

	private = list(
		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			# Define strata
			strata_ids = m_vec
			reservoir_idx = which(strata_ids == 0L)
			if (length(reservoir_idx) > 0L){
				strata_ids[reservoir_idx] = max(strata_ids) + seq_along(reservoir_idx)
			}

			y_ord = as.integer(factor(private$y, ordered = TRUE))
			K = max(y_ord)
			if (K < 2L){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			n_alpha = K - 1L
			n = private$n
			num_strata = max(strata_ids)
			X = private$get_X()

			# Data expansion for continuation ratio (multivariate)
			expanded = expand_continuation_ratio_data_cpp(as.integer(y_ord),
				as.integer(private$w), as.integer(strata_ids), as.integer(K))
			
			# Subject i contributes trials j = 1, ..., min(yi, n_alpha)
			X_stack_list = list()
			for (i in 1:n) {
				yi = y_ord[i]
				n_trials_i = min(yi, n_alpha)
				if (n_trials_i > 0) {
					X_stack_list[[i]] = matrix(rep(X[i, ], each = n_trials_i), nrow = n_trials_i)
				}
			}
			X_stack = do.call(rbind, X_stack_list)
			
			# Fit conditional logistic regression
			mod = clogit_helper(expanded$y, as.data.frame(X_stack), expanded$w, expanded$strata)
			
			if (is.null(mod) || !is.finite(mod$b[1]) || !is.finite(mod$ssq_b_j) || mod$ssq_b_j <= 0){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(mod$b[1])
			private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z         = TRUE
		}
	)
)
