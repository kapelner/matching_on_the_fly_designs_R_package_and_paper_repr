#' Multivariate Conditional Adjacent-Category Inference for KK Designs
#'
#' Fits a conditional adjacent-category logit model for ordinal responses under a KK
#' matching-on-the-fly design, adjusting for baseline covariates.
#'
#' @export
#' @examples
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
#' infer <- InferenceOrdinalMultiKKCondAdjCatLogitRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
#'
InferenceOrdinalMultiKKCondAdjCatLogitRegr = R6::R6Class(
	"InferenceOrdinalMultiKKCondAdjCatLogitRegr",
	lock_objects = FALSE,
	inherit = InferenceOrdinalUnivKKCondAdjCatLogitRegr,
	public = list(
	),

	private = list(
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			strata_ids = m_vec
			reservoir_idx = which(strata_ids == 0L)
			if (length(reservoir_idx) > 0L){
				strata_ids[reservoir_idx] = max(strata_ids) + seq_along(reservoir_idx)
			}

			y_ord = as.integer(factor(private$y, ordered = TRUE))
			K = max(y_ord)
			if (K < 2L) return(NA_real_)

			expanded = expand_adjacent_category_data_cpp(as.integer(y_ord), as.integer(private$w), as.integer(strata_ids), as.integer(K))
			
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			
			X_stack_list = list()
			for (i in 1:private$n) {
				yi = y_ord[i]
				Xi = X_data[i, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				trials_i = integer()
				for (j in 1 : (K-1L)) {
					if (yi == j || yi == j + 1) trials_i = c(trials_i, j)
				}
				if (length(trials_i) > 0) {
					X_stack_list[[i]] = matrix(rep(Xi, each = length(trials_i)), nrow = length(trials_i))
				}
			}
			X_stack = do.call(rbind, X_stack_list)

			mod = clogit_helper(expanded$y, as.data.frame(X_stack), expanded$w, expanded$strata)
			if (is.null(mod) || !is.finite(mod$b[1])) return(NA_real_)
			as.numeric(mod$b[1])
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
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
				private$cache_nonestimable_estimate("ordinal_multi_cond_adj_cat_logit_too_few_categories")
				return(invisible(NULL))
			}

			n_alpha = K - 1L
			n = private$n
			num_strata = max(strata_ids)
			X = private$get_X()

			# Data expansion for adjacent category (multivariate)
			expanded = expand_adjacent_category_data_cpp(as.integer(y_ord),
				as.integer(private$w), as.integer(strata_ids), as.integer(K))
			
			# Rebuild X_stack manually to match trial inclusion
			X_stack_list = list()
			for (i in 1:n) {
				yi = y_ord[i]
				Xi = X[i, ]
				# Trials j where subject i contributes: yi is j or j+1
				trials_i = integer()
				for (j in 1 : n_alpha) {
					if (yi == j || yi == j + 1) trials_i = c(trials_i, j)
				}
				if (length(trials_i) > 0) {
					X_stack_list[[i]] = matrix(rep(Xi, each = length(trials_i)),
						nrow = length(trials_i))
				}
			}
			X_stack = do.call(rbind, X_stack_list)
			
			X_full = cbind(treatment = expanded$w, X_stack)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					clogit_helper(
						expanded$y,
						as.data.frame(X_fit[, -1, drop = FALSE]),
						X_fit[, 1],
						expanded$strata
					)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) &&
						is.finite(mod$b[1]) &&
						is.finite(mod$ssq_b_j) &&
						mod$ssq_b_j > 0
				}
			)
			mod = attempt$fit
			
			if (!is.null(mod)){
				private$best_Xmm_colnames = setdiff(colnames(attempt$X_fit), "treatment")
			}

			if (is.null(mod) || !is.finite(mod$b[1]) || !is.finite(mod$ssq_b_j) || mod$ssq_b_j <= 0){
				if (private$harden && length(unique(expanded$strata)) > 0){
					private$cached_values$beta_hat_T   = 0
					private$cache_nonestimable_se("ordinal_multi_cond_adj_cat_logit_standard_error_unavailable")
				} else {
					private$cache_nonestimable_estimate("ordinal_multi_cond_adj_cat_logit_fit_unavailable")
				}
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = as.numeric(mod$b[1])
			private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z         = TRUE
		}
	)
)
