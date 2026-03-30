#' Azriel Blocked Incidence Inference
#'
#' Unadjusted blocked-design incidence inference using the simple mean-difference
#' point estimate with a block-stratified standard error.
#'
#' @export
InferenceIncidAzriel = R6::R6Class("InferenceIncidAzriel",
	lock_objects = FALSE,
	inherit = InferenceAllSimpleMeanDiff,
	public = list(
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "incidence")
			if (!private$design_is_supported_blocking(des_obj)) {
				stop(class(self)[1], " requires a blocking design.")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
			if (is.null(private$get_block_ids())) {
				stop(class(self)[1], " requires block identifiers from the completed design.")
			}
		}
	),

	private = list(
		design_is_supported_blocking = function(des_obj){
			is(des_obj, "FixedDesignBlocking") ||
				is(des_obj, "DesignSeqOneByOneSPBR") ||
				is(des_obj, "DesignSeqOneByOneRandomBlockSize")
		},

		get_block_ids = function(){
			m_vec = tryCatch(private$des_obj$get_m(), error = function(e) private$des_obj_priv_int$m)
			if (is.null(m_vec)) {
				strata_cols = private$des_obj_priv_int$strata_cols
				Xraw = private$des_obj_priv_int$Xraw
				if (!is.null(strata_cols) && length(strata_cols) > 0L && nrow(Xraw) == length(private$y)) {
					strata_keys = vapply(seq_len(nrow(Xraw)), function(i) {
						vals = vapply(strata_cols, function(col) {
							val = Xraw[i, ][[col]]
							if (is.na(val)) "NA" else as.character(val)
						}, character(1))
						paste(vals, collapse = "|")
					}, character(1))
					m_vec = match(strata_keys, unique(strata_keys))
				}
			}
			if (is.null(m_vec)) return(NULL)
			m_vec = as.integer(m_vec)
			if (length(m_vec) != length(private$y)) return(NULL)
			m_vec
		},

		get_standard_error = function(){
			if (!is.null(private$cached_values$azriel_s_beta_hat_T)) {
				return(private$cached_values$azriel_s_beta_hat_T)
			}

			block_ids = private$get_block_ids()
			if (is.null(block_ids)) {
				private$cached_values$azriel_s_beta_hat_T = NA_real_
				return(NA_real_)
			}

			keep = is.finite(private$y) & !is.na(private$w) & !is.na(block_ids)
			y = private$y[keep]
			w = private$w[keep]
			block_ids = block_ids[keep]
			n = length(y)
			if (n == 0L) {
				private$cached_values$azriel_s_beta_hat_T = NA_real_
				return(NA_real_)
			}

			var_est = 0
			for (block_id in unique(block_ids)) {
				i_block = which(block_ids == block_id)
				w_block = w[i_block]
				y_block = y[i_block]
				n_block = length(i_block)
				n_t = sum(w_block == 1L)
				n_c = sum(w_block == 0L)
				if (n_block == 0L || n_t == 0L || n_c == 0L) next

				s2_t = if (n_t > 1L) stats::var(y_block[w_block == 1L]) else 0
				s2_c = if (n_c > 1L) stats::var(y_block[w_block == 0L]) else 0
				weight = n_block / n
				var_est = var_est + weight^2 * (s2_t / n_t + s2_c / n_c)
			}

			se = if (is.finite(var_est) && var_est >= 0) sqrt(var_est) else NA_real_
			private$cached_values$azriel_s_beta_hat_T = se
			se
		},

		get_degrees_of_freedom = function(){
			NA_real_
		}
	)
)
