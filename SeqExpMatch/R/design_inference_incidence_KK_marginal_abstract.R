# Abstract class for all-subject marginal incidence inference in KK designs
#
# @keywords internal
DesignInferenceAbstractKKMarginalIncid = R6::R6Class("DesignInferenceAbstractKKMarginalIncid",
	inherit = DesignInferenceKKPassThrough,
	public = list(

		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		}
	),

	private = list(
		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

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

		get_cluster_ids = function(){
			cluster_id = private$cached_values$cluster_id
			if (!is.null(cluster_id)) return(cluster_id)

			m_vec = private$m
			if (is.null(m_vec)){
				m_vec = rep(0L, private$n)
			}
			m_vec = as.integer(m_vec)
			m_vec[is.na(m_vec)] = 0L
			cluster_id = compute_kk_cluster_ids_cpp(m_vec)
			private$cached_values$cluster_id = cluster_id
			cluster_id
		}
	)
)
