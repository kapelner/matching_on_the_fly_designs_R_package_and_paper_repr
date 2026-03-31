#' A class that provides for relevant methods when the designs are KK matching-on-the-fly
#'
#' An abstract class
#'
#' @keywords internal
InferenceKKPassThrough = R6::R6Class("InferenceKKPassThrough",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize
		#' @param des_obj         A DesignSeqOneByOne object whose entire n subjects are assigned
		#'   and response y is recorded within.
		#' @param num_cores                       The number of CPU cores to use to parallelize
		#'   the sampling during randomization-based inference
		#' and bootstrap resampling. The default is 1 for serial computation. For simple
		#' estimators (e.g. mean difference
		#' and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For
		#' complex models (e.g. GLMs),
		#' parallelization falls back to R's \code{parallel::mclapply} which incurs
		#' session-forking overhead.
		#' @param verbose                 A flag indicating whether messages should be displayed
		#'   to the user. Default is \code{TRUE}
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(des_obj, num_cores, verbose)
				if (private$has_match_structure){
					# For fixed binary matching, we need to ensure pairs are computed first
					if (is(des_obj, "FixedDesignBinaryMatch")){
						des_obj$.__enclos_env__$private$ensure_bms_computed()
					}
					private$m = des_obj$.__enclos_env__$private$m
					private$compute_basic_match_data()
				}
		},


		#' @description
		#' Creates the boostrap distribution of the estimate for the treatment effect
		#'
		#' @param B						Number of bootstrap samples. The default is 501.
		#'
		#' @return A vector of length \code{B} with the bootstrap values of the estimates of the
		#'   treatment effect
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneKK14$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#'
		#' seq_des_inf = InferenceContinMultOLSKK$new(seq_des)
		#' beta_hat_T_bs = seq_des_inf$approximate_bootstrap_distribution_beta_hat_T(B = 5)
		#' beta_hat_T_bs
		#' }
		#'
		#' @param show_progress Description for show_progress
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE){
				if (!private$has_match_structure){
					super$approximate_bootstrap_distribution_beta_hat_T(B, show_progress)
				} else {
				assertCount(B, positive = TRUE)

				if (is.null(private$cached_values$KKstats)){
					private$compute_basic_match_data()
				}

				n = private$n
				y = private$y
				dead = private$dead
				w = private$w
				X = private$get_X()
				m_vec = private$m
				if (is.null(m_vec)){
					m_vec = rep(NA_integer_, n)
				}
				m_vec[is.na(m_vec)] = 0
				m = private$cached_values$KKstats$m

				# Identify reservoir and matched observations
				i_reservoir = which(m_vec == 0)
				n_reservoir = length(i_reservoir)

				# Check if subclass provides a C++ OpenMP dispatcher to bypass the slow R loop
				if (private$has_private_method("compute_fast_bootstrap_distr")) {
					fast_distr = private$compute_fast_bootstrap_distr(B, i_reservoir, n_reservoir, m, y, w, m_vec)
					if (!is.null(fast_distr)) {
						return(fast_distr)
					}
				}

				# Pure R KK bootstrap implementation
				if (private$num_cores == 1) {
					pb = NULL
					if (private$verbose) {
						pb = utils::txtProgressBar(min = 0, max = B, style = 3)
						on.exit(close(pb), add = TRUE)
					}
					beta_hat_T_bs = rep(NA_real_, B)

					has_res_stat_serial = private$has_private_method("compute_reservoir_and_match_statistics")
					for (b in 1:B) {
						i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
						if (m > 0) {
							pairs_to_include = sample(seq_len(m), m, replace = TRUE)
							i_matched_b = integer(0); m_vec_b_matched = integer(0)
							for (new_pair_id in seq_len(m)) {
								orig_pid = pairs_to_include[new_pair_id]
								pair_idx = which(m_vec == orig_pid)
								i_matched_b = c(i_matched_b, pair_idx)
								m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
							}
						} else {
							i_matched_b = integer(0); m_vec_b_matched = integer(0)
						}
						i_b = c(i_reservoir_b, i_matched_b)
						m_vec_b = c(rep(0L, n_reservoir), m_vec_b_matched)
						boot_inf_obj = self$duplicate()
						boot_inf_obj$.__enclos_env__$private$num_cores = 1L
						boot_inf_obj$.__enclos_env__$private$y = y[i_b]
						boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
						boot_inf_obj$.__enclos_env__$private$w = w[i_b]
						boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
						boot_inf_obj$.__enclos_env__$private$cached_values = list()
						beta_hat_T_bs[b] = tryCatch({
							boot_inf_obj$.__enclos_env__$private$m = m_vec_b
							boot_inf_obj$.__enclos_env__$private$compute_basic_match_data()
							if (has_res_stat_serial) boot_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
							boot_inf_obj$compute_treatment_estimate(estimate_only = TRUE)
						}, error = function(e) { NA_real_ })
						if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
					}
					return(beta_hat_T_bs)
				} else {
					# Parallel bootstrap execution for KK designs
					cores_to_use = private$num_cores

					# Warm-up guard: run 1 iteration at 1 thread to estimate per-iteration cost.
					# Fork overhead on large R sessions can be 0.1-1s; only parallelize if
					# computation clearly outweighs that cost.
					if (cores_to_use > 1L) {
						i_res_w = sample(i_reservoir, n_reservoir, replace = TRUE)
						if (m > 0) {
							pairs_w = sample(1:m, m, replace = TRUE)
							i_mat_w = integer(0); m_vm_w = integer(0)
							for (np in 1:m) {
								op = pairs_w[np]; pi = which(m_vec == op)
								i_mat_w = c(i_mat_w, pi); m_vm_w = c(m_vm_w, np, np)
							}
						} else {
							i_mat_w = integer(0); m_vm_w = integer(0)
						}
						i_b_w = c(i_res_w, i_mat_w)
						m_vec_b_w = c(rep(0, n_reservoir), m_vm_w)
						t_warmup_kk = system.time({
							boot_w = self$duplicate()
							boot_w$.__enclos_env__$private$num_cores = 1L
							boot_w$.__enclos_env__$private$y = y[i_b_w]
							boot_w$.__enclos_env__$private$dead = dead[i_b_w]
							boot_w$.__enclos_env__$private$w = w[i_b_w]
							boot_w$.__enclos_env__$private$X = X[i_b_w, , drop = FALSE]
							boot_w$.__enclos_env__$private$cached_values = list()
							boot_w$.__enclos_env__$private$m = m_vec_b_w
							boot_w$.__enclos_env__$private$compute_basic_match_data()
							if (private$object_has_private_method(boot_w, "compute_reservoir_and_match_statistics"))
								boot_w$.__enclos_env__$private$compute_reservoir_and_match_statistics()
							tryCatch(boot_w$compute_treatment_estimate(estimate_only = TRUE), error = function(e) NA_real_)
						})[[3]]
						fork_overhead_estimate = if (isTRUE(private$make_fork_cluster)) 0.01 else 0.5
						if (!(t_warmup_kk * B > fork_overhead_estimate * cores_to_use))
							cores_to_use = 1L
					}

					kk_template = self$duplicate()
					has_res_stat = private$object_has_private_method(kk_template, "compute_reservoir_and_match_statistics")

					# Internal task for bootstrap execution
					kk_task = function(b) {
						try(set_package_threads(1L), silent = TRUE)
						i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
						if (m > 0) {
							pairs_to_include = sample(seq_len(m), m, replace = TRUE)
							i_matched_b = integer(0); m_vec_b_matched = integer(0)
							for (new_pair_id in seq_len(m)) {
								orig_pid = pairs_to_include[new_pair_id]
								pair_idx = which(m_vec == orig_pid)
								i_matched_b = c(i_matched_b, pair_idx)
								m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
							}
						} else {
							i_matched_b = integer(0); m_vec_b_matched = integer(0)
						}
						i_b = c(i_reservoir_b, i_matched_b)
						m_vec_b = c(rep(0L, n_reservoir), m_vec_b_matched)
						boot_inf_obj = kk_template$duplicate()
						boot_inf_obj$.__enclos_env__$private$num_cores = 1L
						boot_inf_obj$.__enclos_env__$private$y = y[i_b]
						boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
						boot_inf_obj$.__enclos_env__$private$w = w[i_b]
						boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
						boot_inf_obj$.__enclos_env__$private$cached_values = list()
						tryCatch({
							boot_inf_obj$.__enclos_env__$private$m = m_vec_b
							boot_inf_obj$.__enclos_env__$private$compute_basic_match_data()
							if (has_res_stat) boot_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
							boot_inf_obj$compute_treatment_estimate(estimate_only = TRUE)
						}, error = function(e) NA_real_)
					}

					if (isTRUE(private$make_fork_cluster) && cores_to_use > 1L) {
						cl = private$get_or_create_fork_cluster()
						beta_hat_T_bs = unlist(parallel::parLapply(cl, 1:B, kk_task))
					} else {
						has_res_stat_local = has_res_stat
						beta_hat_T_bs = unlist(private$par_lapply(1:B, function(b) {
							set_package_threads(1L)
							i_reservoir_b = sample(i_reservoir, n_reservoir, replace = TRUE)
							if (m > 0) {
								pairs_to_include = sample(1:m, m, replace = TRUE)
								i_matched_b = integer(0); m_vec_b_matched = integer(0)
								for (new_pair_id in 1:m) {
									orig_pid = pairs_to_include[new_pair_id]
									pair_idx = which(m_vec == orig_pid)
									i_matched_b = c(i_matched_b, pair_idx)
									m_vec_b_matched = c(m_vec_b_matched, new_pair_id, new_pair_id)
								}
							} else {
								i_matched_b = integer(0); m_vec_b_matched = integer(0)
							}
							i_b = c(i_reservoir_b, i_matched_b)
							m_vec_b = c(rep(0L, n_reservoir), m_vec_b_matched)
							boot_inf_obj = kk_template$duplicate()
							boot_inf_obj$.__enclos_env__$private$num_cores = 1L
							boot_inf_obj$.__enclos_env__$private$y = y[i_b]
							boot_inf_obj$.__enclos_env__$private$dead = dead[i_b]
							boot_inf_obj$.__enclos_env__$private$w = w[i_b]
							boot_inf_obj$.__enclos_env__$private$X = X[i_b, , drop = FALSE]
							boot_inf_obj$.__enclos_env__$private$cached_values = list()
							tryCatch({
								boot_inf_obj$.__enclos_env__$private$m = m_vec_b
								boot_inf_obj$.__enclos_env__$private$compute_basic_match_data()
								if (has_res_stat_local) boot_inf_obj$.__enclos_env__$private$compute_reservoir_and_match_statistics()
								boot_inf_obj$compute_treatment_estimate(estimate_only = TRUE)
							}, error = function(e) NA_real_)
						}, n_cores = cores_to_use, show_progress = private$verbose))
					}
					return(beta_hat_T_bs)
				}
			}
		}
	),
	private = list(

		m = NULL,

		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			private$cached_values$KKstats = .compute_kk_basic_match_data(
				X = private$X,
				n = private$n,
				y = private$y,
				w = private$w,
				m_vec = private$m
			)
		},

		compute_concordant_and_discordant_match_statistics = function(){
			m = private$cached_values$KKstats$m
			y_matched_diffs = private$cached_values$KKstats$y_matched_diffs
			i_m_conc = which(y_matched_diffs == 0)
			i_m_disc = setdiff(seq_len(m), i_m_conc)
			private$cached_values$KKstats$i_m_conc = i_m_conc
			private$cached_values$KKstats$i_m_disc = i_m_disc
			private$cached_values$KKstats$n_m_conc = length(i_m_conc)
			private$cached_values$KKstats$n_m_disc = length(i_m_disc)
			private$cached_values$KKstats$y_matched_diffs_disc = y_matched_diffs[i_m_disc]
			private$cached_values$KKstats$X_matched_diffs_disc =
				private$cached_values$KKstats$X_matched_diffs[i_m_disc, drop = FALSE]
			i_conc = which(private$m %in% i_m_conc)
			private$cached_values$KKstats$X_conc = private$X[i_conc, drop = FALSE]
			private$cached_values$KKstats$y_conc = private$y[i_conc]
			private$cached_values$KKstats$w_conc = private$w[i_conc]
		},

		#not used now, but could be used for random effects models in the future
		compute_model_matrix_with_matching_dummies = function(){
			if (is.null(private$cached_values$data_frame_with_matching_dummies)){
				if (!is.null(private$m) & uniqueN(private$m) > 1){
					mm = model.matrix(~ 0 + factor(private$m))
					mm = mm[, 2 : (ncol(mm) - 1)]
				} else {
					mm = NULL
				}
				private$cached_values$data_frame_with_matching_dummies = cbind(data.frame(w = private$w), mm)
			}
			private$cached_values$data_frame_with_matching_dummies
		}
	)
)
