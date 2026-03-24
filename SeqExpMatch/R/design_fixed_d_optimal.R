#' A D-optimal Search Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed D-optimal experimental design.
#' This design searches for an allocation that maximizes the determinant of the information matrix
#' (equivalent to minimizing the variance of the parameter estimates).
#'
#' @export
FixedDesignDOptimal = R6::R6Class("FixedDesignDOptimal",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a D-optimal search fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment. Must be 0.5 for exchange search.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignDOptimal` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			if (prob_T != 0.5){
				stop("D-optimal exchange search currently only supports even treatment allocation (prob_T = 0.5)")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$uses_covariates = TRUE
		},

		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			# Precompute projection matrix P = Z0 (Z0' Z0)^-1 Z0'
			# where Z0 = [1, X]
			if (is.null(private$P)){
				X = private$X[1:n, , drop = FALSE]
				Z0 = cbind(1, X)
				# Use pseudo-inverse if singular
				# P = Z0 %*% solve(t(Z0) %*% Z0) %*% t(Z0)
				# But more stable:
				Z0_qr = qr(Z0)
				Q = qr.Q(Z0_qr)
				private$P = Q %*% t(Q)
			}

			if (private$num_cores > 1 && requireNamespace("pbmcapply", quietly = TRUE)){
				w_list = pbmcapply::pbmclapply(1:r, function(i) {
					private$run_one_d_optimal_search()
				}, mc.cores = private$num_cores)
				return(do.call(cbind, w_list))
			} else {
				w_mat = matrix(NA_real_, nrow = n, ncol = r)
				for (j in 1:r){
					w_mat[, j] = private$run_one_d_optimal_search()
				}
				return(w_mat)
			}
		}
	),

	private = list(
		P = NULL, # Projection matrix

		run_one_d_optimal_search = function(){
			n = self$get_n()
			# Initial BCRD
			w_curr = sample(c(rep(1, n / 2), rep(0, n / 2)))
			# Convert 0/1 to -1/1 for the math to be symmetric if needed, 
			# but for det(Z'Z) where Z=[1, X, w], using 0/1 is fine.
			# Objective: minimize w' P w
			obj_curr = as.numeric(t(w_curr) %*% private$P %*% w_curr)
			
			repeat {
				best_obj = obj_curr
				best_switch = NULL
				
				t_idxs = which(w_curr == 1)
				c_idxs = which(w_curr == 0)
				
				improved = FALSE
				for (i in t_idxs){
					for (j in c_idxs){
						# Delta in objective if switching i (1->0) and j (0->1)
						# w_new = w_curr - ei + ej
						# w_new' P w_new = (w' - ei' + ej') P (w - ei + ej)
						# = w'Pw - 2w'P ei + 2w'P ej + ei'P ei + ej'P ej - 2ei'P ej
						# = w'Pw - 2 P_row_i %*% w + 2 P_row_j %*% w + P_ii + P_jj - 2 P_ij
						
						# Precompute P %*% w for the whole loop to be O(n^2) instead of O(n^3)
						# Actually, compute it once per outer iteration.
						if (is.null(best_switch)){
							Pw = private$P %*% w_curr
						}
						
						delta = - 2 * Pw[i] + 2 * Pw[j] + private$P[i, i] + private$P[j, j] - 2 * private$P[i, j]
						
						if (delta < -1e-10 && (obj_curr + delta) < best_obj){
							best_obj = obj_curr + delta
							best_switch = c(i, j)
							improved = TRUE
						}
					}
				}
				
				if (improved){
					w_curr[best_switch[1]] = 0
					w_curr[best_switch[2]] = 1
					obj_curr = best_obj
				} else {
					break
				}
			}
			w_curr
		}
	)
)
