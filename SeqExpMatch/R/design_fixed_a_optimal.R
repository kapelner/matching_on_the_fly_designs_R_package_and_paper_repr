#' An A-optimal Search Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed A-optimal experimental design.
#' This design searches for an allocation that minimizes the trace of the inverse information matrix
#' (equivalent to minimizing the average variance of the parameter estimates).
#'
#' @export
FixedDesignAOptimal = R6::R6Class("FixedDesignAOptimal",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize an A-optimal search fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment. Must be 0.5 for exchange search.
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignAOptimal` object
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
				stop("A-optimal exchange search currently only supports even treatment allocation (prob_T = 0.5)")
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

			if (private$num_cores > 1 && requireNamespace("pbmcapply", quietly = TRUE)){
				w_list = pbmcapply::pbmclapply(1:r, function(i) {
					private$run_one_a_optimal_search()
				}, mc.cores = private$num_cores)
				return(do.call(cbind, w_list))
			} else {
				w_mat = matrix(NA_real_, nrow = n, ncol = r)
				for (j in 1:r){
					w_mat[, j] = private$run_one_a_optimal_search()
				}
				return(w_mat)
			}
		}
	),

	private = list(
		run_one_a_optimal_search = function(){
			n = self$get_n()
			X = private$X[1:n, , drop = FALSE]
			# Augmented matrix including intercept
			Z0 = cbind(1, X)
			
			# Initial BCRD
			w_curr = sample(c(rep(1, n / 2), rep(0, n / 2)))
			
			compute_obj = function(w){
				Z = cbind(Z0, w)
				# Tr( (Z'Z)^-1 )
				# We use pseudo-inverse if singular
				ZtZ = t(Z) %*% Z
				inv = tryCatch(solve(ZtZ), error = function(e) MASS::ginv(ZtZ))
				sum(diag(inv))
			}
			
			obj_curr = compute_obj(w_curr)
			
			repeat {
				best_obj = obj_curr
				best_switch = NULL
				
				t_idxs = which(w_curr == 1)
				c_idxs = which(w_curr == 0)
				
				improved = FALSE
				# A-optimal update is more expensive than D-optimal because we don't have
				# a simple quadratic form for the trace of the inverse.
				# However, we can use the Woodbury identity if needed.
				# For n small, we'll just recompute.
				for (i in t_idxs){
					for (j in c_idxs){
						w_cand = w_curr
						w_cand[i] = 0
						w_cand[j] = 1
						
						obj_cand = compute_obj(w_cand)
						
						if (obj_cand < best_obj - 1e-10){
							best_obj = obj_cand
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
