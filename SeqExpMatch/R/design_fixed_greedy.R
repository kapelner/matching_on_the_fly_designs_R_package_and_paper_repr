#' A Greedy Search Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed greedy experimental design.
#' This design starts with a random allocation and greedily switches pairs to optimize covariate balance.
#'
#' @export
FixedDesignGreedy = R6::R6Class("FixedDesignGreedy",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a greedy search fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment. Must be 0.5.
		#' @param objective 	The objective function to use. Default is "mahal_dist".
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignGreedy` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				objective = "mahal_dist",
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			if (prob_T != 0.5){
				stop("Greedy designs currently only support even treatment allocation (prob_T = 0.5)")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$objective = objective
			private$uses_covariates = TRUE
		},

		draw_ws_according_to_design = function(r = 100){
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			# Precompute S_inv if needed
			if (private$objective == "mahal_dist" && is.null(private$S_inv)){
				X = private$X[1:n, , drop = FALSE]
				S = var(X)
				if (abs(det(S)) < 1e-10){
					S = S + diag(1e-6, ncol(X))
				}
				private$S_inv = solve(S)
			}

			if (private$num_cores > 1 && requireNamespace("pbmcapply", quietly = TRUE)){
				w_list = pbmcapply::pbmclapply(1:r, function(i) {
					private$run_one_greedy_search()
				}, mc.cores = private$num_cores)
				return(do.call(cbind, w_list))
			} else {
				w_mat = matrix(NA_real_, nrow = n, ncol = r)
				for (j in 1:r){
					w_mat[, j] = private$run_one_greedy_search()
				}
				return(w_mat)
			}
		}
	),

	private = list(
		objective = NULL,
		S_inv = NULL,

		run_one_greedy_search = function(){
			n = self$get_n()
			X = private$X[1:n, , drop = FALSE]
			
			# Ensure S_inv is available
			if (private$objective == "mahal_dist" && is.null(private$S_inv)){
				S = var(X)
				if (abs(det(S)) < 1e-10){
					S = S + diag(1e-6, ncol(X))
				}
				private$S_inv = solve(S)
			}

			# Initial BCRD
			w_curr = sample(c(rep(1, n / 2), rep(0, n / 2)))
			obj_curr = private$compute_obj(X, w_curr)
			
			repeat {
				best_obj = obj_curr
				best_switch = NULL
				
				t_idxs = which(w_curr == 1)
				c_idxs = which(w_curr == 0)
				
				improved = FALSE
				for (i in t_idxs){
					for (j in c_idxs){
						# Try switching
						w_cand = w_curr
						w_cand[i] = 0
						w_cand[j] = 1
						
						obj_cand = private$compute_obj(X, w_cand)
						
						if (obj_cand < best_obj){
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
		},

		compute_obj = function(X, w){
			if (private$objective == "mahal_dist"){
				x_T_bar = colMeans(X[w == 1, , drop = FALSE])
				x_C_bar = colMeans(X[w == 0, , drop = FALSE])
				diff = x_T_bar - x_C_bar
				return(as.numeric(diff %*% private$S_inv %*% diff))
			} else if (private$objective == "abs_sum_diff"){
				x_T_bar = colMeans(X[w == 1, , drop = FALSE])
				x_C_bar = colMeans(X[w == 0, , drop = FALSE])
				return(sum(abs(x_T_bar - x_C_bar)))
			}
			stop("Unsupported objective")
		}
	)
)
