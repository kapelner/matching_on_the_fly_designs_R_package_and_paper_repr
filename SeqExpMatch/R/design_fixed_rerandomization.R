#' A Rerandomization Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed rerandomization experimental design.
#' This design generates random allocations and only accepts those that meet a covariate balance criterion.
#'
#' @export
FixedDesignRerandomization = R6::R6Class("FixedDesignRerandomization",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a rerandomization fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment.
		#' @param obj_val_cutoff_to_include 	The maximum allowable objective value.
		#' @param objective 	The objective function to use. Default is "mahal_dist".
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignRerandomization` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				obj_val_cutoff_to_include = NULL,
				objective = "mahal_dist",
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$obj_val_cutoff_to_include = obj_val_cutoff_to_include
			private$objective = objective
			private$uses_covariates = TRUE
		},

		redraw_w_according_to_design = function(){
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				# Fallback if no covariates: just BCRD
				n_T = round(n * private$prob_T)
				private$w[1:n] = sample(c(rep(1, n_T), rep(0, n - n_T)))
				return()
			}

			# Precompute inverse covariance for Mahalanobis if needed
			if (private$objective == "mahal_dist" && is.null(private$S_inv)){
				X = private$X[1:n, , drop = FALSE]
				S = var(X)
				if (abs(det(S)) < 1e-10){
					S = S + diag(1e-6, ncol(X))
				}
				private$S_inv = solve(S)
			}

			X = private$X[1:n, , drop = FALSE]
			n_T = round(n * private$prob_T)
			
			attempts = 0
			max_attempts = 10000 # Safety break
			
			repeat {
				attempts = attempts + 1
				# Generate candidate w
				if (private$prob_T == 0.5){
					# BCRD is better for rerandomization usually
					w_cand = sample(c(rep(1, n_T), rep(0, n - n_T)))
				} else {
					w_cand = rbinom(n, 1, private$prob_T)
				}
				
				# Compute objective
				obj_val = private$compute_obj(X, w_cand)
				
				if (is.null(private$obj_val_cutoff_to_include) || obj_val <= private$obj_val_cutoff_to_include || attempts >= max_attempts){
					private$w[1:n] = w_cand
					if (private$verbose && attempts >= max_attempts){
						warning("Rerandomization reached max attempts without finding a design below cutoff.")
					}
					break
				}
			}
		}
	),

	private = list(
		obj_val_cutoff_to_include = NULL,
		objective = NULL,
		S_inv = NULL,

		compute_obj = function(X, w){
			if (private$objective == "mahal_dist"){
				# Mahal dist = (x_T_bar - x_C_bar)' S^-1 (x_T_bar - x_C_bar)
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
