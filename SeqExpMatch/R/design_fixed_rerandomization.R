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
		#' @param obj_val_cutoff 	The maximum allowable objective value.
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
				obj_val_cutoff = NULL,
				objective = "mahal_dist",
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$obj_val_cutoff = obj_val_cutoff
			private$objective = objective
			private$uses_covariates = TRUE
		},


		#' @description
		#' Draw multiple treatment assignment vectors.
		#'
		#' @param r 	The number of designs to draw.
		#'
		#' @return 		A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				n_T = round(n * private$prob_T)
				return(replicate(r, sample(c(rep(1, n_T), rep(0, n - n_T)))))
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
					private$generate_one_rerandomized_w()
				}, mc.cores = private$num_cores)
				return(do.call(cbind, w_list))
			} else {
				w_mat = matrix(NA_real_, nrow = n, ncol = r)
				for (j in 1:r){
					w_mat[, j] = private$generate_one_rerandomized_w()
				}
				return(w_mat)
			}
		}
	),

	private = list(
		obj_val_cutoff = NULL,
		objective = NULL,
		S_inv = NULL,

		generate_one_rerandomized_w = function(){
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				n_T = round(n * private$prob_T)
				return(sample(c(rep(1, n_T), rep(0, n - n_T))))
			}

			# Ensure S_inv is available
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
			max_attempts = 10000 
			
			repeat {
				attempts = attempts + 1
				if (private$prob_T == 0.5){
					w_cand = sample(c(rep(1, n_T), rep(0, n - n_T)))
				} else {
					w_cand = rbinom(n, 1, private$prob_T)
				}
				
				obj_val = private$compute_obj(X, w_cand)
				
				if (is.null(private$obj_val_cutoff) || obj_val <= private$obj_val_cutoff || attempts >= max_attempts){
					if (private$verbose && attempts >= max_attempts){
						warning("Rerandomization reached max attempts without finding a design below cutoff.")
					}
					return(w_cand)
				}
			}
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
