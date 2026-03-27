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
		#' @param response_type "continuous", "incidence", "proportion", "count", "survival", or "ordinal".
		#' @param prob_T Probability of treatment assignment (default 0.5).
		#' @param include_is_missing_as_a_new_feature Flag for missingness indicators.
		#' @param n Sample size (if fixed).
		#' @param num_cores Number of CPU cores.
		#' @param verbose Flag for verbosity.
		#'
		#' @return A new `FixedDesignAOptimal` object
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


		#' @description
		#' Draw treatment assignments according to the A-optimal search fixed design.
		#'
		#' @param r Number of designs to draw.
		#'
		#' @return A matrix of size n x r.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			if (is.null(private$P)){
				X = private$X[1:n, , drop = FALSE]
				Z0 = cbind(1, X)
				
				# P = Z0 (Z0'Z0)^-1 Z0'
				Z0_qr = qr(Z0)
				Q = qr.Q(Z0_qr)
				private$P = Q %*% t(Q)
				
				# H = Z0 (Z0'Z0)^-2 Z0'
				# (Z0'Z0)^-1 = R^-1 Q' Q (R')^-1 = R^-1 (R')^-1
				# (Z0'Z0)^-2 = R^-1 (R')^-1 R^-1 (R')^-1
				# Actually, let M = (Z0'Z0)^-1
				# H = Z0 M^2 Z0'
				R = qr.R(Z0_qr)
				# M = solve(t(R) %*% R) # If R is small (p+1 x p+1)
				M = tryCatch(solve(t(R) %*% R), error = function(e) MASS::ginv(t(R) %*% R))
				H_kernel = M %*% M
				private$H = Z0 %*% H_kernel %*% t(Z0)
			}

			# Use C++ speedup
			res = a_optimal_search_cpp(private$P, private$H, as.integer(r), as.integer(round(n * private$prob_T)))
			w_mat = res
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),

	private = list(
		P = NULL, # Projection matrix
		H = NULL  # Trace-inverse kernel matrix
	)
)
