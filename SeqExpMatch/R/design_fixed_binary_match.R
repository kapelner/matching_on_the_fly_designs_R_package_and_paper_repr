#' A Binary Match Fixed Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a fixed binary match experimental design.
#' This design pairs subjects based on covariate distances and randomizes within pairs.
#'
#' @export
FixedDesignBinaryMatch = R6::R6Class("FixedDesignBinaryMatch",
	inherit = FixedDesign,
	public = list(
		#' @description
		#' Initialize a binary match fixed experimental design
		#'
		#' @param response_type 	The data type of response values.
		#' @param prob_T	The probability of the treatment assignment. Must be 0.5.
		#' @param mahal_match 	Match using Mahalanobis distance. Default is \code{FALSE} (Euclidean).
		#' @param include_is_missing_as_a_new_feature	Flag for missingness indicators.
		#' @param n			The sample size.
		#' @param num_cores	The number of CPU cores.
		#' @param verbose	Flag for verbosity.
		#'
		#' @return 			A new `FixedDesignBinaryMatch` object
		#'
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				mahal_match = FALSE,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			if (prob_T != 0.5){
				stop("Binary match designs only support even treatment allocation (prob_T = 0.5)")
			}
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			private$mahal_match = mahal_match
			private$uses_covariates = TRUE
		},

		redraw_w_according_to_design = function(){
			private$ensure_pairs_computed()
			
			n = self$get_n()
			if (is.null(private$indices_pairs)){
				# Fallback if no covariates: just BCRD
				private$w[1:n] = sample(c(rep(1, n / 2), rep(0, n / 2)))
				return()
			}

			private$w[1:n] = private$generate_one_binary_match_w()
		},

		draw_ws_according_to_design = function(r = 100){
			private$ensure_pairs_computed()
			n = self$get_n()
			if (is.null(private$indices_pairs)){
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			# This is very fast so parallelization is likely not needed, but good for large r
			if (private$num_cores > 1 && r > 1000 && requireNamespace("pbmcapply", quietly = TRUE)){
				w_list = pbmcapply::pbmclapply(1:r, function(i) {
					private$generate_one_binary_match_w()
				}, mc.cores = private$num_cores)
				return(do.call(cbind, w_list))
			} else {
				w_mat = matrix(NA_real_, nrow = n, ncol = r)
				for (j in 1:r){
					w_mat[, j] = private$generate_one_binary_match_w()
				}
				return(w_mat)
			}
		}
	),

	private = list(
		mahal_match = NULL,
		indices_pairs = NULL,

		ensure_pairs_computed = function(){
			n = self$get_n()
			if (is.null(private$indices_pairs) && !is.null(private$X) && ncol(private$X) > 0){
				X = private$X[1:n, , drop = FALSE]
				if (private$mahal_match){
					S_X = var(X)
					if (abs(det(S_X)) < 1e-10){
						S_X = S_X + diag(1e-6, ncol(X))
					}
					S_X_inv = solve(S_X)
					D = matrix(0, n, n)
					for (i in 1 : (n - 1)){
						for (j in (i + 1) : n){
							xdiff = X[i, ] - X[j, ]
							D[i, j] = D[j, i] = xdiff %*% (S_X_inv %*% xdiff)
						}
					}
				} else {
					D = as.matrix(dist(X))^2
				}
				
				diag(D) = .Machine$double.xmax
				dist_obj = nbpMatching::distancematrix(D)
				match_obj = nbpMatching::nonbimatch(dist_obj)
				matches = match_obj$matches
				private$indices_pairs = matrix(as.integer(as.matrix(matches[, c("Group1.Row", "Group2.Row")])), ncol = 2)
			}
		},

		generate_one_binary_match_w = function(){
			n = self$get_n()
			new_w = rep(NA_real_, n)
			num_pairs = nrow(private$indices_pairs)
			for (i in 1 : num_pairs){
				pair = private$indices_pairs[i, ]
				if (runif(1) < 0.5){
					new_w[pair[1]] = 1
					new_w[pair[2]] = 0
				} else {
					new_w[pair[1]] = 0
					new_w[pair[2]] = 1
				}
			}
			new_w
		}
	)
)
