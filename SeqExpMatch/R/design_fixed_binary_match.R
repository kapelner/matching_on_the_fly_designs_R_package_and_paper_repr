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
			n = self$get_n()
			if (n %% 2 != 0){
				stop("Binary match designs require an even number of subjects")
			}
			
			if (is.null(private$X) || ncol(private$X) == 0){
				# Fallback if no covariates: just BCRD
				private$w[1:n] = sample(c(rep(1, n / 2), rep(0, n / 2)))
				return()
			}

			if (is.null(private$indices_pairs)){
				# Compute pairs once
				X = private$X[1:n, , drop = FALSE]
				if (private$mahal_match){
					S_X = var(X)
					# Add small ridge if singular
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
				
				# Ensure diagonal is infinity
				diag(D) = .Machine$double.xmax
				
				# Use nbpMatching
				# nbpMatching::nonbimatch expects a distance matrix object
				dist_obj = nbpMatching::distancematrix(D)
				match_obj = nbpMatching::nonbimatch(dist_obj)
				
				# Extract pairs
				matches = match_obj$matches
				private$indices_pairs = as.matrix(matches[, c("Group1.Row", "Group2.Row")])
				# Ensure 1-based indexing if not already (nbpMatching uses row names/indices)
				# Actually row names in matches are from the original D matrix
				private$indices_pairs = matrix(as.integer(private$indices_pairs), ncol = 2)
			}

			# Randomize within pairs
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
			private$w[1:n] = new_w
		}
	),

	private = list(
		mahal_match = NULL,
		indices_pairs = NULL
	)
)
