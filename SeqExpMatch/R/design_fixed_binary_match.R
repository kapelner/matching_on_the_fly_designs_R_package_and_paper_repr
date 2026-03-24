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
			if (is.null(private$m)){
				# Fallback if no covariates: just BCRD
				private$w[1:n] = sample(c(rep(1, n / 2), rep(0, n / 2)))
				return()
			}

			# Use C++ for fast generation
			res = generate_permutations_kk_cpp(as.integer(private$m), 1, as.numeric(private$prob_T))
			private$w[1:n] = res$w_mat[, 1]
		},

		draw_ws_according_to_design = function(r = 100){
			private$ensure_pairs_computed()
			n = self$get_n()
			if (is.null(private$m)){
				# Fallback BCRD
				return(replicate(r, sample(c(rep(1, n/2), rep(0, n/2)))))
			}

			# Use C++ for fast generation of many replicates
			res = generate_permutations_kk_cpp(as.integer(private$m), as.integer(r), as.numeric(private$prob_T))
			return(res$w_mat)
		}
	),

	private = list(
		mahal_match = NULL,

		ensure_pairs_computed = function(){
			n = self$get_n()
			if (is.null(private$m) && !is.null(private$X) && ncol(private$X) > 0){
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
				
				# API call to nbpMatching::nonbimatch
				dist_obj = nbpMatching::distancematrix(D)
				match_obj = nbpMatching::nonbimatch(dist_obj)
				matches = match_obj$matches
				
				# Convert to m vector (pair IDs)
				m_vec = rep(NA_integer_, n)
				for (i in 1 : nrow(matches)){
					m_vec[as.integer(matches$Group1.Row[i])] = i
					m_vec[as.integer(matches$Group2.Row[i])] = i
				}
				private$m = m_vec
			}
		}
	)
)
