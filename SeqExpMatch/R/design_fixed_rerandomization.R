#' A Rerandomization Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed rerandomization
#' experimental design.
#' This design generates random allocations and only accepts those that meet a
#' covariate balance criterion.
#' For balanced designs (prob_T = 0.5, even n) uses the \pkg{GreedyExperimentalDesign} package.
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
				
				verbose = FALSE
			) {
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
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
			assertCount(r, positive = TRUE)
			self$assert_all_subjects_arrived()
			n = self$get_n()
			if (is.null(private$X) || ncol(private$X) == 0){
				n_T = round(n * private$prob_T)
				return(replicate(r, sample(c(rep(1, n_T), rep(0, n - n_T)))))
			}

			# Use GED for the balanced even-n case
			use_ged = private$prob_T == 0.5 && n %% 2 == 0 &&
					  requireNamespace("GreedyExperimentalDesign", quietly = TRUE) &&
					  requireNamespace("rJava", quietly = TRUE)

			if (use_ged){
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
				X = private$X[1:n, , drop = FALSE]
				cutoff = if (is.null(private$obj_val_cutoff)) Inf else private$obj_val_cutoff
				search_obj = GreedyExperimentalDesign::initRerandomizationExperimentalDesignObject(
					X                      = X,
					obj_val_cutoff_to_include = cutoff,
					max_designs            = r,
					objective              = private$objective,
					wait                   = TRUE,
					start                  = TRUE,
					num_cores              = self$num_cores,
					verbose                = private$verbose
				)
				res = GreedyExperimentalDesign::resultsRerandomizationSearch(search_obj, include_assignments = TRUE, form = "one_zero")
				w_mat = res$ending_indicTs
				storage.mode(w_mat) = "numeric"
				# Recycle if fewer were found than requested (tight cutoff)
				if (ncol(w_mat) < r){
					w_mat = w_mat[, rep(seq_len(ncol(w_mat)), length.out = r), drop = FALSE]
				}
				return(w_mat[, seq_len(r), drop = FALSE])
			}

			# Fallback pure-R for unbalanced or odd-n cases
			if (private$objective == "mahal_dist" && is.null(private$S_inv)){
				X = private$X[1:n, , drop = FALSE]
				S = var(X)
				if (abs(det(S)) < 1e-10){
					S = S + diag(1e-6, ncol(X))
				}
				private$S_inv = solve(S)
			}

			w_mat = matrix(NA_real_, nrow = n, ncol = r)
			for (j in seq_len(r)){
				w_mat[, j] = private$generate_one_rerandomized_w()
			}
			w_mat
		}
	),

	private = list(
		obj_val_cutoff = NULL,
		objective = NULL,
		S_inv = NULL,

		generate_one_rerandomized_w = function(){
			n = self$get_n()
			X = private$X[1:n, , drop = FALSE]
			n_T = round(n * private$prob_T)

			repeat {
				if (private$prob_T == 0.5){
					w_cand = sample(c(rep(1, n_T), rep(0, n - n_T)))
				} else {
					w_cand = rbinom(n, 1, private$prob_T)
				}
				obj_val = private$compute_obj(X, w_cand)
				if (is.null(private$obj_val_cutoff) || obj_val <= private$obj_val_cutoff){
					return(w_cand)
				}
			}
		},

		compute_obj = function(X, w){
			if (private$objective == "mahal_dist"){
				diff = colMeans(X[w == 1, , drop = FALSE]) - colMeans(X[w == 0, , drop = FALSE])
				return(as.numeric(diff %*% private$S_inv %*% diff))
			} else if (private$objective == "abs_sum_diff"){
				diff = colMeans(X[w == 1, , drop = FALSE]) - colMeans(X[w == 0, , drop = FALSE])
				return(sum(abs(diff)))
			}
			stop("Unsupported objective")
		}
	)
)
