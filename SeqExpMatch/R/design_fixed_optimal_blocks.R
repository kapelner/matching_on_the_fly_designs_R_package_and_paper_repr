#' An Optimal-Blocks Fixed Design
#'
#' An R6 Class encapsulating the data and functionality for a fixed experimental
#' design that first partitions subjects into \code{K} covariate-homogeneous blocks
#' by solving a balanced clustering problem with \pkg{ompr} and glpk, then randomizes
#' treatment within those blocks.
#'
#' @export
FixedDesignOptimalBlocks = R6::R6Class("FixedDesignOptimalBlocks",
	inherit = FixedDesign,
	public = list(
			#' @description
			#' Initialize a fixed optimal-blocks design.
			#' @param K Number of blocks to form.
			#' @param dist Distance specification, either a function or one of
			#'   \code{"euclidean"}, \code{"sum_abs_diff"}, or \code{"mahal"}.
			#' @param response_type The response type for the design.
			#' @param prob_T Treatment assignment probability within each block.
			#' @param include_is_missing_as_a_new_feature Whether to include missingness indicators.
			#' @param n Planned sample size.
			#' @param verbose Whether to print progress messages.
			#' @return A new \code{FixedDesignOptimalBlocks} object.
			initialize = function(
					K,
				dist = "euclidean",
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				
				verbose = FALSE
			) {
			assertCount(K, positive = TRUE)
			if (!(is.function(dist) || (is.character(dist) && length(dist) == 1L))) {
				stop("dist must be a function or one of 'euclidean', 'sum_abs_diff', or 'mahal'.")
			}
			if (is.character(dist)) {
				assertChoice(dist, c("euclidean", "sum_abs_diff", "mahal"))
			}
				assert_optimal_blocks_libraries_installed("FixedDesignOptimalBlocks")
				super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, verbose)
				private$K = as.integer(K)
				private$dist_spec = dist
				private$uses_covariates = TRUE
				if (!is.null(n)) {
					private$assert_feasible_block_sizes(as.integer(n))
				}
			},

		#' @description
		#' Draw treatment assignments according to the optimal-blocks design.
		#' @param r Number of assignment vectors to draw.
		#' @return A numeric matrix with one assignment vector per column.
		draw_ws_according_to_design = function(r = 100){
			self$assert_all_subjects_arrived()
			block_ids = private$get_or_compute_block_ids()
			w_mat = replicate(r, as.numeric(as.character(randomizr::block_ra(blocks = block_ids, prob = private$prob_T))))
			storage.mode(w_mat) = "numeric"
			w_mat
		}
	),

	private = list(
		K = NULL,
		dist_spec = NULL,
		block_ids = NULL,
		distance_matrix = NULL,

		assert_feasible_block_sizes = function(n){
			if (private$K > n) {
				stop(
					"Cannot partition ", n, " subjects into ", private$K,
					" roughly equal blocks. With the floor(n / K) / ceiling(n / K) size rule, ",
					"K must be less than or equal to n."
				)
			}
			invisible(NULL)
		},

		get_or_compute_block_ids = function(){
			if (!is.null(private$block_ids)) {
				return(private$block_ids)
			}

			n = self$get_n()
			private$assert_feasible_block_sizes(n)

			if (is.null(private$X)) {
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
			}
			X = private$X[seq_len(n), , drop = FALSE]
			if (ncol(X) == 0L) {
				private$block_ids = factor(rep(seq_len(private$K), length.out = n))
				return(private$block_ids)
			}

			D = private$get_or_compute_distance_matrix(X)
			private$block_ids = private$solve_optimal_blocks(D)
			private$block_ids
		},

			get_or_compute_distance_matrix = function(X){
				if (!is.null(private$distance_matrix)) {
					return(private$distance_matrix)
				}

				if (is.function(private$dist_spec)) {
					D = optimal_blocks_distance_matrix_cpp(
						X = X,
						dist_code = 0L,
						dist_fn = private$dist_spec
					)
				} else {
					D = switch(private$dist_spec,
						euclidean = optimal_blocks_distance_matrix_cpp(X, dist_code = 1L),
						sum_abs_diff = optimal_blocks_distance_matrix_cpp(X, dist_code = 2L),
						mahal = optimal_blocks_distance_matrix_cpp(X, dist_code = 3L)
					)
				}

			storage.mode(D) = "double"
			private$distance_matrix = D
			D
		},

		solve_optimal_blocks = function(D){
			n = nrow(D)
			K = private$K
			lower_size = floor(n / K)
			upper_size = ceiling(n / K)

			model = ompr::MIPModel()
			model = ompr::add_variable(model, x[i, k], i = 1:n, k = 1:K, type = "binary")
			model = ompr::add_variable(model, z[i, j, k], i = 1:n, j = 1:n, k = 1:K, type = "binary", i < j)
			model = ompr::set_objective(
				model,
				ompr::sum_expr(D[i, j] * z[i, j, k], i = 1:n, j = 1:n, k = 1:K, i < j),
				"min"
			)

			model = ompr::add_constraint(model, ompr::sum_expr(x[i, k], k = 1:K) == 1, i = 1:n)
			model = ompr::add_constraint(model, ompr::sum_expr(x[i, k], i = 1:n) >= lower_size, k = 1:K)
			model = ompr::add_constraint(model, ompr::sum_expr(x[i, k], i = 1:n) <= upper_size, k = 1:K)
			model = ompr::add_constraint(model, x[k, k] == 1, k = 1:K)
			model = ompr::add_constraint(model, z[i, j, k] <= x[i, k], i = 1:n, j = 1:n, k = 1:K, i < j)
			model = ompr::add_constraint(model, z[i, j, k] <= x[j, k], i = 1:n, j = 1:n, k = 1:K, i < j)
			model = ompr::add_constraint(model, z[i, j, k] >= x[i, k] + x[j, k] - 1, i = 1:n, j = 1:n, k = 1:K, i < j)

			result = ompr::solve_model(model, ompr.roi::with_ROI(solver = "glpk", verbose = FALSE))
			solution = ompr::get_solution(result, x[i, k])
			if (nrow(solution) == 0L) {
				stop("ompr failed to produce a block assignment solution.")
			}
			solution = solution[solution$value > 0.5, , drop = FALSE]
			if (nrow(solution) != n) {
				stop("ompr returned an incomplete block assignment.")
			}
			labels = integer(n)
			labels[solution$i] = solution$k
			factor(labels, levels = seq_len(K))
		}
	)
)
