#' A Sequential Design
#'
#' @description
#' An R6 Class encapsulating the data and functionality for a sequential experimental design.
#' This class takes care of data initialization and sequential assignments. The class object
#' should be saved securely after each assignment e.g. on an encrypted cloud server.
#'
#' @export
SeqDesignKK14 = R6::R6Class("SeqDesignKK14",
	inherit = SeqDesign,
	public = list(
		#'
		#' @description
		#' Initialize a matching-on-the-fly sequential experimental design which matches based on
		#' Kapelner and Krieger (2014) or Morrison and Owen (2025)
		#'
		#' @param	response_type 	The data type of response values.
		#' @param	prob_T	The probability of the treatment assignment.
		#' @param include_is_missing_as_a_new_feature     Flag for missingness indicators.
		#' @param	n			The sample size.
		#' @param num_cores The number of CPU cores.
		#' @param verbose A flag for verbosity.
		#' @param lambda   The quantile cutoff.
		#' @param t_0_pct  The percentage where matching begins.
		#' @param morrison        Flag for Morrison formula.
		#' @param p                       The number of covariate features.
		#'
		#' @return	A new `SeqDesignKK14` object
		#'
		initialize = function(
			response_type = "continuous",
			prob_T = 0.5,
			include_is_missing_as_a_new_feature = TRUE,
			n = NULL,
			num_cores = 1,
			verbose = FALSE,
			lambda = NULL,
			t_0_pct = NULL,
			morrison = FALSE,
			p = NULL
		){
			super$initialize(response_type, prob_T, include_is_missing_as_a_new_feature, n, num_cores, verbose)
			self$assert_even_allocation()
			private$assert_KK_and_morrison_parameters_correct(lambda, t_0_pct, morrison, p)
			private$morrison = morrison
			private$uses_covariates = TRUE
			if (morrison){
				private$p = p
				if (self$is_fixed_sample_size()){
	              	private$compute_lambda = function(){private$n^(-1 / (2 * private$p))}
	            } else {
	              	private$compute_lambda = function(){private$t^(-1 / (2 * private$p))}
	            }
			} else {
				private$m = array(NA, n)
				private$lambda = ifelse(is.null(lambda), 0.1, lambda) #10% is the default
				private$compute_lambda = function(){private$lambda}
				private$t_0 = round(ifelse(is.null(t_0_pct), 0.35, t_0_pct) * n) #35% is the default
			}
		},

		#' @description
		#' This returns a list with useful matching statistics.
		matching_statistics = function(){
			if (private$t == 0){
				stop("The experiment has not begun yet")
			}
			num_subjects_matched = sum(private$m != 0, na.rm = TRUE)
			num_subjects_remaining_in_reservoir = private$t - num_subjects_matched
			list(
				num_matches = length(unique(private$m[private$m != 0])) ,
				prop_subjects_matched = num_subjects_matched / private$t,
				num_subjects_remaining_in_reservoir = num_subjects_remaining_in_reservoir,
				prop_subjects_remaining_in_reservoir = num_subjects_remaining_in_reservoir / private$t
			)
		},

		#' @description
		#' Returns the grouping vector `m`
		get_m = function(){ private$m },

		draw_ws_according_to_design = function(r = 100){
			generate_permutations_kk_cpp(as.integer(private$m), as.integer(r), as.numeric(private$prob_T))
		},

		assign_wt = function(){
			wt = 	if (private$too_early_to_match()){
						private$m[private$t] = 0
						rbinom(1, 1, private$prob_T)
					} else {
						all_subject_data = private$compute_all_subject_data()
						S_xs_inv = solve(var(all_subject_data$X_prev) + diag(.Machine$double.eps, all_subject_data$rank_prev), tol = .Machine$double.xmin)
						reservoir_indices = which(private$m == 0)
						sqd_distances_times_two = compute_proportional_mahal_distances_cpp(
						    all_subject_data$xt_prev,
						    all_subject_data$X_prev,
						    reservoir_indices,
						    S_xs_inv
						)
						F_crit =  qf(private$compute_lambda(), all_subject_data$rank_prev, private$t - all_subject_data$rank_prev)
						n = self$get_n()
						T_cutoff_sq = all_subject_data$rank_prev * (n - 1) / (n - all_subject_data$rank_prev) * F_crit
						min_sqd_dist_index = which(sqd_distances_times_two == min(sqd_distances_times_two))[1]

						if (sqd_distances_times_two[min_sqd_dist_index] < T_cutoff_sq){
							match_num = max(private$m, na.rm = TRUE) + 1
							private$m[reservoir_indices[min_sqd_dist_index]] = match_num
							private$m[private$t] = match_num
							1 - private$w[reservoir_indices[min_sqd_dist_index]]
						} else {
							private$m[private$t] = 0
							rbinom(1, 1, private$prob_T)
						}
					}
			if (is.na(private$m[private$t])){
				stop("no match data recorded")
			}
			wt
		}
	),
	private = list(
		m = NULL,
		uses_covariates = TRUE,
		morrison = NULL,
		t_0 = 0,
		lambda = NULL,
		p = NULL,
		compute_lambda = NULL,

		duplicate = function(){
			d = super$duplicate()
			d
		},

		assert_KK_and_morrison_parameters_correct = function(lambda, t_0_pct, morrison, p){
			if (morrison){
				assert_count(p)
			} else {
				self$assert_fixed_sample()
				if (!is.null(lambda)){
					assertNumeric(lambda, lower = 0, upper = 1)
				}
				if (!is.null(t_0_pct)){
					assertNumeric(t_0_pct, lower = .Machine$double.eps, upper = 1)
				}
			}
		},

		too_early_to_match = function(){
			sum(private$m == 0, na.rm = TRUE) == 0 | private$t <= (ncol(private$Xraw) + 2) |
							(!private$morrison & private$t <= private$t_0)
		},

		redraw_w_according_to_design = function(){
			private$w[1:private$t] = redraw_w_kk14_cpp(private$m, private$w)
		}
	)
)
