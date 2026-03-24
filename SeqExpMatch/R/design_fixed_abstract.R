# A Fixed Design
#
# @description
# An abstract R6 Class encapsulating the data and functionality for a fixed experimental design.
# This class takes care of data storage, response handling, and whole-experiment randomization.
#
# @keywords internal
#' @export
FixedDesign = R6::R6Class("FixedDesign",
	public = list(
		#
		# @description
		# Initialize a fixed experimental design
		#
		# @param response_type 	The data type of response values which must be one of the following:
		# 							"continuous" (the default),
		# 							"incidence",
		# 							"proportion",
		# 							"count",
		# 							"survival",
		# 							"ordinal".
		# @param prob_T	The probability of the treatment assignment. This defaults to \code{0.5}.
		# @param include_is_missing_as_a_new_feature	If missing data is present in a variable, should we include another dummy variable for its
		# 												missingness in addition to imputing its value?
		# @param n			The sample size (if fixed). Default is \code{NULL} for not fixed.
		# @param num_cores	The number of CPU cores to use.
		# @param verbose	A flag indicating whether messages should be displayed. Default is \code{FALSE}.
		#
		# @return 			A new `FixedDesign` object
		#
		initialize = function(
				response_type = "continuous",
				prob_T = 0.5,
				include_is_missing_as_a_new_feature = TRUE,
				n = NULL,
				num_cores = 1,
				verbose = FALSE
			) {
			assertChoice(response_type, c("continuous", "incidence", "proportion", "count", "survival", "ordinal"))
			assertNumeric(prob_T, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			assertFlag(include_is_missing_as_a_new_feature)
			assertFlag(verbose)
			assertCount(n, null.ok = TRUE)
			assertCount(num_cores, positive = TRUE)

			if (is.null(n)){
				private$fixed_sample = FALSE
			} else {
				n = as.integer(n)
				private$n = n
				private$fixed_sample = TRUE
			}

			private$prob_T = prob_T
			private$response_type = response_type
			private$include_is_missing_as_a_new_feature = include_is_missing_as_a_new_feature
			private$num_cores = num_cores
			set_package_threads(num_cores)
			private$verbose = verbose

			if (private$fixed_sample){
				private$y = 	rep(NA_real_, n)
				private$w = 	rep(NA_real_, n)
				private$dead =  rep(NA_real_, n)
			}

			if (private$verbose){
				cat(paste0("Initialized a ",
				class(self)[1],
				" experiment with response type ",
				response_type,
				" and ",
				ifelse(private$fixed_sample, "fixed sample", "not fixed sample"),
				 ".\n"))
			}
		},

		# @description
		# Add subject-specific measurements
		#
		add_subject = function(x_new){
			assertClass(x_new, "data.frame")
			x_new = as.data.table(x_new)
			if (nrow(x_new) != 1) stop("You can only add one subject at a time.")
			
			private$Xraw = rbindlist(list(private$Xraw, x_new), fill = TRUE)
			private$t = private$t + 1L

			if (private$t > (ncol(private$Xraw) + 2) & private$uses_covariates){
				private$covariate_impute_if_necessary_and_then_create_model_matrix()
			}
		},

		randomize = function(){
			self$redraw_w_according_to_design()
		},

		supports_resampling = function(){
			class(self)[1] != "FixedDesign"
		},

		redraw_w_according_to_design = function(){
			private$w[1:self$get_n()] = self$draw_ws_according_to_design(1)[, 1]
		},

		draw_ws_according_to_design = function(r = 100){
			stop("Must be implemented by subclass.")
		},

		add_subject_response = function(t, y, dead = 1) {
			if (t > private$t) stop("Subject not yet arrived.")
			if (private$fixed_sample | t <= length(private$y)){
				private$y[t] = y
				private$dead[t] = dead
			} else {
				private$y = c(private$y, y)
				private$dead = c(private$dead, dead)
			}
		},

		add_all_subject_responses = function(ys, deads = NULL) {
			if (is.null(deads)) deads = rep(1, private$t)
			for (t in 1 : private$t) self$add_subject_response(t, ys[t], deads[t])
		},

		add_all_subject_assignments = function(w) {
			private$w = w
		},

		assert_experiment_completed = function(){
			if (private$fixed_sample & private$t < private$n) stop("Experiment incomplete: sample size not met.")
			if (sum(!is.na(private$y)) < private$t) stop("Experiment incomplete: responses missing.")
		},

		any_censoring = function(){ sum(private$dead) < length(private$dead) },

		# Getters
		get_t = function(){ private$t },
		get_X_raw = function(){ private$Xraw },
		get_X_imp = function(){ private$Ximp },
		get_X = function() { private$X },
		get_y = function(){ private$y },
		get_w = function(){ private$w },
		get_n = function(){ ifelse(private$fixed_sample, private$n, private$t) },
		get_dead = function(){ private$dead },
		get_prob_T = function(){ private$prob_T },
		get_response_type = function(){ private$response_type },
		is_fixed_sample_size = function(){ private$fixed_sample },

		duplicate = function(verbose = FALSE){
			d = self$clone()
			d$.__enclos_env__$private$verbose = verbose
			d
		},

		assert_even_allocation = function(){
			if (private$prob_T != 0.5) stop("Requires prob_T = 0.5")
		},

		assert_fixed_sample = function(){
			if (!private$fixed_sample) stop("Requires fixed sample size n")
		}
	),

	private = list(
		t = 0L,
		n = NULL,
		Xraw = data.table(),
		Ximp = data.table(),
		X = NULL,
		w = numeric(),
		y = numeric(),
		dead = numeric(),
		m = NULL,
		prob_T = NULL,
		response_type = NULL,
		fixed_sample = NULL,
		num_cores = NULL,
		include_is_missing_as_a_new_feature = NULL,
		verbose = NULL,
		uses_covariates = FALSE,

		covariate_impute_if_necessary_and_then_create_model_matrix = function(){
			private$Ximp = copy(private$Xraw)

			column_has_missingness = columns_have_missingness_cpp(private$Xraw)
			if (any(column_has_missingness)){
				if (private$include_is_missing_as_a_new_feature){
					missing_cols_idx = which(column_has_missingness)
					if (length(missing_cols_idx) > 0){
						missingness_indicators = create_missingness_indicators_cpp(private$Ximp, missing_cols_idx)
						for (col_name in names(missingness_indicators)) {
							private$Ximp[[col_name]] = missingness_indicators[[col_name]]
						}
					}
				}

				col_types = get_column_types_cpp(private$Ximp)
				idx_cols_to_convert_to_factor = which(col_types == "character")
				private$Ximp[, (idx_cols_to_convert_to_factor) := lapply(.SD, as.factor), .SDcols = idx_cols_to_convert_to_factor]

				private$Ximp = tryCatch({
										if (any(!is.na(private$y))){
											suppressWarnings(missRanger(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), verbose = FALSE, num.threads = private$num_cores)[, 1 : ncol(private$Ximp)])
										} else {
											suppressWarnings(missRanger(private$Ximp, verbose = FALSE, num.threads = private$num_cores))
										}
									}, error = function(e){
										if (any(!is.na(private$y))){
											suppressWarnings(missForest(cbind(private$Ximp, private$y[1 : nrow(private$Ximp)]), num.threads = private$num_cores)$ximp[, 1 : ncol(private$Ximp)])
										} else {
											suppressWarnings(missForest(private$Ximp, num.threads = private$num_cores)$ximp)
										}
									}
								)
			}

			num_unique_values_per_column = count_unique_values_cpp(private$Ximp)
			private$Ximp = private$Ximp[, .SD, .SDcols = which(num_unique_values_per_column > 1)]

			if (ncol(private$Ximp) > 0){
				private$X = model.matrix(~ ., data = private$Ximp)[, -1, drop = FALSE]
				if (is.character(private$X)){
					stop("model.matrix returned a character matrix - this should not happen.")
				}
				private$X = drop_linearly_dependent_cols(private$X)$M
			} else {
				private$X = matrix(NA, nrow = nrow(private$Xraw), ncol = 0)
			}
		},

		compute_all_subject_data = function(){
			i_present_y = which(!is.na(private$y))
			compute_all_subject_data_cpp(as.matrix(private$X[1:private$t, , drop = FALSE]), private$t, as.integer(i_present_y))
		},

		has_private_method = function(method_name) {
			exists(method_name, envir = self$.__enclos_env__$private, inherits = FALSE)
		},

		assign_wt_Bernoulli = function(){
			rbinom(1, 1, private$prob_T)
		}
	)
)
